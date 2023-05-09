function [entr_pt_pain_tbl, entr_pt_pain_stats, inflation, pts_MI_in_bits, pts_MI_perm_p ]...
    ...
    = pain_entr_MI(...
    ...
    cfg, pts_pain_space_struct)

% cfg.pts             = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};
% 
% 
% cfg.tbl_pain_lbls   = {'mayoNRS','painVAS','unpleasantVAS', 'MPQtotal', ...
%                        'PBD_sum', 'PBD_coverage', 'PBD_mean'};
% 
% cfg.nice_plt_lbls   = {'NRS int.','VAS int.', 'VAS unpl.', 'MPQ total',...
%                        'PBD sum', 'PBD cov.', 'PBD mean'};
% 
% cfg.save_dir        = s0_results_dir;
% 
% 


pts                 = cfg.pts;
tbl_pain_lbls       = cfg.tbl_pain_lbls;
nice_plt_lbls       = cfg.nice_plt_lbls;

pt_lbls             = cfg.pt_lbls;


%% calc N bins
% ignoring NRS (for now) as it is a 11-point scale (i.e., max of 11 states)
lbls_for_N_bin_calc        = tbl_pain_lbls(~contains(tbl_pain_lbls, 'NRS'));

n_bins          = zeros(length(lbls_for_N_bin_calc), length(pts));
save_dir        = [cfg.save_dir, 'information_theory_analysis/'];

for i = 1:length(pts)
    redcap           = pts_pain_space_struct.(['stage0', pts{i}]).redcap_kept;
    data             = redcap{:, lbls_for_N_bin_calc};

    %orient data
    [n_rows,n_cols]=size(data);
    if n_cols>n_rows 
        data=data';
        [n_rows,~]=size(data);
    end

    %Freedman-Diaconis rule for histogram bin estimation (Freedman-Diaconis, 1981)
    n_bins(:, i) = ceil(range(data) ./ (2.*iqr(data).*n_rows.^(-1/3)));

end
% format explictly as tbl
%n_bins_tbl = array2table(n_bins, "VariableNames", pts, "RowNames", pain_lbls);

n_bins_used = ceil(mean(n_bins,"all"));

fprintf('%.0f bins via Freedman-Diaconis rule (across ALL pts and pain metrics (excluding NRS))\n',...
         n_bins_used);

%% w/ N bins decided -> calculate entropy per pt and pain metrics
entr_pt_pain = zeros(length(tbl_pain_lbls), length(pts));

for i = 1:length(pts)
    redcap   = pts_pain_space_struct.(['stage0', pts{i}]).redcap_kept;
    data     = redcap(:, tbl_pain_lbls);

    % as entropy changes w/ N, only include where ALL metrics were collected
    data     = data{all(~isnan(data.Variables), 2), :};


    [entr_pt_pain(:, i), inflation, ~, ~] = entropy1D(data, n_bins_used);

    figure('Units', 'Inches', 'Position', [0, 0, 4, 3.25]);
    sgtitle(pt_lbls{i},'Fontsize', 10);

    for j=1:length(nice_plt_lbls)
        subplot(4, 2, j)  
       
        histogram(data(:,j), n_bins_used, 'Normalization','probability');

        title(nice_plt_lbls{j},'FontWeight', 'normal', 'Interpreter','none'); 

        ylim([0, .6])
        set(gca, 'FontName', 'Arial','Fontsize',8);
        grid on; grid minor
    end

    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    
    exportgraphics(gcf, ...
        sprintf('%sintermediate_steps/histograms_alone/stage0%s_hist_Freed_Diac_rule.png',...
        save_dir, pts{i}));
end

entr_pt_pain_tbl = array2table(entr_pt_pain', ...
                   "RowNames", pts, "VariableNames",tbl_pain_lbls);

entr_pt_pain_tbl{'grp_mean', :} = round(mean(entr_pt_pain_tbl{pts,:}),2);
entr_pt_pain_tbl{'grp_std',:}   = round(std(entr_pt_pain_tbl{pts,:}),2);

%% visualize entropy as boxplots per pain metric w/ colored lines indicating pt
figure('Units', 'Inches', 'Position', [3, 5, 10, 5]);

[xpos, ypos, ~, ~, ~]=UnivarScatter(entr_pt_pain_tbl{pts,tbl_pain_lbls},...
    'Compression', 100,... distance btwn boxes
    'Width', 1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 8);
hold on;

xticklabels(nice_plt_lbls);  ylabel('Entropy (bits)');

h_mean = plot(xpos.',ypos.','Linewidth',1,'Color','k');

colors = brewermap(length(pts),'Dark2');

set(h_mean, {'color'}, num2cell(colors,2));

L = legend(h_mean,pt_lbls,'Box','off', 'Location','northoutside', ...
    'Orientation', 'horizontal');

L.ItemTokenSize(1) = 40;

set(gca,'Fontsize',14,'TickLabelInterpreter','none');

exportgraphics(gcf, [save_dir, '/intermediate_steps/entropy_across_pts_pain.png']);

%% one-way repeated measures ANOVA of entropy per pain metric across pts
% w/ Tukey's t-test for individual comparisons
pts_by_pain_entropy = entr_pt_pain_tbl(:, cfg.tbl_pain_lbls);


[~,~,stats] = anova1(pts_by_pain_entropy.Variables);


rm = fitrm(pts_by_pain_entropy, ...
           sprintf('%s-%s~1', tbl_pain_lbls{1}, tbl_pain_lbls{2}));

ranova_tbl = ranova(rm);

if ranova_tbl.pValue(1) <=0.05
    c = multcompare(stats);

    c = c( c(:,end) <= 0.05, :);

    multcomp_tbl = table;

     for i=1:length(cfg.tbl_pain_lbls)

         multcomp_tbl.Var1(c(:,1) == i) = cfg.nice_plt_lbls(i);
         multcomp_tbl.Var2(c(:,2) == i) = cfg.nice_plt_lbls(i);

     end

     multcomp_tbl.lower_CI = round(c(:,4),3);
     multcomp_tbl.upper_CI = round(c(:,5),3);

     multcomp_tbl.p_val = round(c(:,6),3);

end

entr_pt_pain_stats.ranova_tbl       = ranova_tbl;
entr_pt_pain_stats.multcomp_tbl     = multcomp_tbl;

%% btwn all pain metrics, calculate the mutual information
% leveraging Neuroscience Information Theory Matlab Toolbox 
% (Timme, and Lapish, 2018)

pts_MI_in_bits    = nan(length(pts), length(tbl_pain_lbls), length(tbl_pain_lbls)); 
pts_MI_perm_p     = pts_MI_in_bits;

pts_Pxy           = cell(length(pts),1);
pts_Px_given_y    = pts_Pxy;
pts_bin_edges     = nan(length(pts), length(tbl_pain_lbls), n_bins_used +1);
pts_Px            = nan(length(pts), length(tbl_pain_lbls), n_bins_used);

pts_i_x_states    = cell(length(pts),1);

for i = 1:length(pts)
    redcap           = pts_pain_space_struct.(['stage0', pts{i}]).redcap_kept;
    data             = redcap(:, tbl_pain_lbls);

    data        = data(all(~isnan(data.Variables), 2), :);

    MI          = nan(length(tbl_pain_lbls), length(tbl_pain_lbls));      
    MI_perm_p   = MI;
    RHO         = MI;
    PVAL        = MI;

    Pxy         = cell(size(MI));
    Px_given_y  = Pxy;

    bin_edges     = nan(length(tbl_pain_lbls), n_bins_used +1);

    Px            =  nan(length(tbl_pain_lbls), n_bins_used);

    i_x_states    = cell(length(tbl_pain_lbls),1);

    N_comparisons = 0;


    for j =1:length(tbl_pain_lbls)-1
        for h = j+1:length(tbl_pain_lbls)
           
            [MI(h, j), MI_perm_p(h,j),...
             Pxy{h,j}, Px_given_y{h,j}, ...
             bin_edges(j ,:), Px(j,:), i_x_states{j}] ...
                ...
                = quickMI(...
                ...
                data{:, tbl_pain_lbls{j}},...
                data{:, tbl_pain_lbls{h}},...
                'nbins', n_bins_used);

            N_comparisons = N_comparisons + 1;

        [RHO(i,j,h), PVAL(i,j,h)] = corr(data{:, tbl_pain_lbls{j}},data{:, tbl_pain_lbls{h}},'Type','Spearman', 'rows','pairwise');
    

        end
    end

    [~, ~, ~, ~, bin_edges(length(tbl_pain_lbls) ,:), Px(length(tbl_pain_lbls),:), i_x_states{length(tbl_pain_lbls)}] ...
    ...
    = quickMI(...
    ...
    data{:, tbl_pain_lbls{length(tbl_pain_lbls)}},  data{:, tbl_pain_lbls{1}},  'nbins', n_bins_used);


    % from lower triangluar matrix to symmetric matrix for readablity
    MI          = triu(MI',1) + triu(MI')';
    MI_perm_p   = triu(MI_perm_p',1) + triu(MI_perm_p')';

    % format MI and MI permuatation p-value explicitly as tbls
%     MI_tbl         = array2table(MI, "VariableNames", pain_lbls, ...
%                                      "RowNames",pain_lbls);
% 
%     MI_perm_p_tbl  = array2table(MI_perm_p, "VariableNames", pain_lbls,...
%                                             "RowNames",pain_lbls);

    pts_MI_in_bits(i, :, :)         = MI;
    pts_MI_perm_p(i, :, :)  = MI_perm_p;

    pts_Pxy{i}       = Pxy;

    pts_Px_given_y{i} = Px_given_y;

    pts_bin_edges(i, :, :) = bin_edges;

    pts_Px(i, :, :)     = Px;

    pts_i_x_states{i} = i_x_states;
    
    clear MI MI_perm_p Pxy Px_given_y bin_edges Px
end
%% Spearman's correlation--01/27/23 RBL--not fully implemented
%{
preset_colormap(@brewermap, "reds")
figure('Units', 'Inches', 'Position', [0, 0, 6, 9])

for i= 1 :length(pts)

   
    rho   = squeeze(RHO(i, :, :));
    pval  = squeeze(PVAL(i,:,:));

    % Bonferroni correction on the subject-level
    corr_p = round(0.05/N_comparisons, 5);

    rho(pval > corr_p) = nan;

    rho      = triu(rho) + tril(nan(size(rho)));

    plt_rho = zeros(size(rho ,1)+1, size(rho ,2)+1);
    plt_rho(1:size(rho ,1), 1:size(rho ,2)) = rho';

    mat_pos = 1:length(tbl_pain_lbls)+1;


    subplot(3,2,i);
    s = pcolor(mat_pos, flip(mat_pos), plt_rho); hold on

    clim([0,1])

    s.EdgeColor = [0.7 0.7 0.7];     s.LineWidth = 1.5;

    colormap(flip(preset_colormap))
    cb_hdl = colorbar; 



end
%}
%% conditional probability surfaces amoungst pain metrics
close all
prop_mat   = pts_Px_given_y;
sg_tit_txt = 'conditional probability [0,1] between pain metrics';


preset_colormap(@brewermap, "YlGnBu")

cb_max = zeros(length(pts), 1);
for i=1 :length(pts)
    cell_mat  = prop_mat{i};
    max_Pxy   = cellfun(@(x) max(x, [], 'all'), cell_mat, 'UniformOutput', false);
    i_emp     = cellfun(@(x) isempty(x), max_Pxy);
    cb_max(i) = max(cell2mat( max_Pxy(~i_emp)));
end

[n_rows, n_cols] = size(cell_mat);
hist_max = round(1.1*max(pts_Px,[], [2,3]),1);

for i= 1 : length(pts)
    fig_prob_mat = figure('Units', 'Inches', 'Position', [0, 0, 5.5, 5]);
    clim_var =[0, 1];

    for i_row=1  : n_rows
        for i_col = 1 : n_cols

            % rc is grid position of subplot
            rc = [i_row, i_col];
            pos = (rc(1) * n_cols) - (n_cols - rc(2));
            
            subplot(n_rows,n_cols,pos);

            if i_row == i_col

                % stem plot such aligned w/ middle of prob cell
                stem(1.5:n_bins_used+1.5, [squeeze(pts_Px(i, i_col, :)); 0] , ...
                    'Color','k','Marker', 'none', 'LineWidth', 1.75)

                plt_prob =squeeze(pts_Px(i, :, :));

                ylim([0, hist_max(i)]);  yticklabels([]); 
                
                xlim([1, n_bins_used+1])

                bin_edges = squeeze(pts_bin_edges(i, i_col, :));

                bin_edges = bin_edges(~isinf(bin_edges) & ~isnan(bin_edges));

                bin_edge_str = [num2str(min(bin_edges),'%0.0f');...
                               repmat({''}, n_bins_used-1, 1);...
                               num2str(max(bin_edges),'%0.0f')]';

                xticks(1: n_bins_used+1)

                xticklabels(bin_edge_str)

                set(gca, 'XTickLabelRotation',0, 'FontSize', 8,'TickLength', [0 0]);

            elseif i_row > i_col
                tmp_mat = prop_mat{i}{i_row, i_col}';

                plt_mat = zeros(size(tmp_mat,1)+1, size(tmp_mat,2)+1);

                plt_mat(1:end-1, 1:end-1) = tmp_mat;

                plt_mat(isnan(plt_mat)) = 0;

                plt_mat_ax = pcolor(plt_mat);     hold on;    clim(clim_var); 

                %plt_mat_ax.EdgeColor = [0.85, 0.85, 0.85];
                plt_mat_ax.LineStyle = 'none';

                yticklabels([]);  
                xticklabels([]);   
            end

            if i_row == n_rows
                x_lbl = xlabel(nice_plt_lbls(i_col), 'Interpreter','none'); hold on
                set(gca, 'Fontsize', 8)

                 if i_row ~= i_col
                    x_lbl.Position(2) = -4.8;
                    set(gca, 'Fontsize', 8)

                 else
                    x_lbl.Position(2) = x_lbl.Position(2) *1.1;

                 end
            end
            if i_col == 1
                y_lbl = ylabel(nice_plt_lbls(i_row), 'Interpreter','none', 'Rotation', 0); hold on

                y_lbl.Units       = "normalized";
                y_lbl.Position(1) = -0.7;
                set(gca, 'Fontsize', 8)
        
            end 
        end
    end

    % delete upper triangle
    haxes    = reshape(fig_prob_mat.Children, n_rows, n_rows);
    i_delete = triu(ones(n_rows, n_cols),1);
    
    delete(haxes(i_delete == true))
    
    % use unidirectional colormap
    colormap(preset_colormap)
    cb_hdl = colorbar(haxes(2,1));   cb_hdl.Limits = clim_var;     
    cb_hdl.Position = [.925, .115, .02, .82];
    
    cb_hdl.FontSize = 8;             cb_hdl.Ticks  = 0:.5:1;
 
    sgtitle([pt_lbls{i},' ',sg_tit_txt], 'FontSize', 12);
 

    exportgraphics(gcf, [save_dir,'intermediate_steps/stage0',pts{i}, '_CondProb_matrices.png'])
end


%% visualize pair-wise MI as matrix
close all
%
%%%
max_poss_infor  = log2(n_bins_used);
pts_MI_in_per   = 100*pts_MI_in_bits ./ max_poss_infor;
%%%
%
%
preset_colormap(@brewermap, "YlOrRd");

cb_lim = [min(pts_MI_in_per,[], 'all'), max(pts_MI_in_per,[], 'all')];

pts_plt_MI = pts_MI_in_per;

for i= 1 :length(pts)

    figure('Units', 'Inches', 'Position', [0, 0, 4, 3])

    MI      = squeeze(pts_MI_in_per(i, :, :));
    perm_MI = squeeze(pts_MI_perm_p(i,:,:));

    % Bonferroni correction on the subject-level
%     corr_p = round(0.05/N_comparisons, 5);
%     MI(perm_MI >= corr_p ) = nan;

    MI(perm_MI > 0.05 ) = 0;

    MI      = tril(MI) + triu(nan(size(MI)));

    plt_MI = nan(size(MI,1), size(MI,2));
    plt_MI(1:size(MI,1)-1, 1:size(MI,2)) = MI(2:end,:);

    mat_pos = 0:length(tbl_pain_lbls)-1;

    s = pcolor(mat_pos, flip(mat_pos), plt_MI);

    s.EdgeColor = [1 1 1];      s.LineWidth = 1.5;
    
    clim(cb_lim);
     
    lbl_pos = (0:length(tbl_pain_lbls))+0.5;

    xticks(lbl_pos);    xticklabels(nice_plt_lbls);       
    yticks(lbl_pos);    yticklabels(flip(nice_plt_lbls));

    set(gca,'Fontsize',9,'TickLabelInterpreter','none','Color','w',...
        'TickLength', [0,0]);


    t=title(pt_lbls{i}, 'FontSize', 12);

    pts_plt_MI(i, :, :) = MI;

    colormap([[.85 .85 .85];preset_colormap])
    set(gcf, 'InvertHardCopy', 'off')
    
    
    cb_hdl = colorbar;          cb_hdl.Limits = cb_lim; 
    
    cb_hdl.Label.String         = "Percent of max possible MI";
    cb_hdl.Label.Rotation       = 270;
    cb_hdl.Label.Position(1)    = 3.3;
    cb_hdl.FontSize             = 10;
    cb_hdl.Label.FontSize       = 12;

    exportgraphics(gcf, [save_dir,'stage0',pts{i}, '_MI_pain_matrix.png'])

end

%% visualize the group-level MI mean and prevelance of sig. subject-level comparisons
grp_fig =  figure('Units', 'Inches', 'Position', [0, 0, 4, 3]);


MI   = squeeze(mean(pts_plt_MI,'omitnan'));
prop = flip(squeeze(sum(pts_plt_MI ~=0,1)));


plt_MI = nan(size(MI,1), size(MI,2));
plt_MI(1:size(MI,1)-1, 1:size(MI,2)) = MI(2:end,:);

mat_pos = 0:length(tbl_pain_lbls)-1;

s = pcolor(mat_pos, flip(mat_pos), plt_MI);         clim(cb_lim);

s.EdgeColor = [1 1 1];          s.LineWidth = 1.5;
 

lbl_pos = (0:length(tbl_pain_lbls))+0.5;

xticks(lbl_pos);    xticklabels(nice_plt_lbls);       
yticks(lbl_pos);    yticklabels(flip(nice_plt_lbls));


prop_MI = nan(size(prop,1), size(prop,2));
prop_MI(1:size(MI,1)-1, 2:size(MI,2)) = prop(1:end-1, 2:end);

set(gcf, 'InvertHardCopy', 'off')
set(gca,'Fontsize', 9,'TickLabelInterpreter','none',...
    'Color','w', 'TickLength', [0,0]);

colormap(preset_colormap);
cb_hdl = colorbar;          cb_hdl.Limits = cb_lim; 

cb_hdl.Label.String         = "Percent of max possible MI";
cb_hdl.Label.Rotation       = 270;
cb_hdl.Label.Position(1)    = 3.3;
cb_hdl.FontSize             = 10;
cb_hdl.Label.FontSize       = 12;

t = title("Group-level mean of sig. subjects");

t.FontSize    = 12;
%%%
 for i=1:length(grp_fig.CurrentAxes.XTick)-2

    for j=1:length(grp_fig.CurrentAxes.YTick)-1 -i

        text(.95*grp_fig.CurrentAxes.XTick(i),...
             grp_fig.CurrentAxes.YTick(j), ...
             [num2str(prop(j,i)), '/5'],...
             'FontSize', 8);
    end
end


exportgraphics(gcf, [save_dir,'stage0','GroupLevel_MI_pain_matrix.png']);

end