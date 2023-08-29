function [entropy_struct, pts_MI_in_bits, pts_MI_perm_p ]...
    ...
    = pain_entr_MI(...
    ...
    cfg, pts_pain_space_struct)

%% calc N bins
% ignoring NRS (for now) as it is a 11-point scale (i.e., max of 11 states)
lbls_for_N_bin_calc        = cfg.tbl_pain_lbls(~contains(cfg.tbl_pain_lbls, 'NRS'));

n_bins          = zeros(length(lbls_for_N_bin_calc), length(cfg.pts));
save_dir        = fullfile(cfg.save_dir, 'information_theory_analysis/');

if ~isfolder(save_dir);     mkdir(save_dir);    end

for i = 1:length(cfg.pts)
    redcap           = pts_pain_space_struct.(['stage0', cfg.pts{i}]).redcap_kept;

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
%n_bins_tbl = array2table(n_bins, "VariableNames", cfg.pts, "RowNames", pain_lbls);

n_bins_used = ceil(mean(n_bins,"all"));

fprintf('%.0f bins via Freedman-Diaconis rule (across ALL cfg.pts and pain metrics (excluding NRS))\n',...
         n_bins_used);

entropy_struct.N_bins = n_bins_used;
%% w/ N bins decided -> calculate entropy per pt and pain metrics

hist_dir = fullfile(save_dir, 'intermediate_steps/histograms_alone');

if ~isfolder(hist_dir);   mkdir(hist_dir);   end

entr_pt_pain = zeros(length(cfg.pts), length(cfg.tbl_pain_lbls));
N_surveys    = entr_pt_pain;
inflation    = entr_pt_pain;

for i_pt = 1:length(cfg.pts)
    redcap   = pts_pain_space_struct.(['stage0', cfg.pts{i_pt}]).redcap_kept;
    data     = redcap(:, cfg.tbl_pain_lbls);

    % as entropy changes w/ N, only include where ALL metrics were collected
    data     = data{all(~isnan(data.Variables), 2), :};

    [entr_pt_pain(i_pt,:), inflation(i_pt, :), ~, ~] = entropy1D(data, n_bins_used);

    figure('Units', 'Inches', 'Position', [0, 0, 4, 3.25]);
    sgtitle(cfg.pt_lbls{i_pt},'Fontsize', 10);

    for j_pain = 1:length(cfg.nice_plt_lbls)
        subplot(4, 2, j_pain)  
       
        histogram(data(:,j_pain), n_bins_used, 'Normalization','probability');

        N_surveys(i_pt, j_pain) = sum(~isnan(data(:,j_pain)));

        title(sprintf('%s, %g reports', cfg.nice_plt_lbls{j_pain}, N_surveys(i_pt, j_pain)),...
            'FontWeight', 'normal', 'Interpreter','none'); 

        ylim([0, .6])
        set(gca, 'FontName', 'Arial','Fontsize',8);
        grid on; grid minor
    end
    exportgraphics(gcf, ...
                fullfile(hist_dir,  ...
                    sprintf('stage0%s_hist_Freed_Diac_rule.png', cfg.pts{i_pt}))...
                    );
end

% inflation is based on N entries and bin size
entropy_struct.inflation = array2table(inflation, ...
                   "RowNames", cfg.pts, "VariableNames",cfg.tbl_pain_lbls);

entropy_struct.inflation{'grp_mean', :} = round(mean(entropy_struct.inflation{cfg.pts,:}),2);
entropy_struct.inflation{'grp_std',:}   = round(std(entropy_struct.inflation{cfg.pts,:}),2);


% entropy per patient's pain metric
entropy_struct.bits = array2table(entr_pt_pain, ...
                   "RowNames", cfg.pts, "VariableNames",cfg.tbl_pain_lbls);

entropy_struct.bits{'grp_mean', :} = round(mean(entropy_struct.bits{cfg.pts,:}),2);
entropy_struct.bits{'grp_std',:}   = round(std(entropy_struct.bits{cfg.pts,:}),2);

% Number of reports per patient's pain metric
entropy_struct.N_reports = array2table(N_surveys, ...
                   "RowNames", cfg.pts, "VariableNames",cfg.tbl_pain_lbls);

entropy_struct.N_reports{'grp_mean', :} = round(mean(entropy_struct.N_reports{cfg.pts,:}),2);
entropy_struct.N_reports{'grp_std',:}   = round(std(entropy_struct.N_reports{cfg.pts,:}),2);



%% visualize entropy as boxplots per pain metric w/ colored lines indicating pt
figure('Units', 'Inches', 'Position', [3, 5, 10, 5]);

[xpos, ypos, ~, ~, ~]=UnivarScatter(entropy_struct.bits{cfg.pts,cfg.tbl_pain_lbls},...
    'Compression', 100,... distance btwn boxes
    'Width', 1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 8);
hold on;

xticklabels(cfg.nice_plt_lbls);  ylabel('Entropy (bits)');

h_mean = plot(xpos.',ypos.','Linewidth',1,'Color','k');

% kept for generalized use of color rather than this submission
%colors = brewermap(length(cfg.pts),'Dark2');


colors = [68 86 148;...
         162 58 46;...
         1 102 94;...
         84 48 5;...
         157 60 219] ./ 255;

set(h_mean, {'color'}, num2cell(colors,2));

L = legend(h_mean,cfg.pt_lbls,'Box','off', 'Location','northoutside', ...
    'Orientation', 'horizontal');

L.ItemTokenSize(1) = 40;

set(gca,'Fontsize',14,'TickLabelInterpreter','none');

exportgraphics(gcf, [save_dir, '/intermediate_steps/entropy_across_pts_pain.png']);

%% one-way repeated measures ANOVA of entropy per pain metric across cfg.pts
% w/ Tukey's t-test for individual comparisons
tmp_tbl     = entropy_struct.bits(:, cfg.tbl_pain_lbls);
[~,~,stats] = anova1(tmp_tbl.Variables);

rm = fitrm(tmp_tbl, sprintf('%s-%s~1', cfg.tbl_pain_lbls{1}, cfg.tbl_pain_lbls{2}));

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

    entropy_struct.pt_pain_stats.ranova_tbl       = ranova_tbl;
    entropy_struct.pt_pain_stats.multcomp_tbl     = multcomp_tbl;

end


%% btwn all pain metrics, calculate the mutual information
% leveraging Neuroscience Information Theory Matlab Toolbox 
% (Timme, and Lapish, 2018)

pts_MI_in_bits    = nan(length(cfg.pts), length(cfg.tbl_pain_lbls), length(cfg.tbl_pain_lbls)); 
pts_MI_perm_p     = pts_MI_in_bits;

pts_Pxy           = cell(length(cfg.pts),1);
pts_Px_given_y    = pts_Pxy;
pts_bin_edges     = nan(length(cfg.pts), length(cfg.tbl_pain_lbls), n_bins_used +1);
pts_Px            = nan(length(cfg.pts), length(cfg.tbl_pain_lbls), n_bins_used);

pts_i_x_states    = cell(length(cfg.pts),1);

for i = 1:length(cfg.pts)
    redcap           = pts_pain_space_struct.(['stage0', cfg.pts{i}]).redcap_kept;
    data             = redcap(:, cfg.tbl_pain_lbls);

    %data        = data(all(~isnan(data.Variables), 2), :);

    MI.bits          = nan(length(cfg.tbl_pain_lbls), length(cfg.tbl_pain_lbls));      
    MI.perm_p   = MI.bits;
    RHO         = MI.bits;
    PVAL        = MI.bits;

    MI.Pxy         = cell(size(MI.bits));
    MI.Px_given_y  = MI.Pxy;

    MI.bin_edges     = nan(length(cfg.tbl_pain_lbls), n_bins_used +1);

    Px            =  nan(length(cfg.tbl_pain_lbls), n_bins_used);

    i_x_states    = cell(length(cfg.tbl_pain_lbls),1);

    N_comparisons = 0;

    for j_pain =1:length(cfg.tbl_pain_lbls)-1
        for h = j_pain+1:length(cfg.tbl_pain_lbls)

            j_pain_vec = data{:, cfg.tbl_pain_lbls{j_pain}};
            h_pain_vec = data{:, cfg.tbl_pain_lbls{h}};

            i_keep = find(~isnan(j_pain_vec) & ~isnan(h_pain_vec));
            
            [MI.bits(h, j_pain), MI.perm_p(h,j_pain),...
             MI.Pxy{h,j_pain}, MI.Px_given_y{h,j_pain}, ...
             MI.bin_edges(j_pain ,:), Px(j_pain,:), i_x_states{j_pain}] ...
                ...
                = quickMI(...
                ...
                j_pain_vec(i_keep), h_pain_vec(i_keep),...
                'nbins', n_bins_used);

            N_comparisons = N_comparisons + 1;

        [RHO(i,j_pain,h), PVAL(i,j_pain,h)] = corr(data{:, cfg.tbl_pain_lbls{j_pain}},data{:, cfg.tbl_pain_lbls{h}},'Type','Spearman', 'rows','pairwise');
    

        end
    end

     j_pain_vec = data{:, cfg.tbl_pain_lbls{length(cfg.tbl_pain_lbls)}};
     h_pain_vec = data{:, cfg.tbl_pain_lbls{1}};

     i_keep = find(~isnan(j_pain_vec) & ~isnan(h_pain_vec));

    [~, ~, ~, ~, MI.bin_edges(length(cfg.tbl_pain_lbls) ,:), Px(length(cfg.tbl_pain_lbls),:), i_x_states{length(cfg.tbl_pain_lbls)}] ...
    ...
    = quickMI(...
    ...
   j_pain_vec(i_keep),  h_pain_vec(i_keep),  'nbins', n_bins_used);


    % from lower triangluar matrix to symmetric matrix for readablity
    MI.bits          = triu(MI.bits',1) + triu(MI.bits')';
    MI.perm_p   = triu(MI.perm_p',1) + triu(MI.perm_p')';

    % format MI.bits and MI.bits permuatation p-value explicitly as tbls
%     MI_tbl         = array2table(MI.bits, "VariableNames", pain_lbls, ...
%                                      "RowNames",pain_lbls);
% 
%     MI_perm_p_tbl  = array2table(MI.perm_p, "VariableNames", pain_lbls,...
%                                             "RowNames",pain_lbls);

    pts_MI_in_bits(i, :, :)         = MI.bits;
    pts_MI_perm_p(i, :, :)  = MI.perm_p;

    pts_Pxy{i}       = MI.Pxy;

    pts_Px_given_y{i} = MI.Px_given_y;

    pts_bin_edges(i, :, :) = MI.bin_edges;

    pts_Px(i, :, :)     = Px;

    pts_i_x_states{i} = i_x_states;
    
    clear MI.bits MI.perm_p MI.MI.Pxy MI.Px_given_y MI.bin_edges Px
end
%% Spearman's correlation--01/27/23 RBL--not fully implemented
%{
preset_colormap(@brewermap, "reds")
figure('Units', 'Inches', 'Position', [0, 0, 6, 9])

for i= 1 :length(cfg.pts)

   
    rho   = squeeze(RHO(i, :, :));
    pval  = squeeze(PVAL(i,:,:));

    % Bonferroni correction on the subject-level
    corr_p = round(0.05/N_comparisons, 5);

    rho(pval > corr_p) = nan;

    rho      = triu(rho) + tril(nan(size(rho)));

    plt_rho = zeros(size(rho ,1)+1, size(rho ,2)+1);
    plt_rho(1:size(rho ,1), 1:size(rho ,2)) = rho';

    mat_pos = 1:length(cfg.tbl_pain_lbls)+1;


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

cb_max = zeros(length(cfg.pts), 1);
for i=1 :length(cfg.pts)
    cell_mat  = prop_mat{i};
    max_Pxy   = cellfun(@(x) max(x, [], 'all'), cell_mat, 'UniformOutput', false);
    i_emp     = cellfun(@(x) isempty(x), max_Pxy);
    cb_max(i) = max(cell2mat(max_Pxy(~i_emp)));
end

[n_rows, n_cols] = size(cell_mat);
hist_max = round(1.1*max(pts_Px,[], [2,3]),1);

for i= 1 : length(cfg.pts)
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

                MI.bin_edges = squeeze(pts_bin_edges(i, i_col, :));

                MI.bin_edges = MI.bin_edges(~isinf(MI.bin_edges) & ~isnan(MI.bin_edges));

                bin_edge_str = [num2str(min(MI.bin_edges),'%0.0f');...
                               repmat({''}, n_bins_used-1, 1);...
                               num2str(max(MI.bin_edges),'%0.0f')]';

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
                x_lbl = xlabel(cfg.nice_plt_lbls(i_col), 'Interpreter','none'); hold on
                set(gca, 'Fontsize', 8)

                 if i_row ~= i_col
                    x_lbl.Position(2) = -4.8;
                    set(gca, 'Fontsize', 8)

                 else
                    x_lbl.Position(2) = x_lbl.Position(2) *1.1;

                 end
            end
            if i_col == 1
                y_lbl = ylabel(cfg.nice_plt_lbls(i_row), 'Interpreter','none', 'Rotation', 0); hold on

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
 
    sgtitle([cfg.pt_lbls{i},' ',sg_tit_txt], 'FontSize', 12);

    exportgraphics(gcf, [save_dir,'intermediate_steps/stage0',cfg.pts{i}, '_CondProb_matrices.png'])
end


%% visualize pair-wise MI.bits as matrix
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

for i= 1 :length(cfg.pts)

    figure('Units', 'Inches', 'Position', [0, 0, 4, 3])

    MI.bits      = squeeze(pts_MI_in_per(i, :, :));
    perm_MI = squeeze(pts_MI_perm_p(i,:,:));

    % Bonferroni correction on the subject-level
%     corr_p = round(0.05/N_comparisons, 5);
%     MI.bits(perm_MI >= corr_p ) = nan;

    MI.bits(perm_MI > 0.05 ) = 0;

    MI.bits      = tril(MI.bits) + triu(nan(size(MI.bits)));

    plt_MI = nan(size(MI.bits,1), size(MI.bits,2));
    plt_MI(1:size(MI.bits,1)-1, 1:size(MI.bits,2)) = MI.bits(2:end,:);

    mat_pos = 0:length(cfg.tbl_pain_lbls)-1;

    s = pcolor(mat_pos, flip(mat_pos), plt_MI);

    s.EdgeColor = [1 1 1];      s.LineWidth = 1.5;
    
    clim(cb_lim);
     
    lbl_pos = (0:length(cfg.tbl_pain_lbls))+0.5;

    xticks(lbl_pos);    xticklabels(cfg.nice_plt_lbls);       
    yticks(lbl_pos);    yticklabels(flip(cfg.nice_plt_lbls));

    set(gca,'Fontsize',9,'TickLabelInterpreter','none','Color','w',...
        'TickLength', [0,0]);


    title(cfg.pt_lbls{i}, 'FontSize', 12);

    pts_plt_MI(i, :, :) = MI.bits;

    colormap([[.85 .85 .85];preset_colormap])
    set(gcf, 'InvertHardCopy', 'off')
    
    
    cb_hdl = colorbar;          cb_hdl.Limits = cb_lim; 
    
    cb_hdl.Label.String         = "Percent of max possible MI.bits";
    cb_hdl.Label.Rotation       = 270;
    cb_hdl.Label.Position(1)    = 3.3;
    cb_hdl.FontSize             = 10;
    cb_hdl.Label.FontSize       = 12;

    exportgraphics(gcf, fullfile(save_dir,['stage0',cfg.pts{i}, '_MI_pain_matrix.png']));

end

%% visualize the group-level MI.bits mean and prevelance of sig. subject-level comparisons
grp_fig =  figure('Units', 'Inches', 'Position', [0, 0, 4, 3]);


MI.bits   = squeeze(mean(pts_plt_MI,'omitnan'));
prop = flip(squeeze(sum(pts_plt_MI ~=0,1)));


plt_MI = nan(size(MI.bits,1), size(MI.bits,2));
plt_MI(1:size(MI.bits,1)-1, 1:size(MI.bits,2)) = MI.bits(2:end,:);

mat_pos = 0:length(cfg.tbl_pain_lbls)-1;

s = pcolor(mat_pos, flip(mat_pos), plt_MI);         clim(cb_lim);

s.EdgeColor = [1 1 1];          s.LineWidth = 1.5;
 

lbl_pos = (0:length(cfg.tbl_pain_lbls))+0.5;

xticks(lbl_pos);    xticklabels(cfg.nice_plt_lbls);       
yticks(lbl_pos);    yticklabels(flip(cfg.nice_plt_lbls));


prop_MI = nan(size(prop,1), size(prop,2));
prop_MI(1:size(MI.bits,1)-1, 2:size(MI.bits,2)) = prop(1:end-1, 2:end);

set(gcf, 'InvertHardCopy', 'off')
set(gca,'Fontsize', 9,'TickLabelInterpreter','none',...
    'Color','w', 'TickLength', [0,0]);

colormap(preset_colormap);
cb_hdl = colorbar;          cb_hdl.Limits = cb_lim; 

cb_hdl.Label.String         = "Percent of max possible MI.bits";
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

exportgraphics(gcf, [save_dir,'stage0GroupLevel_MI_pain_matrix.png']);

% %%
% metric_any_MI = nan(length(cfg.pts), length(cfg.tbl_pain_lbls));
% 
% high_entr_vars = {'painVAS', 'MPQ', 'PBD'};
% 
% i_PBD = contains(cfg.tbl_pain_lbls, 'PBD');
% 
% for i= 1 :length(cfg.pts)
% 
%     % make only lower triangle w/ p-values, but the rest w/ nans
%     perm_MI = squeeze(pts_MI_perm_p(i,:,:));
% 
%     for j = 1 : length(cfg.tbl_pain_lbls)
% 
%        i_high_entr_oi=  contains(high_entr_vars, cfg.tbl_pain_lbls(j));
% 
%         if any(i_high_entr_oi)
%            
%             i_comp = contains(cfg.tbl_pain_lbls, high_entr_vars(~i_high_entr_oi));
% 
%             metric_any_MI(i, j) = sum(perm_MI(i_comp, j) <= 0.05);
%             
%         end
%     end
% end
% 
% metric_any_MI(metric_any_MI>2) = 2;
% 

end