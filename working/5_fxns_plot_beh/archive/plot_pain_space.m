function [pt_pain_space, vararout] ...
    ...
    = plot_pain_space(cfg, REDcap)

% cfg             = [];
% cfg.dates       = 'AllTime';
% cfg.pca         = true;
% pts             = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};
% cfg.plt_VAS     = true;
% 
% i               = 1;
% 
% cfg.pt_id       = [ pts{i}];       
% redcap          = REDcap.(cfg.pt_id);

redcap                = REDcap.(cfg.pt_id);
%% pull and parse pt REDcap table

[redcap, date_range] = date_parser(cfg, redcap);

% ignore unanswered surveys
i_mpq_answered  = ~((redcap.painVAS >=30 | redcap.mayoNRS >=3 | ...
                     isnan(redcap.mayoNRS) | isnan(redcap.painVAS))...
                     & redcap.MPQtotal ==0);

redcap              = redcap(i_mpq_answered,:);


redcap.MPQaff       = sum([redcap.MPQsickening, redcap.MPQfearful, ...
                           redcap.MPQcruel,  redcap.MPQtiring], 2, 'omitnan');
redcap.MPQsom       = redcap.MPQtotal - redcap.MPQaff;


i_pain_vas           = contains(redcap.Properties.VariableNames, 'VAS') &...
                       ~contains(redcap.Properties.VariableNames, {'mood', 'relief','worst','best'});

redcap_vas           = redcap{:, i_pain_vas};


%% only KEEP if all VAS reports do NOT equal 50 (unless mayo NRS == 5)
i_trim_VAS_any       = any(redcap_vas == 50, 2) & redcap.mayoNRS ~= 5;

redcap_trim_VAS_any  = redcap(~i_trim_VAS_any,:);
redcap_rmv_VAS_any   = redcap(i_trim_VAS_any,:);
                        

% only KEEP reports where 3 VAS reports do NOT equal 50
i_trim_VAS_all        = all(redcap_vas ~= 50, 2) & ~i_trim_VAS_any;

redcap_trim_VAS_all   = redcap(i_trim_VAS_all, :);
redcap_rmv_VAS_all    = redcap(~i_trim_VAS_all, :);


per_any_VAS_remain  = 100 * (1 - (height(redcap_trim_VAS_any) / height(redcap)));
per_all_VAS_remain  = 100 * (1 - (height(redcap_trim_VAS_all) / height(redcap)));


pain_lbl         = {'unpleasantVAS','painVAS','MPQtotal'};
ds               =  datestr(date_range,'dd-mmm-yyyy');

if cfg.plt_VAS

    fig_raw = figure('Units', 'Inches', 'Position', [0, 0, 7, 5]);
    title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);
    hold on
    
    scatter3(redcap, pain_lbl{:},'filled', 'ColorVariable', 'mayoNRS');
    
    
    xlim([-10 110]); ylim([-10 110]); zlim([-2 45]);
    
    c = colorbar;
    c.Label.String = 'mayoNRS';                 c.Limits = [1,10];
    
    pain_lbl = {'unpleasantVAS','painVAS','MPQtotal'};
    
    scatter3(redcap_rmv_VAS_all, pain_lbl{:}, ...
        'filled', 'MarkerFaceColor', 'r', 'Marker', 'hexagram', 'SizeData', 125);
    
    scatter3(redcap_rmv_VAS_any, pain_lbl{:},...
        'filled','MarkerFaceColor', 'k', 'SizeData', 30);
     xticks(0:10:100); yticks(0:10:100);
    
    legend(...
            '',...
           ['Percent(ALL VASs = 50 or ANY VASs = 50, but NRS ~= 5): ', num2str(round(per_all_VAS_remain,1)),'%'],...
           ['Percent(ANY VASs = 50, but NRS ~= 5): ', num2str(round(per_any_VAS_remain,1)),'%'],...
           '','',...
           'FontSize',9, 'Location', 'northwest')
    
    format_plot()
end

%% considered plotting PBM hue along scatter3 of painVAS, MPQtotal, mayoNRS
%{
if cfg.plt_PBM

    pbm_fig = figure('Units', 'Inches', 'Position', [0, 0, 15, 10]);

    title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);
    hold on
    
    scatter3(redcap, 'painVAS','MPQtotal', 'mayoNRS', ...
        'filled', 'ColorVariable', 'PBMpixelvalue');
    
    
    xlim([-10 110]); ylim([-2 45]); zlim([-1 11]);
    
    c = colorbar;
    c.Label.String = '^{\Sigma Hue}/_{N pixels}';       
    
    pain_lbl = {'unpleasantVAS','painVAS','MPQtotal'};
    
    scatter3(redcap_rmv_VAS_all, pain_lbl{:}, ...
        'filled', 'MarkerFaceColor', 'r', 'Marker', 'hexagram', 'SizeData', 125);
    
    scatter3(redcap_rmv_VAS_any, pain_lbl{:},...
        'filled','MarkerFaceColor', 'k', 'SizeData', 30);
     xticks(0:10:100); yticks(0:10:100);
    
    legend(...
            '',...
           ['Percent(ALL VASs = 50 or ANY VASs = 50, but NRS ~= 5): ', num2str(round(per_all_VAS_remain,1)),'%'],...
           ['Percent(ANY VASs = 50, but NRS ~= 5): ', num2str(round(per_any_VAS_remain,1)),'%'],...
           '','',...
           'FontSize',13, 'Location', 'northwest')
    
    format_plot()
end
%}

%% z-score and visualize
if contains(cfg.pt_id, 'stage0')

    redcap_kept     = redcap;
else

    redcap_kept     = redcap(all(redcap_vas ~= 50, 2), :);
    
end
% z-score columns that contain pain metrics
pain_metrics    = redcap_kept.Properties.VariableNames(...
                    contains(redcap_kept.Properties.VariableNames, ...
                        {'NRS', 'VAS', 'MPQtotal'}) ...
                        ...
                   &  ~contains(redcap_kept.Properties.VariableNames, 'mood'));

i_all_nan      = all(isnan(redcap_kept{:, pain_metrics}),2);
redcap_kept    = redcap_kept(~i_all_nan, :);

redcap_z       = redcap_kept;

zscor_xnan     = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
for i = 1 : length(pain_metrics)

   redcap_z.(pain_metrics{i}) = zscor_xnan(redcap_kept.(pain_metrics{i}));

end

fig_z_space = figure('Units', 'Inches', 'Position', [0, 0, 7, 5]);

pain_lbl         = {'MPQtotal','painVAS'};

scatter(redcap_z, pain_lbl{:},'filled', 'ColorVariable', 'mayoNRS');

title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:),...
    newline, ' z-scored pain reports'], 'Fontsize',12);
hold on

c = colorbar;
c.Label.String = 'mayoNRS';    c.Limits = [-3,3];      format_plot();

legend off

%% PCA for dimension reductionality
if cfg.pca == true
    
    pain_oi  = redcap_z(:, pain_metrics);         % <-- not redundant code

    beh_pc = [];

    [beh_pc.coeff, beh_pc.score, beh_pc.latent,...
     beh_pc.tsquared, beh_pc.explained, beh_pc.mu] ...
     ...
        = pca(pain_oi.Variables, 'algorithm','als');
    
    vararout{1} = beh_pc;

    % organized PC comps neatly alongside raw pain metrics
    pc_lbls = cellfun(@(x) ['PC',x,'_score'],...
                cellstr(num2str((1:length(beh_pc.explained))')),...
                'UniformOutput', false) ;

    for i = 1: length(pc_lbls)
        redcap_kept.(pc_lbls{i}) = beh_pc.score(:, i);
    end


    % have lines of axes of above scatter3
    varnames = pain_oi.Properties.VariableNames;
    pain_lbl = {'MPQtotal','painVAS', 'mayoNRS'};
    i_var    = cellfun(@(x) find(strcmp(varnames, x)), pain_lbl);


    for i = 1 : 3
        % Plot the line, the original data, and their projection to the line.
        t = [min(beh_pc.score(:,i)), max(beh_pc.score(:,i))];
    
        endpts = [beh_pc.mu + t(1)*beh_pc.coeff(:,i)'; ...
                  beh_pc.mu + t(2)*beh_pc.coeff(:,i)'] *...
                  beh_pc.explained(i)/100;

        plot3(endpts(:,i_var(1)), endpts(:,i_var(2)), endpts(:,i_var(3)),'LineWidth', 4, 'Color', 'k');
    
        hold on
    end
    
    
    legend({'',['PC 1', ' Variance: ', num2str(beh_pc.explained(1)), ' %'], ...
               ['PC 2', ' Variance: ', num2str(beh_pc.explained(2)), ' %'],...
               ['PC 3', ' Variance: ', num2str(beh_pc.explained(3)), ' %']});
    
    
    xlim([-3 3]); ylim([-3 3]); zlim([-3 3]);
    
    c = colorbar;
    c.Label.String = 'mayoNRS';                 c.Limits = [-3,3];
    
    format_plot();
end

%% cluster based on density peaks (Rodriguez and Laio, 2014 Science)
%set(0,'DefaultFigureVisible','on')
switch cfg.pt_id
    
    case 'stage0RCS05'

        if cfg.VAS_only 

            cfg.pt_id = [cfg.pt_id, '_VAS_only'];

            pain_metrics = {'painVAS','unpleasantVAS'};
    
            [fig_dec, beh_cl.i_cl, beh_cl.i_halo]  = cluster_dp(redcap_z{:,pain_metrics});
        else
             [fig_dec, beh_cl.i_cl, beh_cl.i_halo]  = cluster_dp(cfg, redcap_z{:, pain_metrics});
        end

    otherwise
        try
            [fig_dec, beh_cl.i_cl, beh_cl.i_halo]  = cluster_dp(cfg, redcap_z{:, pain_metrics});
        catch

            fig_dec = gcf;

            beh_cl.i_cl     = ones(1, height(redcap_z));
            beh_cl.i_halo   = beh_cl.i_cl;


        end

end
%set(0,'DefaultFigureVisible','off')
%%

sgtitle(cfg.pt_id, 'Fontsize', 16, 'Interpreter', 'none')


redcap_kept.i_clusters     = beh_cl.i_cl';
redcap_kept.i_halo         = beh_cl.i_halo';

if length(unique(beh_cl.i_cl')) == 2

    if mean(redcap_kept.painVAS(beh_cl.i_cl == 1), 'omitnan') > ...
       mean(redcap_kept.painVAS(beh_cl.i_cl == 2), 'omitnan')

        i_2 = beh_cl.i_cl == 2;
        i_1 = beh_cl.i_cl == 1;

        beh_cl.i_cl(i_2) = 1;
        beh_cl.i_cl(i_1) = 2;
    end
    leg_txt = {'low pain cluster', 'high pain cluster'};

else
    leg_txt = cell(length(unique(beh_cl.i_cl')),1);

    for i = 1 : length(unique(beh_cl.i_cl'))

        leg_txt(i) = {['cluster ', num2str(i)]};

    end
end

fig_clus = figure('Units', 'Inches', 'Position', [0, 0, 7, 5]);


cmap = colormap(brewermap(5,'Dark2'));

switch cfg.pt_id
    case 'stage0RCS05_VAS_only'
        for i = 1 : length(unique(beh_cl.i_cl))

            scatter(redcap_kept(beh_cl.i_cl == i ,:),'painVAS','unpleasantVAS','filled',...
                 'MarkerFaceColor', cmap(i,:),'MarkerEdgeColor', cmap(i,:));
            hold on

        end
        

        text(1.5, 0, ['Apparent VAS', newline,...
            'affective-somatosensory dissociation'], ...
            'FontSize',14);

    case 'RCS06'
        for i = 1 : length(unique(beh_cl.i_cl))

        scatter3(redcap_kept(beh_cl.i_cl == i ,:),'nocVAS','MPQtotal', 'mayoNRS', 'filled',...
             'MarkerFaceColor', cmap(i,:),'MarkerEdgeColor', cmap(i,:));
        hold on
        end

           xlim([0 100]);       ylim([0 45]);        zlim([0 10]);
    otherwise
        for i = 1 : length(unique(beh_cl.i_cl))

        scatter3(redcap_kept(beh_cl.i_cl == i ,:),'painVAS','MPQtotal','mayoNRS','filled',...
             'MarkerFaceColor', cmap(i,:),'MarkerEdgeColor', cmap(i,:));
        hold on
        end
            xlim([0 100]);       ylim([0 45]);        zlim([0 10]);

       
end


legend(leg_txt, 'Fontsize', 12, 'Location','best');      legend boxoff
format_plot();

%    scatter3(redcap_z(beh_cl.i_halo == 0 ,:),'painVAS','unpleasantVAS','MPQtotal','filled',...
%      'MarkerFaceColor', 'k','MarkerEdgeColor', 'k');
%     

title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:),...
    newline, ' Pain State Clustering'], 'Fontsize',12,'Interpreter', 'none');
hold on


% legend({'High Pain', 'Low Pain'});

%xlim([-3 3]); ylim([-3 3]); 

redcap_kept.i_clusters     = beh_cl.i_cl';
redcap_kept.i_halo         = beh_cl.i_halo';


i_kept          = find(~i_trim_VAS_any);
beh_pc.i_raw    = i_kept;
redcap_kept.i_clusters     = beh_cl.i_cl';
redcap_kept.i_halo         = beh_cl.i_halo';


%%

if cfg.clus_fp_plots    % so-called "fingerprints" showing mean +- std of given plotting
    plt_def = {'TickLabelInterpreter', 'none', 'FontSize', 12, 'TickLength', [0,0]};

    fig_fp = figure('Units', 'Inches', 'Position', [0, 0, 5, 7]);

    tiledlayout(length(unique(beh_cl.i_cl)), 2, 'Padding','compact');

    for i = unique(beh_cl.i_cl)
    
        fp_tbl.mean{i} = mean(z_psy_phy_tbl{beh_cl.i_cl == i, cfg.var_oi}, 1, 'omitnan')';
    
        fp_tbl.std{i}  = std(z_psy_phy_tbl{beh_cl.i_cl == i, cfg.var_oi}, 1, 'omitnan')';
    
        plt_mean = [[fp_tbl.mean{i}, zeros(size(cfg.var_oi,2), 1)]; ...
                     zeros(1, size(cfg.var_oi,1) +1)];
    
    
        nexttile;       pcolor(plt_mean);
        
        title(compose('Cluster %g Mean', i));
        colormap(gca, cmap);            clim([-3, 3]); 
        tmp = colorbar;       tmp.Label.String = 'z-score';  tmp.Label.FontSize = 12;
        
        yticks((1:length(cfg.var_oi))+.5);
    
        yticklabels(cfg.var_oi); 
        
        set(gca, plt_def{:}); hold on;
    
        xticklabels([]);
    %%%
        plt_std = [[fp_tbl.std{i}, zeros(size(cfg.var_oi,2), 1)]; ...
                     zeros(1, size(cfg.var_oi,1) +1)];
    
        nexttile;       pcolor(plt_std); 
    
        title(compose('Cluster %g STD', i));
        colormap(gca, flipud(gray));               clim([0, 1.5]); 
        
        tmp = colorbar;        tmp.Label.String = 'z-score';    tmp.Label.FontSize = 12;
    
        yticklabels([]);    xticklabels([]);
    
        set(gca,  plt_def{:})
    end


    sgtitle(pt_id, 'Fontsize',14, 'Interpreter','none');
   
    fp_tbl.Properties.RowNames = compose('Cluster %g',unique(beh_cl.i_cl));

end

%%

pt_pain_space.fig_raw      = fig_raw;
pt_pain_space.fig_z_space   = fig_z_space;
pt_pain_space.fig_dec       = fig_dec;
pt_pain_space.fig_clus        = fig_clus;
pt_pain_space.redcap_z      = redcap_z;
pt_pain_space.redcap_kept   = redcap_kept;





if contains(cfg.pt_id, 'VAS_only')

    pt_id = cfg.pt_id(end-13:end);
    
else
    pt_id = cfg.pt_id(end-4:end);
end

if cfg.save_xlsx 
    % save as source data as .csv
    if ~exist(cfg.source_dir, 'dir');        mkdir(cfg.source_dir);   end
    
    % if csv exists, delete it first to prevent unexpected merging w/ previous
    % version (RBL saw this 1/10/23)
    source_xlsx = [cfg.source_dir, pt_id, '_kept.xlsx'];
    
    if exist(source_xlsx, 'file') == 2;        delete(source_xlsx);   end
    
    writetable(redcap_kept, source_xlsx);
end

    cfg.fig_dir    = [cfg.fig_dir, '/', pt_id,'/'];

% save figures as .pngs
if ~exist(cfg.fig_dir, 'dir');      mkdir(cfg.fig_dir);      end

saveas(fig_raw,    [cfg.fig_dir, '1_raw_metrics.png']);

saveas(fig_z_space, [cfg.fig_dir, '2_z_space.png']);


saveas(fig_dec,     [cfg.fig_dir, '3_clus_decision.png']);
saveas(fig_clus,      [cfg.fig_dir, '4_clus_labels.png']);

saveas(fig_fp,      [cfg.fig_dir, '5_clus_fingerprints.png']);



function format_plot()  

    set(gca,'fontSize',14, 'TickLength', [0 0]); 
    grid on;    box off;

    
end
end