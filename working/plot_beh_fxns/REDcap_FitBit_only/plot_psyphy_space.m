function plot_psyphy_space(pt_id, cfg, psy_phy_tbl)
    

% rename timestamp to 'time' 
vars         = psy_phy_tbl.Properties.VariableNames;
classes      = varfun(@class, psy_phy_tbl, 'OutputFormat','cell');

psy_phy_tbl = renamevars(psy_phy_tbl, ...
                    vars(strcmp(classes, 'datetime')), 'time');

psy_phy_tbl(all(isnan(psy_phy_tbl{:,cfg.var_oi}), 2), :) = [];


cmap = flip(brewermap([], "RdYlBu"));


%%% pull and parse
ds      = datestr([psy_phy_tbl.time(1), psy_phy_tbl.time(end)], 'dd-mmm-yyyy');

raw_fig = figure('Units', 'Inches', 'Position', [0, 0, 10, 7]);

title([pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16, 'Interpreter','none'); hold on

scatter3(psy_phy_tbl, cfg.view_scatter3{:},'filled', 'ColorVariable',  cfg.color_var);

colormap(gca, cmap);  c = colorbar; c.Label.String = cfg.color_var;

format_plot(cfg)
    
%% z-score metrics
z_psy_phy_tbl       = psy_phy_tbl;
zscor_xnan     = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

for i = 1 : length(cfg.var_oi)

   z_psy_phy_tbl.(cfg.var_oi{i}) = zscor_xnan(psy_phy_tbl.(cfg.var_oi{i}));

end

z_space_fig = figure('Units', 'Inches', 'Position', [0, 0, 10, 7]);


title([pt_id, ' (z-scored)', newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16, 'Interpreter','none'); hold on

scatter3(z_psy_phy_tbl, cfg.view_scatter3{:},'filled', 'ColorVariable',  cfg.color_var);


colormap(gca, cmap);  c = colorbar; c.Label.String = cfg.color_var;

clim([-3 3]);   xlim([-3 3]);   ylim([-3 3]);   ylim([-3 3]);

format_plot(cfg)

%% PCA for dimension reductionality
if cfg.pca == true
    
    psy_phy_oi  = z_psy_phy_tbl(:,cfg.var_oi);

    beh_pc = [];

    [beh_pc.coeff, beh_pc.score, beh_pc.latent,...
     beh_pc.tsquared, beh_pc.explained, beh_pc.mu] ...
     ...
        = pca(psy_phy_oi.Variables, 'algorithm','als');
    
    % organized PC comps neatly alongside raw pain metrics
    pc_lbls = cellfun(@(x) ['PC',x,'_score'],...
                replace(cellstr(num2str((1:length(beh_pc.explained))')),' ',''),...
                'UniformOutput', false) ;

    for i = 1: length(pc_lbls)
        psy_phy_oi.(pc_lbls{i}) = beh_pc.score(:, i);
    end


    % have lines of axes of above scatter3
    varnames = psy_phy_oi.Properties.VariableNames;
    i_var    = cellfun(@(x) find(strcmp(varnames, x)), cfg.view_scatter3);


    for i = 1 : 3
        % Plot the line, the original data, and their projection to the line.
        t = [min(beh_pc.score(:,i)), max(beh_pc.score(:,i))];
    
        endpts = [beh_pc.mu + t(1)*beh_pc.coeff(:,i)'; ...
                  beh_pc.mu + t(2)*beh_pc.coeff(:,i)'] *...
                  beh_pc.explained(i)/100;

        plot3(endpts(:,i_var(1)), endpts(:,i_var(2)), endpts(:,i_var(3)),'LineWidth', 4, 'Color', 'k');
    
        hold on
    end
    
    pc_var_txt = {['PC 1', ' Variance: ', num2str(beh_pc.explained(1)), ' %'], ...
               ['PC 2', ' Variance: ', num2str(beh_pc.explained(2)), ' %'],...
               ['PC 3', ' Variance: ', num2str(beh_pc.explained(3)), ' %']};
    
    legend([{''},pc_var_txt]); hold on

    cum_sum_vec = cumsum(beh_pc.explained);

    i_95_var  = find(cum_sum_vec >= 95);

    i_comp_oi = 1:i_95_var(1);

    clus_tbl = psy_phy_oi(:, compose('PC%g_score', i_comp_oi));
    
else
    % if not using PCA
    clus_tbl = z_psy_phy_tbl(:,cfg.var_oi);

end


    
%% cluster based on density peaks (Rodriguez and Laio, 2014 Science)
save_dir = fullfile(cfg.proc, pt_id, cfg.proc_subdir, '/');

if ~exist(save_dir, 'dir');      mkdir(save_dir);      end

saveas(raw_fig,     [save_dir, '1_raw_metrics.png']);
saveas(z_space_fig, [save_dir, '2_z_space.png']);


try
    [dec_fig, beh_cl.i_cl, beh_cl.i_halo]  = cluster_dp(cfg, clus_tbl.Variables);
    saveas(dec_fig,     [save_dir, '3_clus_decision.png']);
catch
    fprintf('%s | no clusters found\n', pt_id)
    return

end
%%
psy_phy_tbl.i_clusters     = beh_cl.i_cl';
psy_phy_tbl.i_halo         = beh_cl.i_halo';


cl_fig = figure('Units', 'Inches', 'Position', [0, 0, 8, 6]);


cmap_clus = colormap(brewermap([],"Dark2"));

for i = 1 : length(unique(beh_cl.i_cl))
    
    scatter3(psy_phy_tbl(beh_cl.i_cl == i ,:), ...
             cfg.view_scatter3{:},...
             'filled', 'MarkerFaceColor', cmap_clus(i,:),'MarkerEdgeColor', cmap_clus(i,:));
hold on
end

format_plot(cfg);

legend(compose('Cluster %g',unique(beh_cl.i_cl)), 'Location','best')

if cfg.pca

    txt_title = sprintf('%s\nclustered from first %g PCs (%.2f%% variance)',...
        pt_id, i_95_var(1),  cum_sum_vec(i_95_var(1)));

else
     txt_title = sprintf('%s\nclustered from z-scored metrics per se', pt_id);


end
title(txt_title,'Fontsize',16, 'Interpreter','none'); hold on
%% plot so-called psychophysio fingerprints (mean +- std of metrics w/n clusters

plt_def = {'TickLabelInterpreter', 'none', 'FontSize', 12, 'TickLength', [0,0]};

fp_tbl = table;

fp_fig = figure('Units', 'Inches', 'Position', [0, 0, 5, 7]);

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

%%
%%
% if .xlsx of clusterd data exists, overwrite it
source_xlsx = [save_dir, pt_id, '_clustered.xlsx'];

writetable(psy_phy_tbl, source_xlsx);

%%% now save figures as .pngs
saveas(cl_fig,      [save_dir, '4_clus_labels.png']);
saveas(fp_fig,      [save_dir, '5_clus_fingerprints.png']);



%%

function format_plot(cfg)  

    set(gca,'fontSize',16, 'TickLength', [0 0]); 
    grid on;    grid MINOR;      box off;

    xlabel(cfg.view_scatter3(1), 'Interpreter', 'none')
    ylabel(cfg.view_scatter3(2), 'Interpreter', 'none')
    zlabel(cfg.view_scatter3(3), 'Interpreter', 'none')
    
end
end