% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
github_dir      = '/Users/Leriche/Github/';


cd([github_dir, 'Analysis-rcs-data/working']);         

addpath(genpath([github_dir, 'Analysis-rcs-data/']));


% where DropBox desktop is saved locally
dropbox_dir     = '/Users/Leriche/Dropbox (UCSF Department of Neurological Surgery)/personal/collab/Presidio/behavioral_data/';

%%
csv_PR_tbl         = struct2table(dir([dropbox_dir,'*.csv']));

for i = 1: height(csv_PR_tbl)

    pt_id        = csv_PR_tbl{i, 'name'}(1:end-4);

    tmp_pr_beh   = readtable([ csv_PR_tbl{i, 'folder'},'/' csv_PR_tbl{i, 'name'}]);
    vars         = tmp_pr_beh.Properties.VariableNames;


    classes      = varfun(@class, tmp_pr_beh, 'OutputFormat','cell');
       
    i_double     = find(strcmp(classes, 'double'));
    i_rmv        = all(isnan(tmp_pr_beh{:,i_double}));

    all_nan_vars = vars(i_double(i_rmv));

    tmp_pr_beh   = removevars(tmp_pr_beh, all_nan_vars);


    classes      = varfun(@class, tmp_pr_beh, 'OutputFormat','cell');


    tmp_pr_beh{:,strcmp(classes, 'datetime')} =...
    tmp_pr_beh{:,strcmp(classes, 'datetime')} + calyears(2000);

    PR_beh.(pt_id) = tmp_pr_beh;

end
%%
cfg                = [];
cfg.var_oi         = {'vas_anxiety', 'vas_depression', ...
                                     'hamd2_score',...
                      'vas_pain_daily',...
                      ...
                      'vas_energy'};

cfg.view_scatter3  =  {'vas_anxiety','vas_energy', 'vas_pain_daily'};


%pain_arousal_lbs   = {'vas_energy', 'sss', 'vas_pain_daily'};

cfg.pca            = true;
cfg.CBDP_method    = 'top_two';


cfg.source_dir     = [dropbox_dir,'clustered/'];

pr_beh    = PR_beh.(pt_id); 


%%% pull and parse
% % only KEEP if all VAS reports do NOT equal 50 
% i_trim_VAS_any       = any(redcap_vas == 50, 2) & redcap.mayoNRS ~= 5;
% 
% redcap_trim_VAS_any  = redcap(~i_trim_VAS_any,:);
% redcap_rmv_VAS_any   = redcap(i_trim_VAS_any,:);
%                         
% 
% % only KEEP reports where 3 VAS reports do NOT equal 50
%i_trim_VAS_all        = all(redcap_vas ~= 50, 2) & ~i_trim_VAS_any;
% 
% redcap_trim_VAS_all   = redcap(i_trim_VAS_all, :);
% redcap_rmv_VAS_all    = redcap(~i_trim_VAS_all, :);
% 
% 
% per_any_VAS_remain  = 100 * (1 - (height(redcap_trim_VAS_any) / height(redcap)));
% per_all_VAS_remain  = 100 * (1 - (height(redcap_trim_VAS_all) / height(redcap)));


%ds               =    datestr(date_range,'dd-mmm-yyyy');

ds   = datestr([pr_beh.multi_daily_scales_timestamp(1),...
               pr_beh.multi_daily_scales_timestamp(end)],...
               'dd-mmm-yyyy');

filt_fig = figure('Units', 'Inches', 'Position', [0, 0, 10, 7]);

title([pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16); hold on


scatter3(pr_beh, cfg.view_scatter3{:},'filled', 'ColorVariable',  'hamd2_score');

% patch([50, 50, 50, 50], [0,0, 100,100],[0, 100, 100, 0],...
%     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
% 
% patch([0,0, 100,100], [50, 50, 50, 50], [0, 100, 100, 0],...
%     [.2, .2, .2], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 


% patch([0,0, 100,100], [0, 100, 100, 0], [50, 50, 50, 50],...
%     [.5, .5, .5], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 

plot3([0, 100], [50, 50], [50, 50], 'LineWidth', 2, 'Color', 'k');
plot3([50, 50], [0, 100], [50, 50], 'LineWidth', 2, 'Color', 'k');
plot3([50, 50], [50, 50], [0, 100], 'LineWidth', 2, 'Color', 'k');

xlabel(cfg.view_scatter3(1), 'Interpreter', 'none')
ylabel(cfg.view_scatter3(2), 'Interpreter', 'none')
zlabel(cfg.view_scatter3(3), 'Interpreter', 'none')
hold on

xlim([0 100]); ylim([0 100]); ylim([0 100]);


c = colorbar;
c.Label.String = 'HAM-D2';                 %c.Limits = [1,10];

format_plot()

%% PCA for dimension reductionality
%%% z-score neuropsychiatric metrics
z_pr_beh       = pr_beh;
zscor_xnan     = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

for i = 1 : length(cfg.var_oi)

   z_pr_beh.(cfg.var_oi{i}) = zscor_xnan(pr_beh.(cfg.var_oi{i}));

end

z_space_fig = figure('Units', 'Inches', 'Position', [0, 0, 15, 10]);


title([pt_id, newline, ds(1,:) ' to ' ds(2,:), 'z-scored'], 'Fontsize',16);


scatter3(z_pr_beh, cfg.view_scatter3{:},'filled', 'ColorVariable',  'hamd2_score');

xlabel(cfg.view_scatter3(1), 'Interpreter', 'none')
ylabel(cfg.view_scatter3(2), 'Interpreter', 'none')
zlabel(cfg.view_scatter3(3), 'Interpreter', 'none')
hold on



if cfg.pca == true

    tbl_oi =  z_pr_beh(:,cfg.var_oi); % <-- not redundant code

    beh_pc = [];

    [beh_pc.coeff, beh_pc.score, beh_pc.latent,...
     beh_pc.tsquared, beh_pc.explained, beh_pc.mu] ...
     ...
        = pca(tbl_oi.Variables, 'algorithm','als');
    
    vararout{1} = beh_pc;

    % organized PC comps neatly alongside raw pain metrics
    pc_lbls = cellfun(@(x) ['PC',x,'_score'],...
                cellstr(num2str((1:length(beh_pc.explained))')),...
                'UniformOutput', false) ;

    for i = 1: length(pc_lbls)
        z_pr_beh.(pc_lbls{i}) = beh_pc.score(:, i);
    end


    % have lines of axes of above scatter3
    varnames = tbl_oi.Properties.VariableNames;
    
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
    
    
    legend({'',['PC 1', ' Variance: ', num2str(beh_pc.explained(1)), ' %'], ...
               ['PC 2', ' Variance: ', num2str(beh_pc.explained(2)), ' %'],...
               ['PC 3', ' Variance: ', num2str(beh_pc.explained(3)), ' %']});
    
    
    xlim([-3 3]); ylim([-3 3]); zlim([-3 3]);
    
    c = colorbar;
    c.Label.String = 'HAM-D2';                 c.Limits = [-3 3];
    format_plot();
end



%% cluster based on density peaks (Rodriguez and Laio, 2014 Science)
%set(0,'DefaultFigureVisible','on')
[dec_fig, beh_cl.i_cl, beh_cl.i_halo]  = cluster_dp(cfg, pr_beh{:,cfg.var_oi});

%set(0,'DefaultFigureVisible','off')
%%

sgtitle(pt_id, 'Fontsize', 16, 'Interpreter', 'none')

pr_beh.i_clusters     = beh_cl.i_cl';
pr_beh.i_halo         = beh_cl.i_halo';


if length(unique(beh_cl.i_cl')) == 2

    if mean(pr_beh.vas_pain_daily(beh_cl.i_cl == 1), 'omitnan') > ...
       mean(pr_beh.vas_pain_daily(beh_cl.i_cl == 2), 'omitnan')

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

cl_fig = figure('Units', 'Inches', 'Position', [0, 0, 15, 10]);

cmap = colormap(brewermap([],"Dark2"));

for i = 1 : length(unique(beh_cl.i_cl))
    
    scatter3(pr_beh(beh_cl.i_cl == i ,:), ...
             cfg.view_scatter3{:},...
             'filled', 'MarkerFaceColor', cmap(i,:),'MarkerEdgeColor', cmap(i,:));
hold on
end


xlabel(cfg.view_scatter3(1), 'Interpreter', 'none')
ylabel(cfg.view_scatter3(2), 'Interpreter', 'none')
zlabel(cfg.view_scatter3(3), 'Interpreter', 'none')

hold on

xlim([-10 110]); ylim([-10 110]); ylim([-10 110]);

legend(leg_txt, 'Fontsize', 12);
format_plot();

title([pt_id, newline, ds(1,:) ' to ' ds(2,:),...
    newline, ' Pain State Clustering'], 'Fontsize',12,'Interpreter', 'none');



pr_beh.i_clusters     = beh_cl.i_cl';
pr_beh.i_halo         = beh_cl.i_halo';


% save as source data as .csv
if ~exist(cfg.source_dir, 'dir');        mkdir(cfg.source_dir);   end

% if csv exists, delete it first to prevent unexpected merging w/ previous
% version (RBL saw this 1/10/23)
source_csv = [cfg.source_dir, pt_id, '_clustered.csv'];

if exist(source_csv, 'file') == 2;        delete(source_csv);   end

writetable(pr_beh, source_csv);
% 
% cfg.fig_dir    = [cfg.fig_dir, '/', pt_id,'/'];
% % save figures as .pngs
% if ~exist(cfg.fig_dir, 'dir');      mkdir(cfg.fig_dir);      end
% 
% saveas(z_space_fig, [cfg.fig_dir, 'z_space.png']);
% saveas(filt_fig,    [cfg.fig_dir, 'filt_VAS.png']);
% 
% saveas(dec_fig,     [cfg.fig_dir, 'clus_decision.png']);
% saveas(cl_fig,      [cfg.fig_dir, 'clus_labels.png']);


function format_plot()  

    set(gca,'fontSize',16, 'TickLength', [0 0]); 
    grid on;    grid MINOR;      box off;

    
end