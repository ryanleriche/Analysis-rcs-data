function [pt_pain_space, vararout] ...
    ...
    = plot_pain_space(cfg, REDcap)
% pull relevant pts REDcap table
redcap  = REDcap.(cfg.pt_id);

filt_fig = figure('Units', 'Inches', 'Position', [0, 0, 15, 10]);

[redcap, date_range] = date_parser(cfg, redcap);

ds  =    datestr(date_range,'dd-mmm-yyyy');
title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);
hold on


redcap.MPQaff       = sum([redcap.MPQsickening, redcap.MPQfearful, redcap.MPQcruel,  redcap.MPQtiring],2,'omitnan');
redcap.MPQsom       = redcap.MPQtotal - redcap.MPQaff;


% only KEEP if all VAS reports do NOT equal 50

i_trim_VAS_any          = ~(redcap.unpleasantVAS == 50 | redcap.painVAS == 50 | redcap.worstVAS == 50);

RCSXX_trim_VAS_any      = redcap(i_trim_VAS_any,:);
RCSXX_rmv_VAS_any       = redcap(~i_trim_VAS_any,:);
                        


% only KEEP reports where 3 VAS reports do NOT equal 50, and in cases where
% at least one VAS report equals 50 keep it ONLY if NRS equals 5 (i.e., a true
% neutral report) AND none of the VAS reports are NaN

% only KEEP reports where 3 VAS reports do NOT equal 50, and in cases where
% at least one VAS report equals 50 keep it ONLY if NRS equals 5 (i.e., a true
% neutral report)

i_trim_VAS_all    = ~(redcap.unpleasantVAS == 50 & redcap.painVAS == 50 & redcap.worstVAS == 50)...
                    &...
                        ~(...
                        (redcap.unpleasantVAS == 50 | redcap.painVAS == 50 | redcap.worstVAS == 50)...
                        & redcap.mayoNRS ~= 5);


RCSXX_trim_VAS_all   = redcap(i_trim_VAS_all , :);

RCSXX_rmv_VAS_all    = redcap(~i_trim_VAS_all , :);


prop_any_VAS_remain  = height(RCSXX_trim_VAS_any) / height(redcap);

prop_all_VAS_remain  = height(RCSXX_trim_VAS_all) / height(redcap);


scatter3(redcap,'unpleasantVAS','painVAS','MPQtotal','filled', ...
    'ColorVariable', 'mayoNRS');


xlim([1 100]); ylim([1 100]); zlim([0 45]);

c = colorbar;
c.Label.String = 'mayoNRS';                 c.Limits = [1,10];

text(50, 15, 40, ...
    [...
        'Prop(remain w/o ANY VAS = 50): ', num2str(prop_any_VAS_remain),...
        newline,...
        'Prop(remain w/o ALL VAS = 50): ', num2str(prop_all_VAS_remain)...
    ],...
    'FontSize',14);


scatter3(RCSXX_rmv_VAS_any,'unpleasantVAS','painVAS','MPQtotal','filled','MarkerFaceColor', [.7 .7 .7]);

scatter3(RCSXX_rmv_VAS_all,'unpleasantVAS','painVAS','MPQtotal','filled','MarkerFaceColor','k');

format_plot()

  
%%
% ignore unanswered surveys
i_missing          = isnan(redcap.mayoNRS) | isnan(redcap.painVAS);
i_mpq_answered     = ~((redcap.painVAS >=40 | redcap.mayoNRS >=4) & redcap.MPQtotal ==0);


i_kept             = i_trim_VAS_any & ~i_missing & i_mpq_answered;


redcap_kept         = redcap(i_kept, :);
redcap_z            = redcap_kept; 
    


zscor_xnan   = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));




% z-score columns that contain pain metrics
pain_metrics = redcap_z.Properties.VariableNames(...
                contains(redcap_z.Properties.VariableNames, {'NRS', 'VAS', 'MPQ'}));

for i = 1 : length(pain_metrics)

   redcap_z.(pain_metrics{i}) = zscor_xnan(redcap_z.(pain_metrics{i}));

end

    
z_space_fig = figure('Units', 'Inches', 'Position', [0, 0, 15, 10]);

scatter3(redcap_z,'unpleasantVAS','painVAS','MPQtotal','filled', ...
    'ColorVariable', 'mayoNRS');

title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:),...
    newline, ' z-scored pain reports'], 'Fontsize',12);
hold on

c = colorbar;
c.Label.String = 'mayoNRS';                 c.Limits = [-3,3];
format_plot()

% considered use PCs, but clustering doesn't account for comp variance
if cfg.pca == true
    beh_pc = [];
    [beh_pc.coeff, beh_pc.score, beh_pc.latent, beh_pc.tsquared, beh_pc.explained, beh_pc.mu] ...
        = pca(redcap_z{:,2:7},'algorithm','als');
    
    % n_comps = find(cumsum(beh_pc.explained) >= 95, 1, 'first');
    
    vararout{1} = beh_pc;
    for i = 1 : 3
        % Plot the line, the original data, and their projection to the line.
        t = [min(beh_pc.score(:,i))-.2, max(beh_pc.score(:,i))+.2];
    
        endpts = [beh_pc.mu + t(1)*beh_pc.coeff(:,i)'; beh_pc.mu + t(2)*beh_pc.coeff(:,i)']*beh_pc.explained(i)/100;
        plot3(endpts(:,1),endpts(:,2),endpts(:,3),'LineWidth', 4, 'Color', 'k');
    
        hold on
    end
    
    
    legend({'',['PC 1', ' Variance: ', num2str(beh_pc.explained(1)), ' %'], ...
                ['PC 2', ' Variance: ', num2str(beh_pc.explained(2)), ' %'],...
                ['PC 3', ' Variance: ', num2str(beh_pc.explained(3)), ' %']});
    
    
    xlim([-4 4]); ylim([-4 4]); zlim([-4 4]);
    
    c = colorbar;
    c.Label.String = 'mayoNRS';                 c.Limits = [-4,4];
    
    format_plot();
end
%%
pain_metrics = redcap_z.Properties.VariableNames(...
                contains(redcap_z.Properties.VariableNames, {'NRS', 'VAS', 'MPQtotal'}));

switch cfg.pt_id
    
    case 'RCS06'
        [dec_fig, beh_cl.i_cl, beh_cl.i_halo]  = cluster_dp(redcap_z{:,[2:8, 27,28]});

    case 'RCS07'
        [dec_fig, beh_cl.i_cl, beh_cl.i_halo]  = cluster_dp(redcap_z{:,[2:8, 12]});

    case 'stage0RCS05'

        pain_metrics = redcap_z.Properties.VariableNames(...
                contains(redcap_z.Properties.VariableNames, {'painVAS','unpleasantVAS'}));

        [dec_fig, beh_cl.i_cl, beh_cl.i_halo]  = cluster_dp(redcap_z{:,pain_metrics});

    otherwise
        [dec_fig, beh_cl.i_cl, beh_cl.i_halo]  = cluster_dp(redcap_z{:,pain_metrics});

end

sgtitle(cfg.pt_id, 'Fontsize', 16)


redcap_z.i_clusters     = beh_cl.i_cl';
redcap_z.i_halo         = beh_cl.i_halo';

%%
cl_fig = figure('Units', 'Inches', 'Position', [0, 0, 15, 10]);


cmap = colormap(brewermap([],"Dark2"));
switch cfg.pt_id
    case 'stage0RCS05'
        for i = 1 : length(unique(beh_cl.i_cl))

            scatter(redcap_z(beh_cl.i_cl == i ,:),'painVAS','unpleasantVAS','filled',...
                 'MarkerFaceColor', cmap(i,:),'MarkerEdgeColor', cmap(i,:));
            hold on

        end
        

        text(1.5, 0, ['Apparent VAS', newline,...
            'affective-somatosensory dissociation'], ...
            'FontSize',14);

        format_plot();

    otherwise
        for i = 1 : length(unique(beh_cl.i_cl))

        scatter3(redcap_z(beh_cl.i_cl == i ,:),'painVAS','unpleasantVAS','MPQtotal','filled',...
             'MarkerFaceColor', cmap(i,:),'MarkerEdgeColor', cmap(i,:));
        hold on
        end

        zlim([-3 3]);
end





%    scatter3(redcap_z(beh_cl.i_halo == 0 ,:),'painVAS','unpleasantVAS','MPQtotal','filled',...
%      'MarkerFaceColor', 'k','MarkerEdgeColor', 'k');
%     

title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:),...
    newline, ' Pain State Clustering'], 'Fontsize',12);
hold on

% legend({'High Pain', 'Low Pain'});

xlim([-3 3]); ylim([-3 3]); 

format_plot();

redcap_kept.i_clusters     = beh_cl.i_cl';
redcap_kept.i_halo         = beh_cl.i_halo';




i_kept          = find(i_trim_VAS_any);
beh_cl.i_raw    = i_kept;
beh_pc.i_raw    = i_kept;


pt_pain_space.filt_fig      = filt_fig;
pt_pain_space.z_space_fig   = z_space_fig;
pt_pain_space.dec_fig       = dec_fig;
pt_pain_space.cl_fig        = cl_fig;
pt_pain_space.redcap_z      = redcap_z;
pt_pain_space.redcap_kept   = redcap_kept;


% save Stage 0 clusters (figures and REDcap) seperately from Stages 1, 2, and 3
if contains(cfg.pt_id, 'stage0')

    save_dir = [cd,'/plot_beh/figs/beh_only/', cfg.pt_id(end-4:end), '/stage_0/'];

    writetable(redcap_kept,...
       ['/Users/Leriche/Dropbox (UCSF Department of Neurological Surgery)/',...
       'SUBNETS Dropbox/Chronic Pain - Activa and Summit 2.0/DATA ANALYSIS/',...
       'Stage 0 ALL PATIENTS Redcap Records/arm1s_only_clustered/',...
        cfg.pt_id(end-4:end), '_redcap_records_kept'], 'FileType','spreadsheet');
else

    save_dir = [cd,'/plot_beh/figs/beh_only/', cfg.pt_id,'/cluster_anal/'];

    writetable(redcap_kept,...
       ['/Users/Leriche/Dropbox (UCSF Department of Neurological Surgery)',...
        '/UFlorida_UCSF_RCS_collab/Pain Reports/beh_clustered/',...
        cfg.pt_id, '_kept'], 'FileType','spreadsheet');

end

if ~exist(save_dir, 'dir')  
 
    mkdir(save_dir)
end

saveas(z_space_fig, [save_dir, 'z_space.png']);
saveas(filt_fig,    [save_dir, 'filt_VAS.png']);

saveas(dec_fig,     [save_dir, 'clus_decision.png']);
saveas(cl_fig,      [save_dir, 'clus_labels.png']);


function format_plot()  

    set(gca,'fontSize',14, 'TickLength', [0 0]); 
    grid on;    grid MINOR;      box off; 
    
end
end