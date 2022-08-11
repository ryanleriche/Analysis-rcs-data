function [filt_fig, z_space_fig, dec_fig, cl_fig,...
    z_RCSXX, vararout] = plot_versus(cfg, RCSXX)

filt_fig = figure('Units', 'Inches', 'Position', [0, 0, 15, 10]);

[~, RCSXX, date_range] = date_parser(cfg, RCSXX, []);

ds  =    datestr(date_range,'dd-mmm-yyyy');
title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);
hold on


RCSXX.MPQaff       = sum([RCSXX.MPQsickening, RCSXX.MPQfearful, RCSXX.MPQcruel],2,'omitnan');
RCSXX.MPQsom       = RCSXX.MPQtotal - RCSXX.MPQaff;


% only KEEP if all VAS reports do NOT equal 50

i_trim_VAS_any           = ~(RCSXX.unpleasantVAS == 50 | RCSXX.painVAS == 50 | RCSXX.worstVAS == 50);

RCSXX_trim_VAS_any      = RCSXX(i_trim_VAS_any,:);
RCSXX_rmv_VAS_any       = RCSXX(~i_trim_VAS_any,:);
                        


% only KEEP reports where 3 VAS reports do NOT equal 50, and in cases where
% at least one VAS report equals 50 keep it ONLY if NRS equals 5 (i.e., a true
% neutral report) AND none of the VAS reports are NaN

i_trim_VAS_all   = ~(RCSXX.unpleasantVAS == 50 & RCSXX.painVAS == 50 & RCSXX.worstVAS == 50)...
                    &...
                        ~(...
                        (RCSXX.unpleasantVAS == 50 | RCSXX.painVAS == 50 | RCSXX.worstVAS == 50)...
                        & RCSXX.mayoNRS ~= 5 ...
                        ...
                        & sum(isnan([RCSXX.unpleasantVAS, RCSXX.painVAS == 50, RCSXX.worstVAS]),2) == 0);


RCSXX_trim_VAS_all   = RCSXX(i_trim_VAS_all , :);

RCSXX_rmv_VAS_all    = RCSXX(~i_trim_VAS_all , :);


prop_any_VAS_remain  = height(RCSXX_trim_VAS_any) / height(RCSXX);

prop_all_VAS_remain  = height(RCSXX_trim_VAS_all) / height(RCSXX);


scatter3(RCSXX,'unpleasantVAS','painVAS','MPQtotal','filled', ...
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

pain_metrics = RCSXX.Properties.VariableNames;

% if unpleasantVAS, painVAS, and worstVAS are ALL 50s, ignore those reports

z_RCSXX      = RCSXX_trim_VAS_any;

zscor_xnan   = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));


for i = 2 : length(pain_metrics)

   z_RCSXX.(pain_metrics{i}) = zscor_xnan(z_RCSXX.(pain_metrics{i}));

end

    
z_space_fig = figure('Units', 'Inches', 'Position', [0, 0, 15, 10]);

scatter3(z_RCSXX,'unpleasantVAS','painVAS','MPQtotal','filled', ...
    'ColorVariable', 'mayoNRS');

title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:),...
    newline, ' z-scored pain reports'], 'Fontsize',12);
hold on

% considered use PCs, but clustering doesn't account for comp variance
if cfg.pca == true
    beh_pc = [];
    [beh_pc.coeff, beh_pc.score, beh_pc.latent, beh_pc.tsquared, beh_pc.explained, beh_pc.mu] ...
        = pca(z_RCSXX{:,2:7},'algorithm','als');
    
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

if strcmp(cfg.pt_id, 'RCS06')
    [dec_fig, beh_cl.i_cl, beh_cl.i_halo]       = cluster_dp(RCSXX_trim_VAS_any{:,[2:9, 11]});
else
    [dec_fig, beh_cl.i_cl, beh_cl.i_halo]       = cluster_dp(z_RCSXX{:,2:7});
end
sgtitle(cfg.pt_id, 'Fontsize', 16)


z_RCSXX.i_clusters     = beh_cl.i_cl';
z_RCSXX.i_halo         = beh_cl.i_halo';

cl_fig = figure('Units', 'Inches', 'Position', [0, 0, 15, 10]);


cmap = colormap(brewermap([],"Dark2"));

for i = 1 : length(unique(beh_cl.i_cl))

    scatter3(z_RCSXX(beh_cl.i_cl == i ,:),'unpleasantVAS','painVAS','MPQtotal','filled',...
         'MarkerFaceColor', cmap(i,:),'MarkerEdgeColor', cmap(i,:));
    hold on

 

end

   scatter3(z_RCSXX(beh_cl.i_halo == 0 ,:),'unpleasantVAS','painVAS','MPQtotal','filled',...
     'MarkerFaceColor', 'k','MarkerEdgeColor', 'k');
    

title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:),...
    newline, ' Pain State Clustering'], 'Fontsize',12);
hold on

% legend({'High Pain', 'Low Pain'});

xlim([-4 4]); ylim([-4 4]); zlim([-4 4]);

format_plot();

i_kept          = find(i_trim_VAS_any);
beh_cl.i_raw    = i_kept;
beh_pc.i_raw    = i_kept;

function format_plot()  

    set(gca,'fontSize',14, 'TickLength', [0 0]); 
    grid on;    grid MINOR;      box off; 
    
end
end