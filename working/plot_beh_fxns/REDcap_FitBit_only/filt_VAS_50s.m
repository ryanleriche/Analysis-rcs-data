function [r_cap] ...
    = filt_VAS_50s(pt_id, cfg, redcap) 


%% ignore unanswered surveys
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
%% plot criteria of VAS 50s
pain_lbl         = {'unpleasantVAS','painVAS','MPQtotal'};

ds      = datestr([redcap.time(1), redcap.time(end)], 'dd-mmm-yyyy');

fig_raw = figure('Units', 'Inches', 'Position', [0, 0, 7, 5]);
title([pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);
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
%%
r_cap.tbl      = redcap_trim_VAS_any;

vars           = r_cap.tbl.Properties.VariableNames;
i_vars         = contains(vars,  {'NRS', 'VAS', 'MPQtotal'}) & ~contains(vars, {'mood'});
r_cap.vars_oi  = vars(i_vars);


%% remove variables w/ only NaNs
i_all_nan = find(all(isnan(r_cap.tbl{:, r_cap.vars_oi}), 1));

r_cap.vars_oi(i_all_nan) = [];

%%
save_dir   = fullfile(cfg.proc, pt_id);
% save figures as .pngs
if ~isfolder(save_dir);     mkdir(save_dir);      end

saveas(fig_raw, fullfile(save_dir, '0_filter_VAS.png'));


%%
function format_plot()  

    set(gca,'fontSize',14, 'TickLength', [0 0]); 
    grid on;    box off;
  
end
end