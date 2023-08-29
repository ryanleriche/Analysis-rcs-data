function plot_timeline_s3(cfg, REDcap, fluct_sum_stats, varargin)
   
%%

redcap                  = REDcap.(cfg.pt_id);

% bring in pain fluctuation study as reference of baseline pain
if contains(cfg.pt_id, 'stage0')
    PFS_sum                 = fluct_sum_stats.(cfg.pt_id(7:end));
else
    PFS_sum                 = fluct_sum_stats.(cfg.pt_id);
end

figure('Units', 'Inches', 'Position', [0, 0, 13, 10])

[~, date_range] = date_parser(cfg, redcap);

ds =        datestr(date_range,'dd-mmm-yyyy');

redcap = redcap(isbetween(redcap.time, date_range(1), date_range(2)), :);


r_cap_start = dateshift(redcap.time(1), "start",'day');
r_cap_end   = dateshift(redcap.time(end), "end",'day');


r_cap_daily      = table;
r_cap_daily.time = (r_cap_start:duration('24:00:00'):r_cap_end)';


r_cap_vars = redcap.Properties.VariableNames;
r_cap_vars = r_cap_vars(~strcmp(r_cap_vars, 'time'));

r_cap_daily{:, r_cap_vars} = nan;

for i = 1 : height(r_cap_daily) -1

    i_r_cap = find(isbetween(redcap.time, r_cap_daily.time(i), r_cap_daily.time(i+1)));

    if isempty(i_r_cap)
        r_cap_daily{i, r_cap_vars} = nan;
    else
        r_cap_daily{i, r_cap_vars} = mean(redcap{i_r_cap, r_cap_vars}, 1, 'omitnan');
    end

end

cfg.ndays = height(r_cap_daily)-1;
%sum_stats =   calc_sum_stats(cfg, redcap);
c = brewermap(7 ,'Dark2');



%% allows exploratory analysis w/ consistent nice formatting
%% NRS
if cfg.subplot == true
    subplot(311);
    sgtitle([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);


elseif contains(cfg.pt_id, 'stage0')
    title(cfg.pt_lbl);
else
    title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);
end
hold on

switch cfg.dates
    case 'AllTime'

        plot(redcap.time, movmean([redcap.mayoNRS,redcap.worstNRS], 5),...
            'LineWidth', 2, 'HandleVisibility','off');


        plot(redcap.time, [redcap.mayoNRS,redcap.worstNRS], '.', 'MarkerSize',5);

    case {'DateRange', 'PreviousDays'}

        scatter(redcap.time, redcap.mayoNRS, 125, 'filled','MarkerFaceAlpha', 0.6,...
            'MarkerFaceColor', c(1,:),'MarkerFaceColor', c(1,:));    
        
    switch cfg.pt_id
        case 'RCS06'

            scatter(redcap.time, redcap.nocNRS, 75, 'filled','MarkerFaceAlpha', 0.6); 
            scatter(redcap.time, redcap.npNRS, 50, 'filled', 'MarkerFaceAlpha', 0.6);

       
        case 'RCS07'

            scatter(redcap.time, redcap.unpNRS, 50, 'filled','MarkerFaceAlpha', 0.6); 
            scatter(redcap.time, redcap.leftarmNRS, 65, 'filled', 'MarkerFaceAlpha', 0.6);

            scatter(redcap.time, redcap.leftlegNRS, 50 , 'filled', 'MarkerFaceAlpha', 0.6);  
            scatter(redcap.time, redcap.leftfaceNRS, 35, 'filled', 'MarkerFaceAlpha', 0.6);


        case {'RCS02', 'RCS04', 'RCS05'}
            scatter(redcap.time, redcap.worstNRS, 75, 'filled', 'MarkerFaceAlpha', 0.6, ...
                'MarkerFaceColor', c(2,:),'MarkerFaceColor', c(2,:));  
    end
end

ylabel('Numeric Rating Scale');     ylim([0,10]); yticks(1:2:10);

% for reference plot mean±std of pain fluctuation study (prior to
% temporary trial period)
tmp_fluct  = PFS_sum(:, 'mayoNRS');

yline(mean(r_cap_daily.mayoNRS, 'omitnan'), '-', 'Color', c(1,:), 'LineWidth', 1.5)

fluct_t    = sort(repmat([date_range(1), date_range(2)],1,2));
fluct_plt  = [tmp_fluct{'mean',1}-tmp_fluct{'std',1},    tmp_fluct{'mean',1}+tmp_fluct{'std',1},...
                tmp_fluct{'mean',1}+tmp_fluct{'std',1},    tmp_fluct{'mean',1}-tmp_fluct{'std',1}];

patch(fluct_t, fluct_plt,'r','FaceAlpha',0.1,'EdgeColor', 'none');
yline(tmp_fluct{'mean',1}, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off')
 
yline(tmp_fluct{'half_improve',1}, '-', 'LineWidth', 1.5)
yline(tmp_fluct{'third_improve',1}, '--', 'LineWidth', 1.5)


switch cfg.pt_id
     case 'RCS06'

        legend({...
                'NRS Intensity','NRS Nociceptive','NRS Neuropathic',...
                '50%⭣in Pre-trial NRS intensity'...
                }, ...
                'Location','northoutside', 'NumColumns',3);

     case 'RCS07'

        legend({...
                'NRS Intensity','NRS Unpleasantness',...
                'NRS Left Arm', 'NRS Left Leg', 'NRS LeftFace',...
                '50%⭣in Pre-trial NRS intensity'...
                },...
                'Location','northoutside', 'NumColumns',3);
     
     case {'RCS02', 'RCS04', 'RCS05'}

        pre_mean_txt =  sprintf(...
                                 'Pre-trial NRS Int. mean ± std (%.2f ± %.2f)',...
                                 tmp_fluct{'mean',1}, tmp_fluct{'std',1});

        last_N_days_txt = sprintf('Mean NRS Int. last %g days (%.2f)', cfg.ndays, mean(r_cap_daily.mayoNRS, 'omitnan'));

        legend({...
                'NRS Int.' , 'NRS Worst Int.', last_N_days_txt,...
                pre_mean_txt,...
                sprintf('Pre-trial NRS Int. 50%%⭣ (%.2f)', tmp_fluct{'half_improve',1}),...
                sprintf('Pre-trial NRS Int. 30%%⭣ (%.2f)', tmp_fluct{'third_improve',1}),...
                },...
                'Location','northoutside', 'NumColumns',3);
 end
 
format_plot();

%% VAS
if cfg.subplot == true
    subplot(312);
else
    figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
    title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16); 
    hold on

end

% for reference plot mean±std of pain fluctuation study (prior to
% temporary trial period)

tmp_fluct  = PFS_sum(:, 'painVAS');

hold on

switch cfg.dates
    case 'AllTime'

        plot(redcap.time, movmean([redcap.painVAS,redcap.unpleasantVAS,redcap.worstVAS], 5),...
            'LineWidth', 2, 'HandleVisibility','off');

        plot(redcap.time, [redcap.painVAS,redcap.unpleasantVAS,redcap.worstVAS], '.', 'MarkerSize',5);
    
    case {'DateRange', 'PreviousDays'}

        scatter(redcap.time, redcap.painVAS, 100, 'filled', 'MarkerFaceAlpha', 0.6,...
            'MarkerFaceColor', c(1,:),'MarkerFaceColor', c(1,:));    
        

    switch cfg.pt_id
        case 'RCS06'

            scatter(redcap.time, redcap.nocVAS, 75, 'filled', 'MarkerFaceAlpha', 0.6);
            scatter(redcap.time, redcap.npVAS, 50, 'filled', 'MarkerFaceAlpha', 0.6);

        case {'RCS02', 'RCS04', 'RCS05'}

            scatter(redcap.time, redcap.worstVAS, 50, 'filled','MarkerFaceAlpha', 0.6,...
                'MarkerFaceColor', c(2,:),'MarkerFaceColor', c(2,:));    
            
            
            scatter(redcap.time, redcap.unpleasantVAS, 75, 'filled','MarkerFaceAlpha', 0.6,...
                'MarkerFaceColor', c(3,:),'MarkerFaceColor', c(3,:));   

        case  'RCS07'

            scatter(redcap.time, redcap.unpleasantVAS, 75, 'filled','MarkerFaceAlpha', 0.6);
            

            plot(redcap.time, movmean(redcap.moodVAS, [n_back 0], 'omitnan'),...
            'LineWidth', 2,'HandleVisibility','off'); hold on;

            scatter(redcap.time, redcap.moodVAS, 75, 'filled','MarkerFaceAlpha', 0.6);

    end
end

ylabel('Visual Analog Scale');         ylim([0,100]);  yticks(0:20:100);

fluct_t    = sort(repmat([date_range(1), date_range(2)],1,2));
fluct_plt  = [tmp_fluct{'mean',1}-tmp_fluct{'std',1},    tmp_fluct{'mean',1}+tmp_fluct{'std',1},...
            tmp_fluct{'mean',1}+tmp_fluct{'std',1},    tmp_fluct{'mean',1}-tmp_fluct{'std',1}];

yline(mean(r_cap_daily.painVAS, 'omitnan'), '-', 'Color', c(1,:), 'LineWidth', 1.5)

patch(fluct_t, fluct_plt,'r','FaceAlpha',0.1,'EdgeColor', 'none'); 
yline(tmp_fluct{'mean',1}, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off')

yline(tmp_fluct{'half_improve',1}, '-', 'LineWidth', 1.5)
yline(tmp_fluct{'third_improve',1}, '--', 'LineWidth', 1.5)


switch cfg.pt_id
     case 'RCS06'

        legend({...
                 'VAS Intensity', 'VAS Nociceptive','VAS Neuropathic',...
                 '50%⭣in Pre-trial VAS intensity'...
                },...
                'Location','northoutside', 'NumColumns',3);

     case 'RCS07'

        legend({...
                'VAS Intensity', 'VAS Unpleasantness','VAS Mood',...
                '50%⭣in Pre-trial VAS intensity'...
                },...
                'Location','northoutside', 'NumColumns',3);
     
     case {'RCS02', 'RCS04', 'RCS05'}

        pre_mean_txt =  sprintf(...
                         'Pre-trial VAS Int. mean ± std (%.2f ± %.2f)',...
                         tmp_fluct{'mean',1}, tmp_fluct{'std',1});

        last_N_days_txt = sprintf('Mean VAS Int. last %g days (%.2f)', cfg.ndays, mean(r_cap_daily.painVAS, 'omitnan'));

        legend({...
                'VAS Int.', 'VAS Unpleasantness','VAS Worst', last_N_days_txt,...
                pre_mean_txt,...
                sprintf('Pre-trial VAS Int. 50%%⭣ (%.2f)', tmp_fluct{'half_improve',1}),...
                sprintf('Pre-trial VAS Int. 30%%⭣ (%.2f)', tmp_fluct{'third_improve',1}),...
                },...
                'Location','northoutside', 'NumColumns',3);
end

format_plot();

        
%% MPQ

% for reference plot mean±std of pain fluctuation study (prior to
% temporary trial period)

if cfg.subplot == true; subplot(313);
else
    figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
    title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16); 
    hold on
end

MPQaff       = sum([redcap.MPQsickening, redcap.MPQfearful, redcap.MPQcruel, redcap.MPQtiring],2,'omitnan');
MPQsom       = redcap.MPQtotal - MPQaff;
hold on
switch cfg.dates
    case 'AllTime'

        plot(redcap.time, movmean(redcap.MPQtotal , 5),'LineWidth', 4.5,'HandleVisibility','off');

        plot(redcap.time, movmean([MPQsom, MPQaff], 5),'LineWidth', 2,'HandleVisibility','off');

        plot(redcap.time, redcap.MPQtotal , '.', 'MarkerSize',9);
        plot(redcap.time, [MPQsom, MPQaff], '.', 'MarkerSize',5);

    case {'DateRange', 'PreviousDays'}
 
        scatter(redcap.time, redcap.MPQtotal , 75, 'filled','MarkerFaceAlpha', 0.6,...
                'MarkerFaceColor', c(1,:),'MarkerFaceColor', c(1,:));   

        scatter(redcap.time, MPQsom, 50, 'filled', 'MarkerFaceAlpha', 0.6,...
                'MarkerFaceColor', c(2,:),'MarkerFaceColor', c(2,:));   

        scatter(redcap.time, MPQaff, 50, 'filled', 'MarkerFaceAlpha', 0.6,...
                'MarkerFaceColor', c(3,:),'MarkerFaceColor', c(3,:));   

end

ylabel('McGill Pain Questionaire');     ylim([0,45]);    yticks(0:9:45)

yline(mean(r_cap_daily.MPQtotal, 'omitnan'), '-', 'Color', c(1,:), 'LineWidth', 1.5)

% return pain flucutation study as reference range
tmp_fluct  = PFS_sum(:, 'MPQtotal');

fluct_t    = sort(repmat([date_range(1), date_range(2)],1,2));
fluct_plt  = [tmp_fluct{'mean',1}-tmp_fluct{'std',1},    tmp_fluct{'mean',1}+tmp_fluct{'std',1},...
            tmp_fluct{'mean',1}+tmp_fluct{'std',1},    tmp_fluct{'mean',1}-tmp_fluct{'std',1}];

patch(fluct_t, fluct_plt,'r','FaceAlpha',0.1,'EdgeColor', 'none'); hold on;
yline(tmp_fluct{'mean',1}, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off')


yline(tmp_fluct{'half_improve',1}, '-', 'LineWidth', 1.5)
yline(tmp_fluct{'third_improve',1}, '--', 'LineWidth', 1.5)

pre_mean_txt =  sprintf(...
             'Pre-trial MPQ Total Int. mean ± std (%.2f ± %.2f)',...
             tmp_fluct{'mean',1}, tmp_fluct{'std',1});

last_N_days_txt = sprintf('Mean MPQ Total last %g days (%.2f)', cfg.ndays, mean(r_cap_daily.MPQtotal, 'omitnan'));


legend({...
        'MPQ Total',  'MPQ Somatic','MPQ Affective',...
        last_N_days_txt, pre_mean_txt,...
        sprintf('50%%⭣in Pre-trial MPQ Total (%.2f)', tmp_fluct{'half_improve',1}),...
        sprintf('30%%⭣in Pre-trial MPQ Total (%.2f)', tmp_fluct{'third_improve',1}),...
        },...
        'Location','northoutside', 'NumColumns',3);

format_plot();

%% local functions
function format_plot()  
    legend boxoff;    box off;
    set(gca,'FontSize',13, 'xlim', date_range , 'TickLength', [0,0]);
end
end
