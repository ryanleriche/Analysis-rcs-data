function plot_timeline(cfg, redcap_RCSXX, db_RCSXXX)

    figure('Units', 'Inches', 'Position', [0, 0, 15, 10])

    [redcap_RCSXX, date_range] = date_parser(cfg, redcap_RCSXX);
    
    ds =        datestr(date_range,'dd-mmm-yyyy');

%% NRS
    if cfg.subplot == true
        subplot(311)
        sgtitle([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);

    else
        title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);
        hold on
    end


    if strcmp(cfg.dates, 'AllTime') == 1
    
        plot(redcap_RCSXX.time, movmean([redcap_RCSXX.mayoNRS,redcap_RCSXX.worstNRS], 5),...
            'LineWidth', 2);

        hold on; set(gca,'ColorOrderIndex',1);

        plot(redcap_RCSXX.time, [redcap_RCSXX.mayoNRS,redcap_RCSXX.worstNRS], '.', 'MarkerSize',5);
    

    else

   
        scatter(redcap_RCSXX.time, redcap_RCSXX.mayoNRS, 150, 'filled');   hold on;
        scatter(redcap_RCSXX.time, redcap_RCSXX.worstNRS, 75, 'filled');


    end
    
     ylabel('Numeric Rating Scale');     ylim([0,10]); yticks(1:2:10);
         
     legend({'NRS Intensity', 'NRS Worst Intensity'}, 'Location','northeastoutside'); 
     
     hold on
     overlay_stim(gca, cfg.stim_parameter, db_RCSXXX);
     format_plot();
     
%% VAS
    if cfg.subplot == true
        subplot(312)
    else
        figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
        title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16); 
        hold on

    end


    if strcmp(cfg.dates, 'AllTime') == 1
    
        plot(redcap_RCSXX.time, movmean([redcap_RCSXX.painVAS,redcap_RCSXX.unpleasantVAS,redcap_RCSXX.worstVAS], 5),...
            'LineWidth', 2);

             hold on; set(gca,'ColorOrderIndex',1);

        plot(redcap_RCSXX.time, [redcap_RCSXX.painVAS,redcap_RCSXX.unpleasantVAS,redcap_RCSXX.worstVAS], '.', 'MarkerSize',5);
    
    else
        
        scatter(redcap_RCSXX.time, redcap_RCSXX.painVAS, 150, 'filled');   
        hold on;
        scatter(redcap_RCSXX.time, redcap_RCSXX.unpleasantVAS, 100, 'filled');
        scatter(redcap_RCSXX.time, redcap_RCSXX.worstVAS, 75, 'filled');

        
    end

     ylabel('Visual Analog Scale');         ylim([0,100]);  yticks(0:20:100);

     legend({'VAS Intensity',  'VAS Worst Intensity', 'VAS Unpleasantness'}, ...
         'Location','northeastoutside'); 

     overlay_stim(gca, cfg.stim_parameter, db_RCSXXX);
     format_plot();
            
%% MPQ
    if cfg.subplot == true

        subplot(313)
    else
        figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
        title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16); 
        hold on

    end
    
    MPQ_total     = redcap_RCSXX.MPQsum;
    MPQ_aff       = sum([redcap_RCSXX.MPQsickening, redcap_RCSXX.MPQfearful, redcap_RCSXX.MPQcruel],2,'omitnan');
    MPQ_som       = MPQ_total - MPQ_aff;

    if strcmp(cfg.dates, 'AllTime') == 1
    
        plot(redcap_RCSXX.time, movmean(MPQ_total, 5),'LineWidth', 4.5);
        hold on
        plot(redcap_RCSXX.time, movmean([MPQ_som, MPQ_aff], 5),'LineWidth', 2);

         set(gca,'ColorOrderIndex',1);

        plot(redcap_RCSXX.time, MPQ_total, '.', 'MarkerSize',9);
        plot(redcap_RCSXX.time, [MPQ_som, MPQ_aff], '.', 'MarkerSize',5);

    else
 
        scatter(redcap_RCSXX.time, MPQ_total, 100, 'filled');   
        hold on;
        scatter(redcap_RCSXX.time, MPQ_som, 50, 'filled');
        scatter(redcap_RCSXX.time, MPQ_aff, 50, 'filled');

    end
    
    ylabel('McGill Pain Questionaire');     ylim([0,45]); yticks(0:15:45)
    
    legend({ 'MPQ Total (0-45)','MPQ Somatic (0-33)', 'MPQ Affective (0-12)'}, ...
        'Location','northeastoutside'); 

   % overlay_stim(cfg.stim_parameter)
    format_plot();
    overlay_stim(gca, cfg.stim_parameter, db_RCSXXX);

%% local functions

    function overlay_stim(ax, stim, RCSXXX_database)
        
        pompadour          = [104  2  63];
        deep_sea           = [0  70 60];
        deep_pink          = [192  11  111];
        persian_rose       = [239   0  150];
        persian_green      = [0  160  144];
        aqua_marine        = [0  220  181];
        carissma           = [255  149  186];
        dark_olive         = [61  60   4];
        aqua               = [95  255  222 ];
        royal_blue         = [0  60  134];
        indigo             = [89  10 135];
        vivid_purple       = [148   0 230];

        cb_colors   = [deep_pink; pompadour ; persian_rose ; persian_green;...
                      aqua_marine ; carissma ; deep_sea; dark_olive ;...
                      aqua ; royal_blue ; indigo ; vivid_purple ]./255;
        if strcmp(stim, 'contacts')
            
            cont_pairs = unique(RCSXXX_database.contacts);

            c = [colororder; cb_colors];

            for i = 2: length(cont_pairs) - 1            

                i_given_pair = strcmp(cont_pairs{i}, RCSXXX_database.contacts) &...
                                   RCSXXX_database.amp ~=0;

                starts = find(diff(i_given_pair)==1) + 1;
                ends   = find(diff(i_given_pair)==-1);

                % based of transitions plot every shading of a given
                % contact by itself

                if i_given_pair(end)

                    ends = [ends; length(i_given_pair)];

                end

               
                for j = 1:length(starts)

                    x = [RCSXXX_database.time(starts(j)), RCSXXX_database.time(ends(j)),...
                         RCSXXX_database.time(ends(j)), RCSXXX_database.time(starts(j))];
    
                    y = [0,0,ax.YLim(2), ax.YLim(2)];

                    if j == 1

                        patch(x, y, c(i,:),'FaceAlpha', .3);

                        ax.Legend.String{end} = cont_pairs{i};

                    else

                        patch(x, y, c(i,:),'FaceAlpha', .3,...
                            'HandleVisibility','off');
                    end

                    hold on
                end

                
            end
        end
    end

    function format_plot()  

        set(gca,'FontSize',16, 'xlim', date_range , 'TickLength', [0 0]); 
        grid on;    grid MINOR;   legend boxoff;    box off;

        t = datetime(cfg.stage_dates, 'TimeZone', '-07:00');

        if length(cfg.stage_dates) == 1
            xline(t, 'LineWidth', 2, 'Stage 1',...
                'HandleVisibility','off');

        elseif length(cfg.stage_dates) == 2
            xline(t, '-', {'Stage 1', 'Stage 2'},...
                'HandleVisibility','off');

        else
            xline(t, '-', ...
                {'Stage 1', 'Stage 2', 'Stage 3'},...
                'HandleVisibility','off');
        
        end
    end
end
