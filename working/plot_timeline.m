function plot_timeline(cfg, RCSXX)

    figure('Units', 'Inches', 'Position', [0, 0, 15, 10])

    [RCSXX, date_range] = date_parser(cfg, RCSXX);
    
    ds =        datestr(date_range,'dd-mmm-yyyy');
    sgtitle([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);

    subplot(311)

        if strcmp(cfg.dates, 'AllTime') == 1
        
            plot(RCSXX.time, movmean([RCSXX.mayoNRS,RCSXX.worstNRS], 5),...
                'LineWidth', 2);

            hold on; set(gca,'ColorOrderIndex',1);

            plot(RCSXX.time, [RCSXX.mayoNRS,RCSXX.worstNRS], '.', 'MarkerSize',5);
        

        else

       
            scatter(RCSXX.time, RCSXX.mayoNRS, 150, 'filled');   hold on;
            scatter(RCSXX.time, RCSXX.worstNRS, 75, 'filled');

    
        end
        
         ylabel('Numeric Rating Scale');     ylim([0,10]); yticks(1:2:10);
             
         legend({'NRS Intensity', 'NRS Worst Intensity'}, 'Location','northeastoutside'); 
         format_plot();

    
    subplot(312)
    

        if strcmp(cfg.dates, 'AllTime') == 1
        
            plot(RCSXX.time, movmean([RCSXX.painVAS,RCSXX.unpleasantVAS,RCSXX.worstVAS], 5),...
                'LineWidth', 2);

                 hold on; set(gca,'ColorOrderIndex',1);

            plot(RCSXX.time, [RCSXX.painVAS,RCSXX.unpleasantVAS,RCSXX.worstVAS], '.', 'MarkerSize',5);
        
        else
            
            scatter(RCSXX.time, RCSXX.painVAS, 150, 'filled');   
            hold on;
            scatter(RCSXX.time, RCSXX.unpleasantVAS, 100, 'filled');
            scatter(RCSXX.time, RCSXX.worstVAS, 75, 'filled');

            
        end

         ylabel('Visual Analog Scale');         ylim([0,100]);  yticks(0:20:100);
    
         legend({'VAS Intensity',  'VAS Worst Intensity', 'VAS Unpleasantness'}, ...
             'Location','northeastoutside');    
         format_plot();
            

    subplot(313)
        MPQ_total     = RCSXX.MPQsum;
        MPQ_aff       = sum([RCSXX.MPQsickening, RCSXX.MPQfearful, RCSXX.MPQcruel],2,'omitnan');
        MPQ_som       = MPQ_total - MPQ_aff;

        if strcmp(cfg.dates, 'AllTime') == 1
        
            plot(RCSXX.time, movmean(MPQ_total, 5),'LineWidth', 4.5);
            hold on
            plot(RCSXX.time, movmean([MPQ_som, MPQ_aff], 5),'LineWidth', 2);

             set(gca,'ColorOrderIndex',1);

            plot(RCSXX.time, MPQ_total, '.', 'MarkerSize',9);
            plot(RCSXX.time, [MPQ_som, MPQ_aff], '.', 'MarkerSize',5);

        else
     
            scatter(RCSXX.time, MPQ_total, 100, 'filled');   
            hold on;
            scatter(RCSXX.time, MPQ_som, 50, 'filled');
            scatter(RCSXX.time, MPQ_aff, 50, 'filled');

        end
        
        ylabel('McGill Pain Questionaire');     ylim([0,45]); yticks(0:15:45)
        
        legend({ 'MPQ Total (0-45)','MPQ Somatic (0-33)', 'MPQ Affective (0-12)'}, ...
            'Location','northeastoutside');     
        format_plot();

    function format_plot()  

        set(gca,'FontSize',16, 'xlim', date_range , 'TickLength', [0 0]); 
        grid on;    grid MINOR;   legend boxoff;    box off;

        if length(cfg.stage_dates) == 1
            xline(datetime(cfg.stage_dates), 'LineWidth', 2, 'Stage 1',...
                'HandleVisibility','off');

        elseif length(cfg.stage_dates) == 2
            xline(datetime(cfg.stage_dates), '-', {'Stage 1', 'Stage 2'},...
                'HandleVisibility','off');

        else
            xline(datetime(cfg.stage_dates), '-', ...
                {'Stage 1', 'Stage 2', 'Stage 3'},...
                'HandleVisibility','off');
        
        end
    end
end
