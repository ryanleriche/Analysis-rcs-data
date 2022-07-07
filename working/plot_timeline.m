function plot_timeline(RCSXX, epoch, varargin)

    if strcmp(epoch, 'AllTime') == 1

        figure('Units', 'Inches', 'Position', [0, 0, 15, 10])

        date_range           = [RCSXX.time(1), RCSXX.time(end)];

    elseif strcmp(epoch, 'PreviousDays') == 1

        n_days             = varargin{1};

        figure('Units', 'Inches', 'Position', [0, 0, 10, 7.5])
        
        [~, i_date] = min(abs(RCSXX.time - (datetime('now') - days(n_days))));

        date_range         = [RCSXX.time(i_date), RCSXX.time(end)];   


    elseif strcmp(epoch, 'DateRange') == 1

        



    end
    
    ds  =    datestr(date_range,'dd-mmm-yyyy');

    sgtitle([ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);
    


    t = RCSXX.time;

    subplot(311)

        if strcmp(epoch, 'AllTime') == 1
        
            plot(t, movmean([RCSXX.mayoNRS,RCSXX.worstNRS], 5),...
                'LineWidth', 2);
        
        elseif strcmp(epoch, 'PreviousDays') == 1
        
            plot(t, [RCSXX.mayoNRS,RCSXX.worstNRS], '.', 'MarkerSize',15);
        
        end
        
         ylabel('Numeric Rating Scale');     ylim([0,10])
             
         legend({'NRS: Intensity', 'NRS: Worst Intensity'}, 'Location','northeastoutside'); 
         format_plot();

    
    subplot(312)

        if strcmp(epoch, 'AllTime') == 1
        
            plot(t, movmean([RCSXX.painVAS,RCSXX.unpleasantVAS,RCSXX.worstVAS], 5),...
                'LineWidth', 2);
        
        elseif strcmp(epoch, 'PreviousDays') == 1
            
            plot(t, [RCSXX.painVAS,RCSXX.unpleasantVAS,RCSXX.worstVAS], '.', 'MarkerSize',15);
            
        end

         ylabel('Visual Analog Scale');         ylim([0,100]);
    
         legend({'VAS: Intensity',  'VAS: Worst Intensity', 'VAS: Unpleasantness'}, ...
             'Location','northeastoutside');    
         format_plot();
            

    subplot(313)
        MPQ_total     = RCSXX.MPQsum;
        MPQ_aff       = sum([RCSXX.MPQsickening, RCSXX.MPQfearful, RCSXX.MPQcruel],2,'omitnan');
        MPQ_som       = MPQ_total - MPQ_aff;

        if strcmp(epoch, 'AllTime') == 1
        
            plot(t, movmean([MPQ_som, MPQ_aff, MPQ_total, ], 5),...
                'LineWidth', 2);
            
        elseif strcmp(epoch, 'PreviousDays') == 1

            plot(t, [MPQ_som, MPQ_aff, MPQ_total, ], '.', 'MarkerSize',15);
            
        end
        
        ylabel('McGill Pain Questionaire');     ylim([0,45])
        
        legend({'MPQ: Somatic (0-33)', 'MPQ: Affective (0-12)', 'MPQ: Total (0-45)'}, ...
            'Location','northeastoutside');     
        format_plot();

    function format_plot()  

        set(gca,'fontSize',16, 'xlim', date_range , 'TickLength', [0 0]); 
        grid on;    grid MINOR;    legend boxoff;   box off
        
    end
end
