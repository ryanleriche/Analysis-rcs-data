function plot_hist(RCSXX, epoch, varargin)

    if strcmp(epoch, 'all-time') == 1

        figure('Units', 'Inches', 'Position', [0, 0, 15, 10])

        date_range           = [RCSXX.time(1), RCSXX.time(end)];

        

    elseif strcmp(epoch, 'PreviousDays') == 1

        n_days             = varargin{1};

        figure('Units', 'Inches', 'Position', [0, 0, 10, 7.5])
        
        [~, i_date] = min(abs(RCSXX.time - (datetime('now') - days(n_days))));

        date_range         = [RCSXX.time(i_date), RCSXX.time(end)];   
    end

    ds  =    datestr(date_range,'dd-mmm-yyyy');

    sgtitle([ds(1,:) ' to ' ds(2,:)]);
    
    subplot(3,2,[1,2])

        histogram(RCSXX.mayoNRS);    hold on;   histogram(RCSXX.worstNRS);

         ylabel('Numeric Rating Scale');
             
         legend({'NRS: Intensity', 'NRS: Worst Intensity'}, 'Location','northeastoutside'); 

         format_plot(); legend boxoff; 
            
    subplot(3,2,[3,4])

        histogram(RCSXX.painVAS,'BinWidth', 2);    hold on;   histogram(RCSXX.unpleasantVAS,'BinWidth', 2);

        ylabel('Visual Analog Scale');       
        legend({'VAS: Intensity', 'VAS: Unpleasantness'}, ...
             'Location','northeastoutside'); 

        format_plot(); legend boxoff; 

    subplot(325)


        MPQ_total     = RCSXX.MPQsum;
        MPQ_aff       = sum([RCSXX.MPQsickening, RCSXX.MPQfearful, RCSXX.MPQcruel],2,'omitnan');
        MPQ_som       = MPQ_total - MPQ_aff;


        histogram(MPQ_som); 

        ylabel('MPQ: Somatic (0-33)'); 

        format_plot;
        
     subplot(326)
        histogram(nonzeros(MPQ_aff));    xlim([0 45])
        
        ylabel('MPQ: Affective (0-12)','FontWeight', 'bold');

        format_plot();


     

    function format_plot()  

        set(gca,'defaultAxesFontSize',12, 'TickLength', [0 0]); 
        grid on;    grid MINOR;      box off
        
    end
end
