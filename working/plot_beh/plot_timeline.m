function plot_timeline(cfg, redcap_RCSXX, db_beh_RCSXX)

    figure('Units', 'Inches', 'Position', [0, 0, 15, 10])

    [db_beh_RCSXX, redcap_RCSXX, date_range] = date_parser(cfg, redcap_RCSXX, db_beh_RCSXX);
    
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

     if ~strcmp(cfg.stim_parameter, '')
         overlay_stim(gca, cfg.stim_parameter, db_beh_RCSXX);   
     end

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

     if ~strcmp(cfg.stim_parameter, '')
         overlay_stim(gca, cfg.stim_parameter, db_beh_RCSXX);   
     end
     format_plot();
            
%% MPQ
    if cfg.subplot == true

        subplot(313)
    else
        figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
        title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16); 
        hold on

    end
    
    MPQ_aff       = sum([redcap_RCSXX.MPQsickening, redcap_RCSXX.MPQfearful, redcap_RCSXX.MPQcruel],2,'omitnan');
    MPQ_som       = redcap_RCSXX.MPQtotal - MPQ_aff;

    if strcmp(cfg.dates, 'AllTime') == 1
    
        plot(redcap_RCSXX.time, movmean(redcap_RCSXX.MPQtotal , 5),'LineWidth', 4.5);
        hold on
        plot(redcap_RCSXX.time, movmean([MPQ_som, MPQ_aff], 5),'LineWidth', 2);

         set(gca,'ColorOrderIndex',1);

        plot(redcap_RCSXX.time, redcap_RCSXX.MPQtotal , '.', 'MarkerSize',9);
        plot(redcap_RCSXX.time, [MPQ_som, MPQ_aff], '.', 'MarkerSize',5);

    else
 
        scatter(redcap_RCSXX.time, redcap_RCSXX.MPQtotal , 100, 'filled');   
        hold on;
        scatter(redcap_RCSXX.time, MPQ_som, 50, 'filled');
        scatter(redcap_RCSXX.time, MPQ_aff, 50, 'filled');

    end
    
    ylabel('McGill Pain Questionaire');     ylim([0,45]); yticks(0:15:45)
    
    legend({ 'MPQ Total (0-45)','MPQ Somatic (0-33)', 'MPQ Affective (0-12)'}, ...
        'Location','northeastoutside'); 

   % overlay_stim(cfg.stim_parameter)
    format_plot();
    
     if ~strcmp(cfg.stim_parameter, '')
         overlay_stim(gca, cfg.stim_parameter, db_beh_RCSXX);   
     end
%% local functions

function overlay_stim(ax, stim, db_beh_RCSXX)
    if strcmp(stim, '')
        return


    elseif strcmp(stim, 'contacts')

        cont_pairs = unique(db_beh_RCSXX.stimRegOn(~strcmp(db_beh_RCSXX.stimRegOn,'')));


        c1 = brewermap(length(cont_pairs) ,'Set1');


        for i = 1: length(cont_pairs)           

            given_pair = strcmp(cont_pairs{i}, db_beh_RCSXX.stimRegOn);

            starts  = find(diff(given_pair)==1) + 1;
            stops   = find(diff(given_pair)== -1);

            % based of transitions plot every shading of a given
            % contact by itself

            if given_pair(end)

                stops = [stops; length(given_pair)];

            end

           
            for j = 1:length(starts)

                x = [db_beh_RCSXX.timeStart(starts(j)), db_beh_RCSXX.timeStop(stops(j)),...
                     db_beh_RCSXX.timeStop(stops(j)), db_beh_RCSXX.timeStart(starts(j))];

                y = [0,0,ax.YLim(2), ax.YLim(2)];

                if j == 1

                    patch(x, y, c1(i,:),'FaceAlpha', .3);

                    ax.Legend.String{end} = cont_pairs{i};

                else

                    patch(x, y, c1(i,:),'FaceAlpha', .3,...
                        'HandleVisibility','off');
                end

                hold on
            end
        end
    end
end


    %         elseif strcmp(stim, 'freq')
    %             
    %             cont_pairs = unique(db_RCSXXL.contacts);
    % 
    %             c = [colororder; cb_colors];
    % 
    %             for i = 2: length(cont_pairs) - 1            
    % 
    %                 i_given_pair = strcmp(cont_pairs{i}, db_RCSXXL.contacts) &...
    %                                    db_RCSXXL.amp ~=0;
    % 
    %                 starts = find(diff(i_given_pair)==1) + 1;
    %                 ends   = find(diff(i_given_pair)==-1);
    % 
    %                 % based of transitions plot every shading of a given
    %                 % contact by itself
    % 
    %                 if i_given_pair(end)
    % 
    %                     ends = [ends; length(i_given_pair)];
    % 
    %                 end
    % 
    %                
    %                 for j = 1:length(starts)
    % 
    %                     x = [db_RCSXXL.time(starts(j)), db_RCSXXL.time(ends(j)),...
    %                          db_RCSXXL.time(ends(j)), db_RCSXXL.time(starts(j))];
    %     
    %                     y = [0,0,ax.YLim(2), ax.YLim(2)];
    % 
    %                     if j == 1
    % 
    %                         patch(x, y, c(i,:));
    % 
    %                         ax.Legend.String{end} = cont_pairs{i};
    % 
    %                     else
    % 
    %                         patch(x, y, c(i,:),...
    %                             'HandleVisibility','off');
    %                     end
    % 
    %                     hold on
    %                 end
    %             end
            

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
