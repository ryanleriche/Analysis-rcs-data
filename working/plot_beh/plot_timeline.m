function plot_timeline(cfg, redcap_RCSXX, db_RCSXXL, db_RCSXXR)

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
     overlay_stim(gca, cfg.stim_parameter, db_RCSXXL,db_RCSXXR);
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

     overlay_stim(gca, cfg.stim_parameter, db_RCSXXL, db_RCSXXR);
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
    overlay_stim(gca, cfg.stim_parameter, db_RCSXXL, db_RCSXXR);

%% local functions

    function overlay_stim(ax, stim, db_RCSXXL, db_RCSXXR)
        if strcmp(stim, '')
            db_RCSXXL = parse_db(db_RCSXXL);
    
            db_RCSXXR = parse_db(db_RCSXXR);
    
        elseif strcmp(stim, 'contacts')
                
    
                cont_pairs_L = unique(db_RCSXXL.contacts);
    
                cont_pairs_R = unique(db_RCSXXR.contacts);
    
                c1 = brewermap(length(cont_pairs_L) - 1 ,'Set3');
    
                c2 = brewermap(length(cont_pairs_R) - 1 ,'Set1');
    
    
                for i = 2: length(cont_pairs_L) - 1            
    
                    L_given_pair = strcmp(cont_pairs_L{i}, db_RCSXXL.contacts) &...
                                       db_RCSXXL.amp ~=0;
    
                    L_starts = find(diff(L_given_pair)==1) + 1;
                    L_stops   = find(diff(L_given_pair)== -1);
    
                    % based of transitions plot every shading of a given
                    % contact by itself
    
                    if L_given_pair(end)
    
                        L_stops = [L_stops; length(L_given_pair)];
    
                    end
    
                   
                    for j = 1:length(L_starts)
    
                        x = [db_RCSXXL.timeStart(L_starts(j)), db_RCSXXL.timeStop(L_stops(j)),...
                             db_RCSXXL.timeStop(L_stops(j)), db_RCSXXL.timeStart(L_starts(j))];
        
                        y = [0,0,ax.YLim(2), ax.YLim(2)];
    
                        if j == 1
    
                            patch(x, y, c1(i,:),'FaceAlpha', .3);
    
                            ax.Legend.String{end} = [db_RCSXXL.stimName{L_starts(j)}, ': ' cont_pairs_L{i}];
    
                        else
    
                            patch(x, y, c1(i,:),'FaceAlpha', .3,...
                                'HandleVisibility','off');
                        end
    
                        hold on
                    end
    
                    
                end
            
                for i = 2: length(cont_pairs_R) - 1            
    
                    R_given_pair = strcmp(cont_pairs_R{i}, db_RCSXXR.contacts) &...
                                       db_RCSXXR.amp ~=0;
    
                    R_starts = find(diff(R_given_pair)==1) + 1;
                    R_stops   = find(diff(R_given_pair)== -1);
    
                    % based of transitions plot every shading of a given
                    % contact by itself
    
                    if R_given_pair(end)
    
                        R_stops = [R_stops; length(R_given_pair)];
    
                    end
    
                   
                    for j = 1:length(R_starts)
    
                        x = [db_RCSXXR.timeStart(R_starts(j)), db_RCSXXR.timeStop(R_stops(j)),...
                             db_RCSXXR.timeStop(R_stops(j)), db_RCSXXR.timeStart(R_starts(j))];
        
                        y = [0,0,ax.YLim(2), ax.YLim(2)];
    
                        if j == 1
    
                            patch(x, y, c2(i,:),'FaceAlpha', .5);
    
                            ax.Legend.String{end} = [db_RCSXXR.stimName{R_starts(j)}, ': ' cont_pairs_R{i}];
    
                        else
    
                            patch(x, y, c2(i,:),'FaceAlpha', .5,...
                                'HandleVisibility','off');
                        end
    
                        hold on
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
