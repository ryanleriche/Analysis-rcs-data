function plot_timeline(cfg, REDcap, varargin)
    

% cfg                     = [];
% cfg.pt_id               = 'RCS04';
% cfg.stage_dates         = stage_dates{4}; % starts at Stage 1
% cfg.subplot             = true;
% 
% cfg.stim_parameter      = 'all';
% 
% cfg.dates               = 'DateRange';
% cfg.date_range               = {'14-Jul-2022'; '6-Sep-2022'};
% cfg.subplot             = false;
% 
% redcap                  = wrt_stim_REDcap.RCS04;
%%

redcap                  = REDcap.(cfg.pt_id);
[~, stage_dates]        = make_visit_dates;

cfg.stage_dates         = stage_dates{str2double(cfg.pt_id(end))}; % starts at Stage 1

    figure('Units', 'Inches', 'Position', [0, 0, 15, 7])

    [redcap, date_range] = date_parser(cfg, redcap);
    
    ds =        datestr(date_range,'dd-mmm-yyyy');


    sum_stats =   calc_sum_stats(cfg, redcap);
    

% allows exploratory analysis w/ consistent nice formatting
switch cfg.pt_id(7:end)

    case 'mismatch'


        i_trim_VAS_all   = ~(redcap.unpleasantVAS == 50 & redcap.painVAS == 50 & redcap.worstVAS == 50)...
            &...
                ~(...
                (redcap.unpleasantVAS == 50 | redcap.painVAS == 50 | redcap.worstVAS == 50)...
                & redcap.mayoNRS ~= 5 ...
                ...
                & sum(isnan([redcap.unpleasantVAS, redcap.painVAS == 50, redcap.worstVAS]),2) == 0);
        
        
       redcap  = redcap(i_trim_VAS_all , :);

        title([cfg.pt_id(1:5), ': VAS to NRS Timeline',...
            newline,...
            ds(1,:) ' to ' ds(2,:),...
            newline, ...
            'Stages 1, 2, and 3'], 'Fontsize',16);

        yyaxis left

        plot(redcap.time, movmean(redcap.mayoNRS, 5), 'LineWidth', 2); hold on;
        plot(redcap.time, [redcap.mayoNRS], '.', 'MarkerSize',5);
   

        ylabel('Numeric Rating Scale');   ylim([0,10]);   yticks(1:1:10);

        yyaxis right

        plot(redcap.time, movmean(redcap.painVAS, 5), 'LineWidth', 2); hold on;
        plot(redcap.time, [redcap.painVAS], '.', 'MarkerSize',5);
    
        ylabel('Visual Analog Scale');   ylim([0,100]);  yticks(0:10:100);


        legend(''); 
     

        format_plot();

        hold off;


        exportgraphics(gcf,...
            [cd,'/plot_beh/figs/beh_only/RCS02/nrs_vas_timeline.png'],...
            'Resolution',300); 

        %% plot difference of NRS and VAS/10
        figure('Units', 'Inches', 'Position', [0, 0, 25, 10])


        title([cfg.pt_id(1:5), ': VAS, NRS difference Timeline',...
            newline,...
            ds(1,:) ' to ' ds(2,:),...
            newline, ...
            'Stages 1, 2, and 3'], 'Fontsize',16);

        hold on;
        plot(redcap.time, movmean((redcap.painVAS/10 - redcap.mayoNRS), 5), 'LineWidth', 2); 
        plot(redcap.time, (redcap.painVAS/10 - redcap.mayoNRS), '.', 'MarkerSize',8, 'Color','k');
   

        yline([-1, 1], 'LineWidth', 3) 
        ylabel('(VAS/10) - NRS');   ylim([-3.5,5.75]);  legend(''); 

        format_plot();

        exportgraphics(gcf,...
            [cd,'/plot_beh/figs/beh_only/RCS02/nrs_vas_diff_timeline.png'],...
            'Resolution',300); 
%%
    otherwise
%% NRS
    if cfg.subplot == true
        ax=subplot(311);
        sgtitle([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);

    else
        title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);
        hold on
    end


    if strcmp(cfg.dates, 'AllTime') == 1
    
        plot(redcap.time, movmean([redcap.mayoNRS,redcap.worstNRS], 5),...
            'LineWidth', 2);

        hold on; set(gca,'ColorOrderIndex',1);

        plot(redcap.time, [redcap.mayoNRS,redcap.worstNRS], '.', 'MarkerSize',5);
    

    else

        plot(redcap.time, movmean(redcap.mayoNRS, [2 0], 'omitnan'),...
            'LineWidth', 2); hold on;

        set(gca,'ColorOrderIndex',1)

        scatter(redcap.time, redcap.mayoNRS, 125, 'filled','MarkerFaceAlpha', 0.6);    
        
        if strcmp(cfg.pt_id, 'RCS06')
            set(gca,'ColorOrderIndex',3)

            plot(redcap.time, movmean(redcap.nocNRS, [2 0], 'omitnan'),...
                'LineWidth', 2); hold on;
                set(gca,'ColorOrderIndex',3)
            scatter(redcap.time, redcap.nocNRS, 75, 'filled','MarkerFaceAlpha', 0.6); 


            plot(redcap.time, movmean(redcap.npNRS, [2 0], 'omitnan'),...
                'LineWidth', 2); hold on;
                set(gca,'ColorOrderIndex',4)
            scatter(redcap.time, redcap.npNRS, 50, 'filled', 'MarkerFaceAlpha', 0.6);

            legend({'','NRS Intensity','','NRS Nociceptive',...
                    '','NRS Neuropathic'}, 'Location','northeastoutside');

        elseif strcmp(cfg.pt_id, 'RCS07')

            scatter(redcap.time, redcap.unpNRS, 50, 'filled','MarkerFaceAlpha', 0.6); 
            scatter(redcap.time, redcap.leftarmNRS, 65, 'filled', 'MarkerFaceAlpha', 0.6);

            scatter(redcap.time, redcap.leftlegNRS, 50 , 'filled', 'MarkerFaceAlpha', 0.6);  
            scatter(redcap.time, redcap.leftfaceNRS, 35, 'filled', 'MarkerFaceAlpha', 0.6);

            legend({'','NRS Intensity', '','NRS Unpleasantness',...
                    'NRS Left Arm', 'NRS Left Leg', 'NRS LeftFace'},...
                    'Location','northeastoutside');
     
        else
            scatter(redcap.time, redcap.worstNRS, 75, 'filled', 'MarkerFaceAlpha', 0.6);

        end


    end
  
     ylabel('Numeric Rating Scale');     ylim([0,10]); yticks(1:2:10);
         
     if ~strcmp(cfg.pt_id, {'RCS06', 'RCS07'})

        legend({'','NRS Intensity', 'NRS Worst Intensity'}, 'Location','northeastoutside'); 
     
     end
     
     hold on
     
     if cfg.sum_stat_txt
         sum_stat_str =  sprintf('mean: %0.2f ± %0.2f\nrange: %0.0f (max %0.0f – min %0.0f)',...
                                  sum_stats.mayoNRS({'mean','std', 'range', 'max', 'min'}));
         
         text(ax.XTick(round(length(ax.XTick)/3)), 8, sum_stat_str, 'Interpreter','none', 'FontSize', 10)
     end

     if ~strcmp(cfg.stim_parameter, '')
         overlay_stim(gca, cfg.stim_parameter, redcap);   
     end

     format_plot();
     
%% VAS
    if cfg.subplot == true
        ax=subplot(312);
    else
        figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
        title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16); 
        hold on

    end


    if strcmp(cfg.dates, 'AllTime') == 1
    
        plot(redcap.time, movmean([redcap.painVAS,redcap.unpleasantVAS,redcap.worstVAS], 5),...
            'LineWidth', 2);

             hold on; set(gca,'ColorOrderIndex',1);

        plot(redcap.time, [redcap.painVAS,redcap.unpleasantVAS,redcap.worstVAS], '.', 'MarkerSize',5);
    
    else

        scatter(redcap.time, redcap.painVAS, 100, 'filled', 'MarkerFaceAlpha', 0.6);
            hold on;

        if strcmp(cfg.pt_id, 'RCS06')

            set(gca,'ColorOrderIndex',3);
            plot(redcap.time, movmean(redcap.nocVAS, [2 0], 'omitnan'),...
            'LineWidth', 2.5); hold on;

            set(gca,'ColorOrderIndex',3);
            scatter(redcap.time, redcap.nocVAS, 75, 'filled', 'MarkerFaceAlpha', 0.6);

            set(gca,'ColorOrderIndex',4);
            plot(redcap.time, movmean(redcap.npVAS, [2 0], 'omitnan'),...
            'LineWidth', 1.5); hold on;

            set(gca,'ColorOrderIndex',4);
            scatter(redcap.time, redcap.npVAS, 50, 'filled', 'MarkerFaceAlpha', 0.6);


        else
            set(gca,'ColorOrderIndex',1);

             plot(redcap.time, movmean(redcap.painVAS, [2 0], 'omitnan'),...
            'LineWidth', 2.5); hold on;

            scatter(redcap.time, redcap.unpleasantVAS, 75, 'filled','MarkerFaceAlpha', 0.6);
            scatter(redcap.time, redcap.worstVAS, 50, 'filled','MarkerFaceAlpha', 0.6);

        end
    end

    if strcmp(cfg.pt_id, 'RCS07')

        
        plot(redcap.time, movmean(redcap.moodVAS, [2 0], 'omitnan'),...
            'LineWidth', 1.5); hold on;

        set(gca,'ColorOrderIndex',4);
        
        scatter(redcap.time, redcap.moodVAS, 75, 'filled','MarkerFaceAlpha', 0.6);
    end


     ylabel('Visual Analog Scale');         ylim([0,100]);  yticks(0:20:100);

     if strcmp(cfg.pt_id, 'RCS06')

        legend({'VAS Intensity','', 'VAS Nociceptive', '',...
            'VAS Neuropathic'}, ...
            'Location','northeastoutside');

     elseif strcmp(cfg.pt_id, 'RCS07')

        legend({'VAS Intensity', '','VAS Unpleasantness', '','', 'VAS Mood'},...
            'Location','northeastoutside');
     else

        legend({'VAS Intensity', '','VAS Unpleasantness', 'VAS Worst Intensity'}, ...
            'Location','northeastoutside'); 

     end

     if ~strcmp(cfg.stim_parameter, '')
         overlay_stim(gca, cfg.stim_parameter, redcap);   
     end

     if cfg.sum_stat_txt
        sum_stat_str =  sprintf('mean: %0.2f ± %0.2f\nrange: %0.0f (max %0.0f – min %0.0f)',...
                              sum_stats.painVAS({'mean','std', 'range', 'max', 'min'}));
     
        text(ax.XTick(round(length(ax.XTick)/3)), 80, sum_stat_str, 'Interpreter','none', 'FontSize', 10)
     end
     format_plot();
            
%% MPQ
    if cfg.subplot == true

        ax = subplot(313);
    else
        figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
        title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16); 
        hold on

    end
    
    MPQaff       = sum([redcap.MPQsickening, redcap.MPQfearful, redcap.MPQcruel, redcap.MPQtiring],2,'omitnan');
    MPQsom       = redcap.MPQtotal - MPQaff;

    if strcmp(cfg.dates, 'AllTime') == 1
    
        plot(redcap.time, movmean(redcap.MPQtotal , 5),'LineWidth', 4.5);
        hold on
        plot(redcap.time, movmean([MPQsom, MPQaff], 5),'LineWidth', 2);

         set(gca,'ColorOrderIndex',1);

        plot(redcap.time, redcap.MPQtotal , '.', 'MarkerSize',9);
        plot(redcap.time, [MPQsom, MPQaff], '.', 'MarkerSize',5);

    else

        plot(redcap.time, movmean(redcap.MPQtotal, [3 0], 'omitnan'),...
            'LineWidth', 2); hold on;

        set(gca,'ColorOrderIndex',1);
 
        scatter(redcap.time, redcap.MPQtotal , 75, 'filled','MarkerFaceAlpha', 0.6);   
        hold on;
        scatter(redcap.time, MPQsom, 50, 'filled', 'MarkerFaceAlpha', 0.6);
        scatter(redcap.time, MPQaff, 50, 'filled', 'MarkerFaceAlpha', 0.6);

    end
    
    ylabel('McGill Pain Questionaire');     ylim([0,45]); yticks(0:15:45)
    
    legend({ '','MPQ Total (0-45)','MPQ Somatic (0-33)', 'MPQ Affective (0-12)'}, ...
        'Location','northeastoutside'); 

    if cfg.sum_stat_txt

        
        sum_stat_str =  sprintf('mean: %0.2f ± %0.2f\nrange: %0.0f (max %0.0f – min %0.0f)',...
                              sum_stats.MPQtotal({'mean','std', 'range', 'max', 'min'}));
     
        text(ax.XTick(round(length(ax.XTick)/3)), 30, sum_stat_str, 'Interpreter','none', 'FontSize', 10)
    end

   % overlay_stim(cfg.stim_parameter)
    format_plot();
    
     if ~strcmp(cfg.stim_parameter, '')
         overlay_stim(gca, cfg.stim_parameter, redcap);   
     end
end
%% local functions
%{
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
            
%}
  function overlay_stim(ax, stim_params,  redcap)

    if strcmp(stim_params, '')
        return


    elseif strcmp(stim_params, 'all')

        redcap.L_cycleOnTime(isnan(redcap.L_cycleOnTime))    = 0;
        redcap.L_cycleOffTime(isnan(redcap.L_cycleOffTime))  = 0;

        redcap.L_i_ol_DBS = contains(redcap.L_activeGroup, {'A', 'B', 'C'});

        redcap.R_cycleOnTime(isnan(redcap.R_cycleOnTime))    = 0;
        redcap.R_cycleOffTime(isnan(redcap.R_cycleOffTime))  = 0;

        redcap.R_i_ol_DBS = contains(redcap.R_activeGroup, {'A', 'B', 'C'});
        

        redcap.R_stimfreq(redcap.R_stimfreq == 114.9) = 115;
        redcap.R_stimfreq(redcap.R_stimfreq == 115.7) = 115;


        redcap.L_stimfreq(redcap.L_stimfreq == 115.7) = 115;



        % Find the unique rows (which corresponds to unique slices of original array)
        compared_vars = {'R_stimContacts', 'R_stimAmp', 'R_stimfreq',...
            ...
            'L_stimContacts', 'L_stimAmp', 'L_stimfreq'};
        
        [~, i_redcap] = unique(redcap(:, compared_vars),'rows');

        
        u_stim = redcap(i_redcap, :);


        c1 = [brewermap(height(u_stim)-3 ,'Set1');brewermap(height(u_stim)-3 ,'Dark2')];


         for i = 1: height(u_stim)   

            [i_u_stim,~] = ismember(redcap(:, compared_vars),...
                               u_stim(i,compared_vars),'rows');

            starts  = find(diff(i_u_stim) ==  1) + 1;
            stops   = find(diff(i_u_stim) == -1) +1;

            % based of transitions plot every shading of a given
            % contact by itself

            if i == height(u_stim) || isempty(stops)

                stops = [stops; height(redcap)];

            end

            if isempty(starts)

                starts = 1;

            end
           
            for j = 1 : length(starts)

                x = [redcap.time(starts(j)), redcap.time(stops(j)),...
                     redcap.time(stops(j)), redcap.time(starts(j))];

                y = [0,0, ax.YLim(2), ax.YLim(2)];

                if j == 1

                    patch(x, y, c1(i,:),'FaceAlpha', .3);

                    str = '';

                    for k = 1 : length(compared_vars)

                        param = u_stim{i, compared_vars(k)};


                        if contains(compared_vars(k), 'R_stimContacts')
                    
                            str = [str,' ', param{1}];

                        elseif contains(compared_vars(k), 'R_stimAmp')
                    
                            str = [str,' ', num2str(param), 'mA'];

                        elseif contains(compared_vars(k), 'R_stimfreq')

                            str = [str,' ', num2str(param), 'Hz'];


                        elseif contains(compared_vars(k), 'L_stimContacts')

                            str = [str,' ', param{1}];

                        elseif contains(compared_vars(k), 'L_stimAmp')
                    
                            str = [str,' ', num2str(param), 'mA'];

                        elseif contains(compared_vars(k), 'L_stimfreq')

                            str = [str,' ', num2str(param), 'Hz'];

                        end                 
                    end
                

                    ax.Legend.String{end} = str;

                else

                    patch(x, y, c1(i,:),'FaceAlpha', .3,...
                        'HandleVisibility','off');
                end

                hold on
            end
         end

%{
    elseif strcmp(stim_params, 'contacts')

        cont_pairs = unique(redcap.stimRegOn(~strcmp(redcap.stimRegOn,'')));


        c1 = brewermap(length(cont_pairs) ,'Set1');


        for i = 1: length(cont_pairs)           

            given_pair = strcmp(cont_pairs{i}, redcap.stimRegOn);

            starts  = find(diff(given_pair)==1) + 1;
            stops   = find(diff(given_pair)== -1);

            % based of transitions plot every shading of a given
            % contact by itself

            if given_pair(end)

                stops = [stops; length(given_pair)];

            end

           
            for j = 1:length(starts)

                x = [redcap.timeStart(starts(j)), redcap.timeStop(stops(j)),...
                     redcap.timeStop(stops(j)), redcap.timeStart(starts(j))];

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
%}
    end
end

    
    
    
function format_plot()  

    grid on;    grid MINOR;   legend boxoff;    box off;

    set(gca,'FontSize',12, 'xlim', date_range , 'TickLength', [0 0],...
        'GridAlpha',0.4,'MinorGridAlpha',0.7, 'GridColor', 'k', 'MinorGridColor', 'k'); 

    t = datetime(cfg.stage_dates, 'TimeZone', '-07:00');

    if length(cfg.stage_dates) == 1
         xline(t, '-', {'Stage 1'},...
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
