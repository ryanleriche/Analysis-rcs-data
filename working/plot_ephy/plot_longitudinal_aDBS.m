function plot_longitudinal_aDBS(cfg, REDcap, INS_logs_proc, app_SS_tbl, INS_ss_merge_g_changes)

% cfg             = [];
% cfg.dates       = 'AllTime';
% cfg.proc_dir    = [pia_dir, 'processed/aDBS_offline_sessions/'];
% 
% cfg.dates          = 'DateRange';
% cfg.date_range     = {'01-Jan-2023'; '30-May-2023'};
% cfg.pt_sides       = {'RCS02R'}; %, 'RCS05L', 'RCS05R'};
% 
% i =1;
% 

% see cfg.proc_dir for the aDBS longitudinal plots)
set(0,'DefaultFigureVisible','off')

for d = 1:length(cfg.pt_sides)
    % per pt side, pull processed INS logs, streaming sessions, and REDcap
    proc_app        = INS_logs_proc.(cfg.pt_sides{d}).app;
    app_ss_tbl      = app_SS_tbl.(cfg.pt_sides{d});
    
    proc_g_chan     = INS_ss_merge_g_changes.(cfg.pt_sides{d});
    
    redcap          = REDcap.(cfg.pt_sides{d}(1:end-1));
    
    % check for/make saving directory
    save_dir        = [cfg.proc_dir, cfg.pt_sides{d} ,'/'];
    if ~isfolder(save_dir)
        mkdir(save_dir);

    end

    % option to run specific dates
    if strcmp(cfg.dates, 'DateRange') == 1
    
        date_range  = datetime(cfg.date_range, 'TimeZone', 'America/Los_Angeles', 'InputFormat','dd-MMM-uuuu');
        
        i_entries   = find(ge(proc_app.time_INS, date_range(1)) & ...
                         le(proc_app.time_INS, date_range(2)));
        app_oi      = proc_app(i_entries,:);
    
    
        i_entries   = find(ge(app_ss_tbl.timeStart, date_range(1)) & ...
                         le(app_ss_tbl.timeStart, date_range(2)));
        ss_tbl_oi   = app_ss_tbl(i_entries,:);
    
    
        i_entries   = find(ge(proc_g_chan.time_align, date_range(1)) & ...
                         le(proc_g_chan.time_align, date_range(2)));
        g_chan_oi   = proc_g_chan(i_entries, :);
    
    
    elseif strcmp(cfg.dates, 'AllTime') == 1
    
        app_oi      = proc_app;
        ss_tbl_oi   = app_ss_tbl;
    
    end

 %% go through INS logs based on their unique settings
    u_settings  = unique(app_oi.sess_w_same_settings);
    
    for i =   1:length(u_settings)
   
        j              =  u_settings(i);
        i_ss           = find(ss_tbl_oi.sess_w_same_settings == j);
        
        plt_ss_tbl_oi  = ss_tbl_oi(i_ss(1), :);
        plt_app_oi     = app_oi(app_oi.sess_w_same_settings == j,...
                                                :);
    
        if ge(plt_app_oi.time_INS(end) -plt_app_oi.time_INS(1) , duration('00:10:00')) &&...
                le(plt_app_oi.time_INS(end) -plt_app_oi.time_INS(1) , duration('21:00:00:00'))
    
            step_dur        = duration('00:00:10');
            
            start_time      = dateshift(plt_app_oi.time_INS(1), 'start', 'minute');
            stop_time       = dateshift(plt_app_oi.time_INS(end), 'end', 'minute');
            
            t_vec           = start_time:step_dur:stop_time;
            state_vec       = nan(length(t_vec),1);
            
            stim_vec        = state_vec;
            stream_sess_vec = state_vec;
    
            % from linearly-spaced vector from start to end time, define based
            % off of times of INS entries
            
            for h = 2:height(plt_app_oi)
            
                i_t_vec = find(ge(t_vec, plt_app_oi.time_INS(h-1))...
                               &...
                               le(t_vec, plt_app_oi.time_INS(h))...
                               );
            
                state_vec(i_t_vec) = plt_app_oi.oldstate(h);
    
                if plt_app_oi.oldstate(h) == 15
                         stim_vec(i_t_vec)       = NaN;
                else
                    % stim vector according to amplitude setting of OLD state NOT current state
                    stim_vec(i_t_vec)       = plt_app_oi.prog0mA(h);
                end

                stream_sess_vec(i_t_vec)    = plt_app_oi.sess_w_same_settings(h);
            
            end
    
            % any amplitude of stim > 0 mA consider aDBS on
            on_off_vec                = 100*(stim_vec > 0);

            % attempt to define washout based on therapyStatus of
            % AppLog.txt

            %{
            on_off_vec(isnan(stim_vec)) = NaN;

            % therapyStatus Off based on AppLog.txt files only
            i_nan        = isnan(on_off_vec);
            
            i_starts_nan  = find(diff(i_nan)==1);
            i_ends_nan   = find(diff(i_nan)==-1);
            
            
            if i_nan(1)
            i_starts_nan = [1; i_starts_nan]; %#ok<AGROW> 
            end
            
            if i_nan(end)
            i_ends_nan = [i_ends_nan; length(on_off_vec)]; %#ok<AGROW> 
            end
            
            
            plt_ther_off = [t_vec(i_starts_nan); t_vec(i_ends_nan)];
            %}

            [sense_meta, by_pb_meta, ld0_meta, ld1_meta, state_meta]...
                ...
                = parse_aDBS_params(...
                ...
            plt_ss_tbl_oi);


            % therapy status based on group changes and streaming sessions
           
            g_chan_oi.therapyStatus  = replace(cellfun(@(x) x(end-2:end), g_chan_oi.event,...
                                                    'UniformOutput', false),...
                                       '_','');
    
            tmp_tbl         = g_chan_oi(~contains(g_chan_oi.event, {'Lead', 'Transitioning'}),:);
    
            i_rep_events    = diff([inf; findgroups(tmp_tbl.therapyStatus)])~=0;
            tmp_tbl         = tmp_tbl(i_rep_events,:);
    
    
            i_any_on    = find(contains(tmp_tbl.therapyStatus,'On'));
            i_any_off   = find(strcmp(tmp_tbl.therapyStatus, 'Off'));
    
    
            any_off          = tmp_tbl.time_align(i_any_off);
            any_on           = tmp_tbl.time_align(i_any_on);
    
            if i_any_on(1) ==1
                any_on = any_on(2:end);
            end
    
            % when last parsed RCS database shows stim still being On
            if length(any_off) < length(any_on)
                any_on = any_on(2:end); 
            end
    
            % when last parsed RCS database shows stim still off--end patch at
            % plt end
            if length(any_off) > length(any_on)
                any_on(end+1) = stop_time;
            end
         
            wash_out         = table();
            wash_out.start   = any_off;           wash_out.stop    = any_on;

            %% when offline aDBS is longer than 12 hours, view w/n 12 hour
            % chunks w/n folder
    
            longest_dur = duration('24:00:00');
    
            if ge(stop_time - start_time, longest_dur)
    
                t_sub_sess = dateshift(start_time, 'start', 'day')...
                                :longest_dur:...
                             dateshift(stop_time, 'end', 'day');
    
                n_sub_sess = length(t_sub_sess) - 1;
    
            else
    
                n_sub_sess = 1;
                t_sub_sess = [dateshift(start_time, 'start', 'day'), stop_time];
    
            end
    
            offline_sess_name = sprintf('%s-->%s', ...
                                    datestr(dateshift(start_time, 'start', 'day'), 'YY-mm-DD'),...
                                    datestr(dateshift(stop_time, 'start', 'day'), 'YY-mm-DD'));
    
            colors =brewermap(NaN, 'Set2');

            for h = 1:  n_sub_sess
    
                fig_h = figure('Units', 'Inches', 'Position', [0, 0, 22, 10]);
                subplot(5,3, 1:5)
                
                stairs(t_vec, movmean(on_off_vec, [duration('00:01:00')/step_dur,0]));  
                
                hold on
    
                plot(t_vec,   movmean(on_off_vec, [duration('01:00:00')/step_dur,0]), 'k', 'LineWidth',2);
                
                % attempt to show washout from consensus from 
                % EventLog.txt and DeviceSettings.json files
                
                for k = 1:height(wash_out)
                    % only include handle of inital patch for ease of plotting
                    % legend (legend account for EVERY object handle)
                    if k == 1
                        patch([wash_out.start(k), wash_out.start(k), wash_out.stop(k), wash_out.stop(k)],...
                               [0,100,100,0],[0.7, 0.7,0.7], ...
                                'FaceAlpha',0.5,'EdgeColor', 'none');  hold on
    
                    else
                        patch([wash_out.start(k), wash_out.start(k), wash_out.stop(k), wash_out.stop(k)],...
                               [0,100,100,0],[0.7, 0.7,0.7], ...
                            'FaceAlpha',0.5,'EdgeColor', 'none',...
                            'HandleVisibility','off'); 
                    end
                end
                
                
    
                if ge(plt_app_oi.time_INS(end) -plt_app_oi.time_INS(1), duration('03:00:00'))      
                    round_to = 'hour';
                else
                    round_to = 'minute';
                end
        
    
%                 time_ticks = dateshift(...
%                    linspace(t_sub_sess(h), t_sub_sess(h+1), 8),...
%                         "start", round_to);
    
              
    
        
                % attempt to show AppLog.txt files
                %{
                for i_off = 1:size(plt_ther_off,2)

                    if i_off == 1
        
                    patch(sort([plt_ther_off(:,i_off);plt_ther_off(:,i_off)]),...
                               [0,100,100,0],[0.7, 0.7,0.7], ...
                        'FaceAlpha',0.5,'EdgeColor', 'none'); hold on
                    else
                         patch(sort([plt_ther_off(:,i_off);plt_ther_off(:,i_off)]),...
                               [0,100,100,0],[0.7, 0.7,0.7], ...
                        'FaceAlpha',0.5,'EdgeColor', 'none', 'HandleVisibility','off'); hold on
                    end
        
                end
                %}
        
        
                ylabel(['Percent time ON', newline,'(stim amplitude > 0 mA)'], 'FontSize',18)
                
                ylim([-5,105]);        
                
                yyaxis right; 

                plt_rcap = redcap;
                plt_rcap.painVAS  = redcap.painVAS/10;
                plt_rcap.MPQtotal = redcap.MPQtotal/4.5;
           
                pain_metrics =  {'mayoNRS', 'painVAS', 'MPQtotal'};

                for i_pain =1:size(pain_metrics,2)
                    plot(plt_rcap.time, movmean(plt_rcap.(pain_metrics{:,i_pain}), 3),...
                        'LineWidth', 3, 'Color', colors(i_pain+3,:),...
                        ...
                        'Marker', 'none','LineStyle', '-',...
                        'HandleVisibility','off'); hold on

                    scatter(plt_rcap.time, plt_rcap.(pain_metrics{:,i_pain}), 130, ...
                        'filled','MarkerFaceAlpha', 0.6, ...
                        'MarkerEdgeColor', colors(i_pain+3,:), 'MarkerFaceColor', colors(i_pain+3,:)...
                        );   
                end


                y_lbls_txt = cellfun(@(x) sprintf('%g   %g   %g',x), ...
                num2cell([0:10; 0:10:100; 0:4.5:45]', 2), 'UniformOutput', false);
                
                yticks(0:10); ylim([0,10]);  yticklabels(y_lbls_txt);
                
                fig_h.CurrentAxes.YAxis(1).FontSize = 14;     
                fig_h.CurrentAxes.YAxis(2).FontSize = 14; 
                
                fig_h.CurrentAxes.XAxis.FontSize    = 14;
                
                legend({'Lagging Mean of 1 min', 'Lagging Mean of 1 hour', 'aDBS Off', ...
                    'NRS intensity', 'VAS intensity', 'MPQ total'},...
                    'FontSize',14, 'Location', 'northoutside', 'NumColumns', 6);

                grid on;

                set(gca,'xlim', [t_sub_sess(h), t_sub_sess(h+1)],...
                        'GridAlpha',0.3,...
                        'YColor', 'k', 'TickLength', [0,0]);

                xticks(t_sub_sess(h):duration('1:00:00'):t_sub_sess(h+1));

                % reduce axis to 65% of its size so sense, and LD data can fit
                % above
                fig_h.CurrentAxes.Position(2) = fig_h.CurrentAxes.Position(2) *.65;
    
   
                % based on time range shown, add pt, sense, and LD meta
                % data
                dur_range =  t_sub_sess(h+1) -  t_sub_sess(h);
            
                text(t_sub_sess(h) - dur_range/15, fig_h.CurrentAxes.YLim(2)*1.9, cfg.pt_sides{d} , 'FontSize', 32);
                
                text(t_sub_sess(h) - dur_range/15, fig_h.CurrentAxes.YLim(2)*1.4, sense_meta, 'FontSize', 10);
                text(t_sub_sess(h) + dur_range/5.5, fig_h.CurrentAxes.YLim(2)*1.62, by_pb_meta, 'FontSize',10, 'Interpreter','none');
                
                
                text(t_sub_sess(h+1) - dur_range/1.6, fig_h.CurrentAxes.YLim(2)*1.6, ld0_meta, 'FontSize',10, 'Interpreter','none');
                text(t_sub_sess(h+1) - dur_range/2.5, fig_h.CurrentAxes.YLim(2)*1.6, ld1_meta, 'FontSize',10, 'Interpreter','none');
        
                text(t_sub_sess(h+1) - dur_range/8,...
                          fig_h.CurrentAxes.YLim(2) *1.6, ...
                          ...
                          [state_meta, sprintf('    GroupDrateInHz | %g', mode(plt_app_oi.rateHz(plt_app_oi.prog0mA > 0)))],...
                          ...
                          'FontSize',10, 'Interpreter','none');
             
                %% explicilty show state (0-7) to stim (current in mA) relationship
                plt_app_oi.oldstate(plt_app_oi.oldstate == 15) = NaN;

                subplot(5,3, 12:15)
        
                        stairs(plt_app_oi.time_INS(1:end-1), plt_app_oi.prog0mA(2:end),...
                            '-','LineWidth',1.25, 'Color',  'k');    hold on
        
                        yticks(0:0.5:3);
    
                        ylim([0, 3.25]); 
    
                        ylabel('Current (mA)');
    
                        set(gca, 'FontSize', 12);  
                        
                        yyaxis right
            
                            stairs(plt_app_oi.time_INS(1:end-1), plt_app_oi.oldstate(2:end), ...
                                '-','LineWidth',1.5, 'Color',  colors(1, :));  
                            ylim([0,8.5]); yticks(0:8); ylabel('State');    grid on;

                switch cfg.plt_state_dur
                    case 'two_chunks'
                        % specific sizing
                        fig_h.CurrentAxes.Position([1, 2, 3,4]) = [0.125, 0.1, .35, 0.24];
                        
                        set(gca, 'TickLength', [0, 0],'FontSize', 12, 'YColor', colors(1, :), ...
                                      'XLim', [t_sub_sess(h) + duration('1:00:00'), ...
                                      t_sub_sess(h) + + duration('2:00:00')],...
                                      'GridAlpha',0.3);
    
                        xticks(t_sub_sess(h) + duration('1:00:00')...
                            :duration('0:05:00'):...
                            t_sub_sess(h)+ duration('2:00:00'))

                        % repeatl ↑↑↑ code for second sub-session hour chunk
                        subplot(5,3, 15)

                        fig_h.CurrentAxes.Position([1, 2, 3,4]) = [0.557, 0.1, .35, 0.24];
        
                            stairs(plt_app_oi.time_INS(1:end-1), plt_app_oi.prog0mA(2:end),...
                                '-','LineWidth',1.75, 'Color',  'k');    hold on
            
                             yticks(0:0.5:3);
    
                             ylim([0, 3.25]);
        
                            ylabel('Current (mA)');
        
                            set(gca, 'FontSize', 12);  
                            
                            yyaxis right
                
                                stairs(plt_app_oi.time_INS(1:end-1), plt_app_oi.oldstate(2:end), ...
                                    '-','LineWidth',2, 'Color',  colors(1, :));  
                                ylim([0,8.5]); yticks(0:8); ylabel('State');    grid on;

                        set(gca, 'TickLength', [0, 0],'FontSize', 12, 'YColor', colors(1, :), ...
                                'XLim', [t_sub_sess(h) + duration('13:00:00'), ...
                                t_sub_sess(h) + + duration('14:00:00')],...
                                'GridAlpha',0.3);

                        xticks(t_sub_sess(h) + duration('13:00:00')...
                            :duration('0:05:00'):...
                            t_sub_sess(h)+ duration('14:00:00'))

                    case 'sub_session_duration'

                      set(gca,'TickLength', [0, 0],'FontSize', 12, 'YColor', colors(1, :), ...
                              'XLim', [t_sub_sess(h), t_sub_sess(h+1)],...
                              'GridAlpha',0.3);

                        % specific sizing
                        fig_h.CurrentAxes.Position([2, 4]) = [0.1, 0.24];

                        xticks(t_sub_sess(h):duration('0:30:00'):t_sub_sess(h+1)) 
                end
  
               %% make folders based on sub-session times, and save as .pngs
                if ge(stop_time - start_time, longest_dur)
    
                    if ~isfolder([save_dir, offline_sess_name])
    
                        mkdir([save_dir, offline_sess_name])
                    end
    
                    filename =  sprintf('%s/%s-->%s.png', ...
                                    offline_sess_name, string(t_sub_sess(h), 'yy-MM-dd'), string(t_sub_sess(h+1), 'yy-MM-dd'));
                else
                    filename = sprintf('%s.png', offline_sess_name);
    
                end
                % save_dir defined at start from cfg.proc_dir specified in
                % wrapper
                exportgraphics(gcf, [save_dir,filename]);
            end
        end
    end
    close all % remove hidden figures per pt to reduce overhead
end
set(0,'DefaultFigureVisible','on')
end