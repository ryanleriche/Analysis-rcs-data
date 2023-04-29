function  long_DBS_tbl = plot_timeline_DBS(cfg, pt_id, REDcap, INS_logs, app_SS_tbl)

% cfg             = [];
% cfg.dates       = 'AllTime';
% cfg.ephy_anal_dir    = [pia_dir, 'processed/aDBS_offline_sessions/'];
% 
% cfg.dates          = 'DateRange';
% cfg.date_range     = {'01-Jan-2023'; '30-May-2023'};
% cfg.pt_sides       = {'RCS02R'}; %, 'RCS05L', 'RCS05R'};
% 
% i =1;
% 

% see cfg.ephy_anal_dir for the aDBS longitudinal plots)
set(0,'DefaultFigureVisible','off')


redcap          = REDcap.(pt_id);

% per pt side, pull processed INS logs, and streaming sessions
all_pt_hemi     = fieldnames(INS_logs);
sides           = cellfun(@(x) {x(end)}, all_pt_hemi);

i_pt_hemi       = contains(all_pt_hemi, pt_id);

i_pt_R          = i_pt_hemi & strcmp(sides, 'R');
i_pt_L          = i_pt_hemi & strcmp(sides, 'L');

if any(i_pt_L)

    [L.dbs_oi, L.ss_tbl_oi] =...
    ...
    pull_INSLog_SS_oi(...
    ...
    cfg, all_pt_hemi{i_pt_L}, INS_logs, app_SS_tbl);

end

if any(i_pt_R)

    [R.dbs_oi, R.ss_tbl_oi] =...
    ...
    pull_INSLog_SS_oi(...
    ...
    cfg, all_pt_hemi{i_pt_R}, INS_logs, app_SS_tbl);

end

%%
% check for/make saving directory
save_dir        = [cfg.ephy_anal_dir, pt_id ,'/'];

if ~isfolder(save_dir);     mkdir(save_dir);           end


%% go through INS logs based on their unique settings
% initalize longitudinal DBS summary table (long_DBS_tbl)

long_DBS_tbl  = table;
u_settings    = unique(dbs_oi.sess_w_same_settings);

for i =   1   :  length(u_settings)

    j              =  u_settings(i);
    i_ss           = find(ss_tbl_oi.sess_w_same_settings == j);
    
    plt_ss_tbl_oi  = ss_tbl_oi(i_ss(1), :);
    plt_dbs_oi     = dbs_oi(dbs_oi.sess_w_same_settings == j,...
                                            :);

    
    if ge(plt_dbs_oi.time_INS(end) -plt_dbs_oi.time_INS(1) , duration('00:10:00')) &&...
            le(plt_dbs_oi.time_INS(end) -plt_dbs_oi.time_INS(1) , duration('21:00:00:00'))

      
        step_dur        = duration('00:00:10');
               
        start_time      = plt_dbs_oi.time_INS(1);
        stop_time       = plt_dbs_oi.time_INS(end);
        

        t_vec           = start_time:step_dur:stop_time;
        state_vec       = nan(length(t_vec),1);
        
        amp_vec         = state_vec;
        rate_vec        = state_vec;
        stream_sess_vec = state_vec;

      
        % from linearly-spaced vector from start to end time, define based
        % off of times of INS entries
        
        for h = 2:height(plt_dbs_oi)
        
            i_t_vec = find(ge(t_vec, plt_dbs_oi.time_INS(h-1))...
                           &...
                           le(t_vec, plt_dbs_oi.time_INS(h))...
                           );
        
            state_vec(i_t_vec) = plt_dbs_oi.oldstate(h);

            if plt_dbs_oi.oldstate(h) == 15
                     amp_vec(i_t_vec)       = NaN;
                     rate_vec(i_t_vec)      = NaN;
            else
                % stim vector according to amplitude setting of OLD state NOT current state
                amp_vec(i_t_vec)       = plt_dbs_oi.prog0mA(h);
                rate_vec(i_t_vec)      = plt_dbs_oi.rateHz(h);
            end

            stream_sess_vec(i_t_vec)    = plt_dbs_oi.sess_w_same_settings(h);
        
        end

        % any amplitude of stim > 0 mA consider aDBS on
        on_off_vec                = 100*(amp_vec> 0);

        [sense_meta, by_ld0_pb_meta, by_ld1_pb_meta,...
            ld0_meta, ld1_meta, state_meta]...
            ...
            = parse_aDBS_params(...
            ...
        plt_ss_tbl_oi);


        %% therapy status based on group changes and streaming sessions
        tmp_tbl         = g_chan_oi(~contains(g_chan_oi.event, {'Group','Lead', 'Transitioning'}),:);

        i_rep_events    = diff([inf; findgroups(tmp_tbl.event)])~=0;
        tmp_tbl         = tmp_tbl(i_rep_events,:);


        i_any_on    = find(contains(tmp_tbl.event,'On'));
        i_any_off   = find(strcmp(tmp_tbl.event, 'Off'));


        any_off          = tmp_tbl.time_INS(i_any_off);
        any_on           = tmp_tbl.time_INS(i_any_on);

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
            any_on(end+1) = stop_time; %#ok<*AGROW> 
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

        brew_col = brewermap(NaN, 'Set2');
%%% save times, streaming sessions, same session index, and "spontaneous"
%%% aDBS--w/ 0.1 Hz "samping rate" from the raw AppLog.txt outputs

        long_DBS_tbl.timeStart_INS_log(i)    = start_time;
        long_DBS_tbl.timeStop_INS_log(i)     = stop_time;
        
        long_DBS_tbl.sess_name{i}            = ss_tbl_oi.sess_name(i_ss);
        long_DBS_tbl.sess_w_same_settings(i) = u_settings(i);
        
        long_DBS_tbl.avg_percent_on(i)       = mean(on_off_vec, 'omitnan');

        if mean(on_off_vec, 'omitnan') > 0 || mean(on_off_vec, 'omitnan') ==100

            long_DBS_tbl.t_vec{i}                = t_vec;
            long_DBS_tbl.on_off_vec{i}           = on_off_vec;
            long_DBS_tbl.state_vec{i}            = state_vec;
            long_DBS_tbl.amp_vec{i}              = amp_vec;
            long_DBS_tbl.ratevec{i}              = rate_vec;
        end
%% plot day by day and whole time period with given--uninterrupted setting
        for h = 0:  n_sub_sess
            if h ==0
                t_plt_start = t_sub_sess(1);  t_plt_end   = t_sub_sess(end);

                x_tick_dur = duration('12:00:00');
                x_tick_rotation = 0;
            else
                t_plt_start = t_sub_sess(h);   t_plt_end   = t_sub_sess(h+1);

                x_tick_dur = duration('0:30:00');  
                x_tick_rotation = 0;
            end
             
            if mean(on_off_vec, 'omitnan') > 0 || mean(on_off_vec, 'omitnan') ==100
                fig_h = figure('Units', 'Inches', 'Position', [0, 0, 22, 10]);
                subplot(5,3, 1:5)
                
                plt_percent_ON_overtime(t_vec, on_off_vec, step_dur, ...
                                        wash_out, redcap, brew_col,...
                                        t_plt_start, t_plt_end)
    
                xticks(t_plt_start:x_tick_dur:t_plt_end);
    
    
                set(gca,'xlim', [t_plt_start, t_plt_end],...
                        'GridAlpha',0.3,...
                        'YColor', 'k', 'TickLength', [0,0], 'FontSize', 12);
    
                
                % reduce axis to 65% of its size so sense, and LD data can fit
                % above
                fig_h.CurrentAxes.Position(2)        = fig_h.CurrentAxes.Position(2) *.65;
                
                % format xtick labels to fit on multiple lines
                if h == 0 && ge(t_plt_end-t_plt_start, duration('72:00:00'))
                    fig_h.CurrentAxes.XTickLabelRotation = x_tick_rotation;
                    tmp_tick = cellfun(@(x) split(x,','), fig_h.CurrentAxes.XTickLabel, 'UniformOutput', false);
                    new_tick = cellfun(@(x) sprintf('%s\\newline%s\n', x{1},x{2}), tmp_tick, 'UniformOutput', false);
                    xticklabels([new_tick{:}]);
                end
    
    
                % based on time range shown, add pt, sense, and LD meta
                % data
                dur_range =  t_plt_end -  t_plt_start;
            
                text(t_plt_start - dur_range/15, fig_h.CurrentAxes.YLim(2)*1.9, pt_id , 'FontSize', 32);
                
                text(t_plt_start - dur_range/15, fig_h.CurrentAxes.YLim(2)*1.4, sense_meta, 'FontSize', 10);
    
                text(t_plt_start + dur_range/5.5, fig_h.CurrentAxes.YLim(2)*1.8, by_ld0_pb_meta, 'FontSize',10, 'Interpreter','none');
                text(t_plt_start + dur_range/5.5, fig_h.CurrentAxes.YLim(2)*1.45, by_ld1_pb_meta, 'FontSize',10, 'Interpreter','none');
                
                
                
                text(t_plt_end - dur_range/1.6, fig_h.CurrentAxes.YLim(2)*1.7, ld0_meta, 'FontSize',9, 'Interpreter','none');
                text(t_plt_end - dur_range/2.5, fig_h.CurrentAxes.YLim(2)*1.7, ld1_meta, 'FontSize',9, 'Interpreter','none');
        
                text(t_plt_end - dur_range/8,...
                          fig_h.CurrentAxes.YLim(2) *1.6, ...
                          ...
                          [state_meta, sprintf('    GroupDrateInHz | %g', mode(plt_dbs_oi.rateHz(plt_dbs_oi.prog0mA > 0)))],...
                          ...
                          'FontSize',10, 'Interpreter','none');
             
                %%% explicilty show state (0-7) to stim (current in mA) relationship
                plt_dbs_oi.oldstate(plt_dbs_oi.oldstate == 15) = NaN;
    
                subplot(5,3, 12:15)
        
                        stairs(plt_dbs_oi.time_INS(1:end-1), plt_dbs_oi.prog0mA(2:end),...
                            '-','LineWidth',1.25, 'Color',  'k');    hold on
        
                        yticks(0:0.5:3);            ylim([0, 3.25]); 
    
                        ylabel('Current (mA)');     set(gca, 'FontSize', 12);  
                        
                        yyaxis right
            
                            stairs(plt_dbs_oi.time_INS(1:end-1), plt_dbs_oi.oldstate(2:end), ...
                                '-','LineWidth',1.5, 'Color',  brew_col(1, :));  
                            ylim([0,8.5]); yticks(0:8); ylabel('State');    grid on;
    
                switch cfg.plt_state_dur
                    case 'two_chunks'
                        % specific sizing
                        fig_h.CurrentAxes.Position([1, 2, 3,4]) = [0.125, 0.1, .35, 0.24];
                        
                        set(gca, 'TickLength', [0, 0],'FontSize', 12, 'YColor', brew_col(1, :), ...
                                      'XLim', [t_plt_start + duration('1:00:00'), ...
                                      t_plt_start + + duration('2:00:00')],...
                                      'GridAlpha',0.3);
    
                        xticks(t_plt_start + duration('1:00:00')...
                            :duration('0:05:00'):...
                            t_plt_start+ duration('2:00:00'))
    
                        % repeatl ↑↑↑ code for second sub-session hour chunk
                        subplot(5,3, 15)
    
                        fig_h.CurrentAxes.Position([1, 2, 3,4]) = [0.557, 0.1, .35, 0.24];
        
                            stairs(plt_dbs_oi.time_INS(1:end-1), plt_dbs_oi.prog0mA(2:end),...
                                '-','LineWidth',1.75, 'Color',  'k');    hold on
            
                             yticks(0:0.5:3);         ylim([0, 3.25]);
                             ylabel('Current (mA)');   set(gca, 'FontSize', 12);  
                            
                             yyaxis right
                
                                stairs(plt_dbs_oi.time_INS(1:end-1), plt_dbs_oi.oldstate(2:end), ...
                                    '-','LineWidth',2, 'Color',  brew_col(1, :));  
                                ylim([0,8.5]); yticks(0:8); ylabel('State');    grid on;
    
                        set(gca, 'TickLength', [0, 0],'FontSize', 12, 'YColor', brew_col(1, :), ...
                                'XLim', [t_plt_start + duration('13:00:00'), ...
                                t_plt_start + + duration('14:00:00')],...
                                'GridAlpha',0.3);
    
                        xticks(t_plt_start + duration('13:00:00')...
                            :duration('0:05:00'):...
                            t_plt_start+ duration('14:00:00'))
    
                    case 'sub_session_duration'
    
                      set(gca,'TickLength', [0, 0],'FontSize', 12, 'YColor', brew_col(1, :), ...
                              'XLim', [t_plt_start, t_plt_end],...
                              'GridAlpha',0.3);
    
                        % specific sizing
                        fig_h.CurrentAxes.Position([2, 4]) = [0.1, 0.24];
    
                        xticks(t_plt_start:x_tick_dur:t_plt_end) 
                end
    
                 % format xtick labels to fit on multiple lines
                if h == 0 && ge(t_plt_end-t_plt_start, duration('72:00:00'))
                    fig_h.CurrentAxes.XTickLabelRotation = x_tick_rotation;
                    tmp_tick = cellfun(@(x) split(x,','), fig_h.CurrentAxes.XTickLabel, 'UniformOutput', false);
                    new_tick = cellfun(@(x) sprintf('%s\\newline%s\n', x{1},x{2}), tmp_tick, 'UniformOutput', false);
                    xticklabels([new_tick{:}]);
                end
    
               %%% make folders based on sub-session times, and save as .pngs
                
                if h == 0 && ge(t_plt_end-t_plt_start, duration('72:00:00'))
                     
                    filename =  sprintf('%s (whole time-period).png', ...
                                    offline_sess_name);
    
                    long_DBS_tbl.report_dir{i}    =  [save_dir,filename];

                else

                    if ge(stop_time - start_time, longest_dur)
                    
                        if ~isfolder([save_dir, offline_sess_name]);  mkdir([save_dir, offline_sess_name]);    end
                    
                        filename =  sprintf('%s/%s-->%s.png', ...
                        offline_sess_name, string(t_plt_start, 'yy-MM-dd'), string(t_plt_end, 'yy-MM-dd'));
                    else
                        filename = sprintf('%s.png', offline_sess_name);
                    
                    end
                end
                % defined at start from cfg.ephy_anal_dir as specified in wrapper
                if isfile([save_dir,filename]);    delete([save_dir,filename]);   end
                exportgraphics(gcf, [save_dir,filename]);

            else % if aDBS was always or never On
                 long_DBS_tbl.report_dir{i}    = 'No report (0% or 100% duty cycle)';
            end
        end

    else % if aDBS is of uninteresting duration (<10 min | > 21 days)
        long_DBS_tbl.sess_w_same_settings(i) = NaN;
        long_DBS_tbl.avg_percent_on(i)       = NaN;
        long_DBS_tbl.report_dir{i}           = 'No report (not between 10 min or 21 days)';
    end
end

long_DBS_tbl= sortrows(long_DBS_tbl,'timeStart_INS_log','ascend');
close all % remove hidden figures per pt to reduce overhead
set(0,'DefaultFigureVisible','on')
end