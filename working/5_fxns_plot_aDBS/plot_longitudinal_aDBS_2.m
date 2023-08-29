function  aDBS_sum = plot_longitudinal_aDBS_2(cfg, pt_side_id, REDcap, INS_logs, app_SS_tbl)

% see "cfg.anal_dir" for where aDBS longitudinal plots are saved
close all
set(0,'DefaultFigureVisible','off')
%set(0,'DefaultFigureVisible','on')


% per pt side, pull processed INS logs, streaming sessions, and REDcap
proc_app        = sortrows(INS_logs.(pt_side_id).app, 'time_INS');
app_ss_tbl      = app_SS_tbl.(pt_side_id);

proc_g_chan     = INS_logs.(pt_side_id).group_changes;
redcap          = REDcap.(pt_side_id(1:end-1));

% parse down to times of interst based of cfg.Dates field
[app_oi, ss_tbl_oi, g_chan_oi]...
    ...
    = aDBS_toi(...
    ...
cfg, proc_app, app_ss_tbl, proc_g_chan);
%%
% check for/make saving directory of aDBS plot .png files
save_dir        = fullfile(cfg.anal_dir,'aDBS_offline_sessions', pt_side_id);

if ~isfolder(save_dir);     mkdir(save_dir);           end

%%% return all previously made file/folder names
%{
* 'ignoreold' option to only plot new aDBS plots
    * easier than constantly changing the date range

* load previous 'aDBS_sum'

* see which streaming sessions have already been plotted, and thereby all
previous INS log entries since entries are pulled at streaming session start


* plot even when no new streaming session has occured 


*

%}

% %%
% cfg.proc_dir = [cfg.proc_dir, '/databases/'];
% 
% if ~exist(cfg.proc_dir, 'dir');    mkdir(cfg.proc_dir);      end
% 
% outputFileName    = fullfile(cfg.proc_dir,[pt_side_id '_database.mat']);
% 
% 
% if exist(outputFileName,'file') && ~cfg.ignoreold_db 
% 

aDBS_dir = fullfile(cfg.proc_dir, 'aDBS_summaries');

if ~isfolder(aDBS_dir);     mkdir(aDBS_dir);    end

aDBS_sum_path = fullfile(aDBS_dir, [pt_side_id,'_aDBS_sum.mat']);

%% go through INS logs based on their unique settings
% initalize longitudinal aDBS summary table (aDBS_sum)
aDBS_sum       = table;
u_settings     = unique(ss_tbl_oi.sess_setting_changes, 'stable');

for i =  1   : length(u_settings)

    i_ss           = find(ss_tbl_oi.sess_setting_changes == u_settings(i));
    
    % of the streaming sessions and AppLogs w/n user-specified date window,
    % return continuous settings/entries w/ no changes
    plt_ss_tbl_oi  = ss_tbl_oi(i_ss(1), :);
    plt_app_oi     = app_oi(...
                            plt_ss_tbl_oi.sess_setting_changes == app_oi.sess_setting_changes...
                            , :);


    if height(plt_app_oi) < 2;        continue;           end

    start_time      = plt_app_oi.time_INS(1);
    stop_time       = plt_app_oi.time_INS(end);


    % SAVE every unique setting
    aDBS_sum.timeStart_INS_log(i)    = start_time;
    aDBS_sum.timeStop_INS_log(i)     = stop_time;
    
    %%% use linearly-spaced vector from start to end time
    % --> interpolate aDBS state/stim for spontaneous duty cycle reporting
    step_dur        = duration('00:00:10');
    
    [t_vec, state_vec, amp_vec, rate_vec, ~, on_off_vec]...
        ...
        = spon_aDBS_status(...
        ...
     plt_app_oi, start_time, stop_time, step_dur);

    %%% from parsed database entry of streaming session (so-called plt_ss_tbl_oi)
    % return plot-ready character arrays
    meta_aDBS_txt   = parse_aDBS_params(plt_ss_tbl_oi);

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


    stop_time = ss_tbl_oi.timeStart(i_ss(end));


    if le(stop_time - start_time, longest_dur)

        n_sub_sess = 0;
        t_sub_sess = [dateshift(start_time, 'start', 'day'), stop_time];

    else 
        t_sub_sess = dateshift(start_time, 'start', 'day')...
                        :longest_dur:...
                     dateshift(stop_time, 'end', 'day');
    
        n_sub_sess = length(t_sub_sess) - 1;
    end

    offline_sess_name = sprintf('%s-->%s', ...
                            datestr(dateshift(t_sub_sess(1), 'start', 'day'), 'YY-mm-DD'),...
                            datestr(dateshift(t_sub_sess(end), 'start', 'day'), 'YY-mm-DD'));

   
    brew_col = brewermap(NaN, 'Set2');

    %%% save times, streaming sessions, same session index, and "spontaneous"
    %%% aDBS--w/ 0.1 Hz "samping rate" from the raw AppLog.txt outputs

    aDBS_sum.sess_name{i}             = ss_tbl_oi.sess_name(i_ss);
    aDBS_sum.sess_setting_changes(i)  = u_settings(i);
    
    aDBS_sum.avg_percent_on(i)        = mean(on_off_vec, 'omitnan');

    aDBS_sum.t_vec{i}                 = t_vec;
    aDBS_sum.on_off_vec{i}            = on_off_vec;
    aDBS_sum.state_vec{i}             = state_vec;
    aDBS_sum.amp_vec{i}               = amp_vec;
    aDBS_sum.ratevec{i}               = rate_vec;

%% plot day by day and whole time period with given--uninterrupted setting
    for h = 0 :  n_sub_sess
        
        if h == 0;      t_plt_start  = t_sub_sess(1);  t_plt_end   = t_sub_sess(end);
        else;           t_plt_start  = t_sub_sess(h);  t_plt_end   = t_sub_sess(h+1);        end
    
        % for latest data (i.e., last of the unique settings)
        % -> explictly show last streaming session, and current time in
        % 'whole time period' and '24'
        up_to_now_bool = eq(app_oi.time_INS(end), plt_app_oi.time_INS(end))...
                            && h==0;

        if up_to_now_bool

            t_plt_end = datetime('now','TimeZone','local');

            offline_sess_name = ...
                sprintf('%s-->%s', ...
                    datestr(dateshift(start_time, 'start', 'day'), 'YY-mm-DD'),...
                    datestr(dateshift(t_plt_end, 'start', 'day'), 'YY-mm-DD'));

        end

        if ge(t_plt_end-t_plt_start, duration('14:00:00:00'))
            x_tick_dur = duration('24:00:00');


        elseif ge(t_plt_end-t_plt_start, duration('7:00:00:00'))

            x_tick_dur = duration('24:00:00');

        elseif ge(t_plt_end-t_plt_start, duration('3:00:00:00'))

            x_tick_dur = duration('12:00:00');

        else
            x_tick_dur = duration('00:30:00');

        end

        fig_h = figure('Units', 'Inches', 'Position', [0, 0, 22, 10]);
        subplot(5,3, 1:5)


        if eq(app_oi.time_INS(end), plt_app_oi.time_INS(end))

           hand_x = xline(ss_tbl_oi.timeStop(end), '-', ...
                            sprintf('Last streaming\n%s\n%s',...
                                ss_tbl_oi.sess_name{end}, ss_tbl_oi.timeStop(end))...
                                );
           hold on

           hand_x.LineWidth = 1.5;
           hand_x.HandleVisibility = 'off';

           act          = table;
           act.states   = unique(plt_app_oi.oldstate);
           act.prog0mA  = act.states;
           act.rateHz   = act.states;


           for k = 1 : length(act.states)

              act.prog0mA(k) =  mode(plt_app_oi.prog0mA(act.states(k) == plt_app_oi.oldstate));
              act.rateHz(k)  =  mode(plt_app_oi.rateHz(act.states(k) == plt_app_oi.oldstate));

           end
          
           last             = struct;
           last.state       = plt_app_oi.newstate(end);
           last.prog0mA     = act.prog0mA(act.states == last.state);
           last.rateHz      = act.rateHz(act.states == last.state);


           hand_x = xline(app_oi.time_INS(end), '-', ...
                    sprintf('Last INS Log entry\n %s\nstate %g (%g mA, %g Hz)\n%s\n%s',...
                        app_oi.time_INS(end), last.state, last.prog0mA, last.rateHz)...
                        );

           hand_x.LineWidth = 2;
           hand_x.HandleVisibility = 'off';
           hand_x.Color = 'r';
           hand_x.LabelHorizontalAlignment= 'left';

            try
                cur_t_vec = t_vec(end):step_dur:ss_tbl_oi.timeStop(end);
                
                plt_t_vec      = [t_vec, cur_t_vec];
                plt_on_off_vec = [on_off_vec; ...
                                    repmat(last.prog0mA > 0, length(cur_t_vec),1)];
                
                curr_app            = plt_app_oi(end, :);
                curr_app.time_INS   = ss_tbl_oi.timeStart(end);
                curr_app.oldstate   = curr_app.newstate;
                curr_app.prog0mA    = last.prog0mA;
                curr_app.rateHz     = last.rateHz;
                
                plt_app_oi         = [plt_app_oi; curr_app; curr_app];
            catch

                plt_t_vec = t_vec;
                plt_on_off_vec = on_off_vec;
            end
        else

            plt_t_vec = t_vec;
            plt_on_off_vec = on_off_vec;
   
        end
        

        %%% plot lagging 1-hr percent aDBS ON (i.e., whenever current > 0 mA)
        plt_percent_ON_overtime(plt_t_vec, plt_on_off_vec, step_dur, ...
                                wash_out, redcap, brew_col,...
                                t_plt_start, t_plt_end, x_tick_dur, pt_side_id)


        %%% based on time range add pt, sense, and LD meta data
        aDBS_plot_txt(fig_h, t_plt_start, t_plt_end, pt_side_id, meta_aDBS_txt, plt_app_oi, h)
        
     
        %%% explicilty show state (0-8) to stim (current in mA) relationship
        subplot(5,3, 12:15)

            plt_state_mA_overtime(fig_h, plt_app_oi, t_plt_start, t_plt_end, x_tick_dur,h)


       %%% make folders based on sub-session times, and save as .pngs
       % filename of the entirety given, continuous and unique set of aDBS settings 
  
       if h == 0

            filename                  =  sprintf('%s (whole time-period).png', offline_sess_name);
            aDBS_sum.report_dir{i}    =  fullfile(save_dir,filename);

        % filename of individual days
        else

            if ge(stop_time - start_time, longest_dur)
                filename =  sprintf('%s/%s-->%s.png', ...
                                    offline_sess_name, ...
                                    string(t_plt_start, 'yy-MM-dd'), string(t_plt_end, 'yy-MM-dd'));
            else
                filename = sprintf('%s.png', offline_sess_name);
            
            end

            % make new folder of for offline aDBS session
            if ~isfolder(fullfile(save_dir, offline_sess_name))
                mkdir(fullfile(save_dir, offline_sess_name))
            end
       end
        % save as .png
        exportgraphics(gcf, fullfile(save_dir, filename));
    end
%
%
%
end

if ~isempty(aDBS_sum)
    aDBS_sum = sortrows(aDBS_sum,'timeStart_INS_log','ascend');
    save(aDBS_sum_path ,"aDBS_sum",'-v7.3')
 
end

close all % remove hidden figures per pt to reduce overhead
set(0,'DefaultFigureVisible','on')

%% local fxns


end