function plot_longitudinal_aDBS(cfg, REDcap, INS_logs_proc, app_SS_tbl, INS_ss_merge_g_changes)


% cfg = [];
% cfg.pt_id_side     = 'RCS02R';

% cfg.save_dir       = [github_dir, 'Analysis-rcs-data/working/plot_ephy/aDBS_offline_sessions/'];


proc_app        = INS_logs_proc.(cfg.pt_id_side).app;
app_ss_tbl      = app_SS_tbl.(cfg.pt_id_side);

proc_g_chan     = INS_ss_merge_g_changes.(cfg.pt_id_side);

redcap          = REDcap.(cfg.pt_id_side(1:end-1));

%%

if strcmp(cfg.dates, 'DateRange') == 1

    date_range  = datetime(cfg.date_range, 'TimeZone', 'America/Los_Angeles', 'InputFormat','dd-MMM-uuuu');
    
    i_entries   = find(ge(proc_app.time_INS, date_range(1)) & ...
                     le(proc_app.time_INS, date_range(2)));
    app_oi      = proc_app(i_entries,:);


%     i_entries   = find(ge(proc_group.time_INS, date_range(1)) & ...
%                      le(proc_group.time_INS, date_range(2)));
%     group_oi    = proc_group(i_entries,:);


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
%%
save_dir = [cfg.save_dir, cfg.pt_id_side ,'/'];

if ~isfolder(save_dir)

    mkdir(save_dir);

end
%%
u_settings  = unique(app_oi.sess_w_same_settings);



for i=   1:length(u_settings)
%%
    j              =  u_settings(i);
    i_ss           = find(ss_tbl_oi.sess_w_same_settings == j);
    
    plt_ss_tbl_oi  = ss_tbl_oi(i_ss(1), :);
    plt_app_oi     = app_oi(app_oi.sess_w_same_settings == j,...
                                            :);

% 01/26/23 RBL: entertain when unique setting has been tried on different
% weeks/months
%     i_aDBS_ends = [find(ge(t_diff, duration('6:00:00'))); height(plt_app_oi)];
% 
%     if ~isempty(i_aDBS_ends)
% 
%         i_aDBS_starts = [1; zeros(length(i_aDBS_ends)-1,1)];
% 
%         i_aDBS_starts(2:end) = i_aDBS_ends(1:end-1)+1;
% 
% 
%     end


    if ge(plt_app_oi.time_INS(end) -plt_app_oi.time_INS(1) , duration('00:10:00')) &&...
            le(plt_app_oi.time_INS(end) -plt_app_oi.time_INS(1) , duration('21:00:00:00'))

        step_dur        = duration('00:00:10');

        %t_in_state      = diff([NaT(1,1,'TimeZone', 'America/Los_Angeles'); plt_app_oi.time_INS]);
        
        
        start_time      = dateshift(plt_app_oi.time_INS(1), 'start', 'minute');
        stop_time       = dateshift(plt_app_oi.time_INS(end), 'end', 'minute');
        
        t_vec           = start_time:step_dur:stop_time;
        state_vec       = nan(length(t_vec),1);
        
        stim_vec        = state_vec;
        stream_sess_vec = state_vec;

     
        
   
        amp_def  = struct;

        [amp_def.u_states, i_u_states] = unique(plt_app_oi.newstate);

        amp_def.u_states(end+1) = 15;

        amp_def.prog0mA = [plt_app_oi{i_u_states, 'prog0mA'}; NaN];
        amp_def.prog1mA = [plt_app_oi{i_u_states, 'prog1mA'}; NaN];
        amp_def.prog2mA = [plt_app_oi{i_u_states, 'prog2mA'}; NaN];
        amp_def.prog3mA = [plt_app_oi{i_u_states, 'prog3mA'}; NaN];
  
        
        
        
        % from linearly-spaced vector from start to end time, define based
        % off of times of INS entries
        
        for h = 2:height(plt_app_oi)
        
            i_t_vec = find(ge(t_vec, plt_app_oi.time_INS(h-1)) &...
                           le(t_vec, plt_app_oi.time_INS(h)));
        
            state_vec(i_t_vec) = plt_app_oi.oldstate(h);

            % stim vector according to amplitude setting of OLD state NOT
            % current state
            stim_vec(i_t_vec)             =  amp_def.prog0mA(amp_def.u_states== plt_app_oi.oldstate(h));
 
            stream_sess_vec(i_t_vec)      = plt_app_oi.sess_w_same_settings(h);
        
        end

        % any amplitude of stim > 0 mA consider aDBS on
        on_off_vec                = 100*(stim_vec > 0);
 
    
        [sense_meta, by_pb_meta, ld0_meta, ld1_meta]...
            ...
            = parse_aDBS_params(...
            ...
        plt_ss_tbl_oi);
         


        %% parse through g_chan_oi to tmp_tbl to ONLY show therapy status
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

%         D_on           = g_chan_oi.time_align(i_D_on);
% 
%         for h=1:length(D_on)
% 
% 
% 
%         end
%
%         d_ON     = find(strcmp(g_chan_oi.event, 'GroupD_On'));
% 
%         ol_ON    = find(any(strcmp(g_chan_oi.event,...
%                                  {'GroupA_On', 'GroupB_On', 'GroupC_On'})));
% 
% 
%         i_D_on    = find(contains(g_chan_oi.event,'GroupD_On'));


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
            any_on(end+1) = stop_time; %#ok<AGROW> 
        end

       
        wash_out         = table();
        wash_out.start   = any_off;           wash_out.stop    = any_on;


        %%% when offline aDBS is longer than 12 hours, view w/n 12 hour
        % chunks w/n folder

        longest_dur = duration('24:00:00');

        if ge(stop_time - start_time, longest_dur)

            t_sub_sess = start_time:longest_dur:stop_time;

            n_sub_sess = length(t_sub_sess) - 1;

        else

            n_sub_sess = 1;
            t_sub_sess = [start_time, stop_time];

        end

        offline_sess_name = sprintf('%s-->%s', ...
                                datestr(dateshift(start_time, 'start', 'day'), 'YY-mm-DD'),...
                                datestr(dateshift(stop_time, 'start', 'day'), 'YY-mm-DD'));

%%
        colors = brewermap(5 ,'Set1');
    
        for h = 1:  n_sub_sess
            figure('Units', 'Inches', 'Position', [0, 0, 18, 10]);
            
            fig_h = subplot(1,1,1);        fig_h.Position(4) = 0.6;
    
            stairs(t_vec, movmean(on_off_vec, [duration('00:05:00')/step_dur,0]));  
            
            hold on

            plot(t_vec,   movmean(on_off_vec, [duration('01:00:00')/step_dur,0]), 'k', 'LineWidth',2);
            
            
            for k = 1:height(wash_out)  
                patch([wash_out.start(k), wash_out.start(k), wash_out.stop(k), wash_out.stop(k)],...
                           [0,100,100,0],[0.7, 0.7,0.7], ...
                    'FaceAlpha',0.5,'EdgeColor', 'none');  hold on
            end

            if ge(plt_app_oi.time_INS(end) -plt_app_oi.time_INS(1), duration('03:00:00'))      
                round_to = 'hour';
            else
                round_to = 'minute';
            end
    

            time_ticks = dateshift(...
               linspace(t_sub_sess(h), t_sub_sess(h+1), 8),...
                    "start", round_to);

            xticks(time_ticks);       grid on; grid minor;

    
    %{
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
    
            for i_off = 1:size(plt_ther_off,2)
    
                patch(sort([plt_ther_off(:,i_off);plt_ther_off(:,i_off)]),...
                           [0,100,100,0],[0.7, 0.7,0.7], ...
                    'FaceAlpha',0.5,'EdgeColor', 'none'); hold on
    
            end
    
    %}
            ylabel(['Percent time ON', newline,'(stim amplitude > 0 mA)'], 'FontSize',18)
            
            ylim([0,100]);        
            
            yyaxis right;          
        
            fig_h.YAxis(2).Color = colors(1,:);

            scatter(redcap.time, redcap.mayoNRS,'filled',...
                'MarkerFaceAlpha', 0.6, ...
                'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:));  
        
            plot(redcap.time, movmean(redcap.mayoNRS, [3, 0], 'omitnan'), ...
                'LineStyle', '-','LineWidth', 2, 'Color', colors(1,:));
        

            ylabel('NRS intensity', 'FontSize',18);     ylim([0,10]);

            fig_h.Position([1,3]) = [0.075, 0.85];

            %%% plot 
            dur_range =  t_sub_sess(h+1) -  t_sub_sess(h);
        
            text(t_sub_sess(h) - dur_range/15, fig_h.YLim(2)*1.43, cfg.pt_id_side , 'FontSize', 32);
            
            text(t_sub_sess(h) - dur_range/15, fig_h.YLim(2)*1.2, sense_meta, 'FontSize', 12);
            text(t_sub_sess(h) + dur_range/4.5, fig_h.YLim(2)*1.22, by_pb_meta, 'FontSize',12, 'Interpreter','none');
            
            
            text(t_sub_sess(h+1) - dur_range/1.85, fig_h.YLim(2)*1.2, ld0_meta, 'FontSize',12, 'Interpreter','none');
            text(t_sub_sess(h+1) - dur_range/4.25, fig_h.YLim(2)*1.2, ld1_meta, 'FontSize',12, 'Interpreter','none');
            
            

            

            fig_h.YAxis(1).FontSize = 14;   fig_h.YAxis(2).FontSize = 14;
            fig_h.XAxis.FontSize    = 14;
    
            legend({'Lagging Mean of 5 min', 'Lagging Mean of 1 hour', 'DBS Off', '', ''}, 'FontSize',14);

            set(gca,'xlim', [t_sub_sess(h), t_sub_sess(h+1)],...
            'GridAlpha',0.4,'MinorGridAlpha',0.7, 'GridColor', 'k', 'MinorGridColor', 'k'); 
        

            if ge(stop_time - start_time, longest_dur)

                if ~isfolder([save_dir, offline_sess_name])

                    mkdir([save_dir, offline_sess_name])
                end

                filename =  sprintf('%s/%s-->%s.png', ...
                                offline_sess_name, string(t_sub_sess(h), 'yy-MM-dd (hh mm)'), datestr(t_sub_sess(h+1), 'yy-MM-dd (hh mm)'));
            else
                filename = sprintf('%s.png', offline_sess_name);

            end

            exportgraphics(gcf, [save_dir,filename]);
        end
    end
    %catch

%     figure;
%     plot(1:10);   text(5, fig_h.YLim(2)*1.4, cfg.pt_id_side , 'FontSize', 32);
% 
%     exportgraphics(gcf, [save_dir,filename]);
% 

    %end
end  
end
  