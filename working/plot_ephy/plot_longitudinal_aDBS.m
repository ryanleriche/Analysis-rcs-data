cfg             = [];
cfg.pt_id       = 'RCS02';

cfg.dates       = 'AllTime';
%cfg.date_range  = {'20-Dec-2022', '24-Jan-2023'};

cfg.save_dir    = [github_dir, 'Analysis-rcs-data/working/plot_ephy/aDBS_offline_sessions/'];

proc_app        = proc_app_log.RCS02R;
app_ss_tbl      = app_SS_tbl.RCS02R;


redcap          = REDcap.(cfg.pt_id);
[~, date_range] = date_parser(cfg, redcap);


%proc_group      =  proc_group_changes.RCS02R;
%%

if strcmp(cfg.dates, 'DateRange') == 1

    date_range  = datetime(cfg.date_range, 'TimeZone', 'America/Los_Angeles');
    
    i_entries   = find(ge(proc_app.time_INS, date_range(1)) & ...
                     le(proc_app.time_INS, date_range(2)));
    app_oi      = proc_app(i_entries,:);


%     i_entries   = find(ge(proc_group.time_INS, date_range(1)) & ...
%                      le(proc_group.time_INS, date_range(2)));
%     group_oi    = proc_group(i_entries,:);


    i_entries   = find(ge(app_ss_tbl.timeStart, date_range(1)) & ...
                     le(app_ss_tbl.timeStart, date_range(2)));
    ss_tbl_oi   = app_ss_tbl(i_entries,:);


elseif strcmp(cfg.dates, 'AllTime') == 1


    app_oi      = proc_app;
    ss_tbl_oi   = app_ss_tbl;


%     i_off_ss    = [1;...
%                   find(ge(diff(app_oi.time_INS), duration('24:00:00')));...
%                   height(app_oi)];

%     aDBS_by_offline_sess = cell(length(u_settings),1);
% 
%     for i=1:length(u_settings)
% 
%         aDBS_by_offline_sess{i}  = app_oi(app_oi.sess_w_same_settings == u_settings(i),...
%                                             :);
%     end
end
%%
save_dir = [cfg.save_dir, cfg.pt_id,'/'];

if ~isfolder(save_dir)

    mkdir(save_dir);

end

close all
set(0,'DefaultFigureVisible','off')
u_settings  = unique(app_oi.sess_w_same_settings);

for i= 1 : 12%length(u_settings)
    
    i_ss           = find(ss_tbl_oi.sess_w_same_settings == u_settings(i));
    
    plt_ss_tbl_oi  = ss_tbl_oi(i_ss(1), :);

    plt_app_oi     = app_oi(app_oi.sess_w_same_settings == u_settings(i),...
                                            :);

    t_diff = diff(plt_app_oi.time_INS);

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


    if ge(plt_app_oi.time_INS(end) -plt_app_oi.time_INS(1) , duration('00:10:00'))
        t_in_state   = diff([NaT(1,1,'TimeZone', 'America/Los_Angeles'); plt_app_oi.time_INS]);
        
        
        start_time   = dateshift(plt_app_oi.time_INS(1), 'end', 'minute');
        stop_time    = dateshift(plt_app_oi.time_INS(end), 'start', 'minute');
        
        t_vec        = start_time:duration('00:00:30'):stop_time;
        state_vec    = nan(length(t_vec),1);
        
        stim_vec     = state_vec;
        on_off_vec   = state_vec;
        
        stream_sess_vec = state_vec;
        
        
        for j = 1:height(plt_app_oi)-1
        
            i_t_vec = find(ge(t_vec, plt_app_oi.time_INS(j)) & le(t_vec, plt_app_oi.time_INS(j+1)));
        
            state_vec(i_t_vec) = plt_app_oi.oldstate(j);
        
            stim_vec(i_t_vec)  = plt_app_oi.prog0mA(j);
        
            on_off_vec(i_t_vec) = 100*(plt_app_oi.prog0mA(j) > 0);
        
            stream_sess_vec(i_t_vec)      = plt_app_oi.sess_w_same_settings(j);
        
        end
    %%%
    
     
        all_vars = plt_ss_tbl_oi.Properties.VariableNames;
        
        plt_vars = all_vars(...
                        contains(all_vars,...
                         {'TDsampleRates', 'bandFormationConfig',...
                        'fft_interval', 'fft_size', 'powerBin','chanOut', }));
        
       
        static_across_sess = {'TDsampleRates', 'fft_intervalInMilliseconds', 'fft_sizeInSamples', 'bandFormationConfig'};
        
        
        vars = table2cell(plt_ss_tbl_oi(1,  static_across_sess));
        
        sense_meta = sprintf(['TD samp rate   | %.0f Hz',newline,...
                             'FFT interval    | %.0f ms', newline,...
                             'FFT size        | %.0f samples', newline,...
                             'Bit Shift       | %s'],...
                            vars{:});
        
        
        
        by_chans_pb = {'Ch0_chanOut', 'Ch0_powerBinInHz0', 'Ch0_powerBinInHz1',...
                        'Ch1_chanOut', 'Ch1_powerBinInHz2', 'Ch1_powerBinInHz3',...
                        'Ch2_chanOut', 'Ch2_powerBinInHz4', 'Ch2_powerBinInHz5',...
                        'Ch3_chanOut', 'Ch3_powerBinInHz6', 'Ch3_powerBinInHz7'};
        
        vars          = table2cell(plt_ss_tbl_oi(1,  by_chans_pb));
        
        by_ch_pb_meta = sprintf(['Ch0 | %s',newline,...
                                 '              PowerBin0 | %s', newline,...
                                 '              PowerBin1 | %s', newline,...
                                 'Ch1 | %s',newline,...
                                 '              PowerBin2 | %s', newline,...
                                 '              PowerBin3 | %s', newline,...
                                 'Ch2 | %s',newline,...
                                 '              PowerBin4 | %s', newline,...
                                 '              PowerBin5 | %s', newline,...
                                 'Ch3 | %s',newline,...
                                 '              PowerBin6 | %s', newline,...
                                 '              PowerBin7 | %s'
                                 ...
                                 ],...
                                 ...
                                 vars{:});
          
        
        ld0_vars  = all_vars(contains(all_vars, {'LD0'}) & ... 
                           strcmp(varfun(@class, plt_ss_tbl_oi(:,all_vars),...
                                'OutputFormat','cell'), 'double')...
                           );
            
        
        plt_ld_vars = ld0_vars(plt_ss_tbl_oi{:, ld0_vars} > 0);
        
        
        ld_meta_cell    = table2cell(plt_ss_tbl_oi(1,  plt_ld_vars));
        
        by_ld_meta      = cellfun(@(x,y) [x, ' | ', num2str(y)], ...
                                  plt_ld_vars, ld_meta_cell,...
                                  'UniformOutput',false);
        
        str_form = repmat(['%s', newline], 1, length(by_ld_meta));
        
        
        ld0_meta = sprintf(str_form, by_ld_meta{:});
        
        %%%
        ld1_vars  = all_vars(contains(all_vars, {'LD1'}) & ... 
                           strcmp(varfun(@class, plt_ss_tbl_oi(:,all_vars),...
                                'OutputFormat','cell'), 'double')...
                           );
            
        
        plt_ld_vars = ld1_vars(plt_ss_tbl_oi{:, ld1_vars} > 0);
        
        
        ld_meta_cell    = table2cell(plt_ss_tbl_oi(1,  plt_ld_vars));
        
        by_ld_meta      = cellfun(@(x,y) [x, '    | ', num2str(y)], ...
                                  plt_ld_vars, ld_meta_cell,...
                                  'UniformOutput',false);
        
        str_form = repmat(['%s', newline], 1, length(by_ld_meta));
        
        
        ld1_meta = sprintf(str_form, by_ld_meta{:});
        
        
        %%%
    
        colors = brewermap(5 ,'Set1');
    
        figure('Units', 'Inches', 'Position', [0, 0, 18, 10]);
        
        fig_h = subplot(1,1,1);
        fig_h.Position(4) = 0.6;
        
        scatter(t_vec, movmean(on_off_vec, [30,0]), 5, 'k'); hold on
        
        plot(t_vec, movmean(on_off_vec, [120,0]), 'k');
    
        ylabel(['Percent time ON', newline,'(stim amplitude > 0 mA)'], 'FontSize',18)
        
        ylim([0,40]);
    
        yyaxis right
    
        fig_h.YAxis(2).Color = colors(1,:);
        scatter(redcap.time, redcap.mayoNRS,'filled',...
            'MarkerFaceAlpha', 0.6, ...
            'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:));  
    
        plot(redcap.time, movmean(redcap.mayoNRS, [3, 0], 'omitnan'), ...
            'LineStyle', '-','LineWidth', 2, 'Color', colors(1,:));
    
        ylim([0,10]);
        
        ylabel('NRS intensity', 'FontSize',18)
        fig_h.Position(3) = 0.85;
        fig_h.Position(1) = 0.075;
        %%% plot 
        dur_range = stop_time - start_time;
    
        text(start_time - dur_range/15, fig_h.YLim(2)*1.35, cfg.pt_id, 'FontSize', 32);
        
        text(start_time - dur_range/15, fig_h.YLim(2)*1.2, sense_meta, 'FontSize', 14);
        text(start_time + dur_range/5, fig_h.YLim(2)*1.25, by_ch_pb_meta, 'FontSize',13);
        
        
        text(stop_time - dur_range/1.85, fig_h.YLim(2)*1.2, ld0_meta, 'FontSize',16, 'Interpreter','none');
        text(stop_time - dur_range/4.25, fig_h.YLim(2)*1.2, ld1_meta, 'FontSize',16, 'Interpreter','none');
        
        
        if ge(plt_app_oi.time_INS(end) -plt_app_oi.time_INS(1) , duration('03:00:00'))      
            round_to = 'minute';
        else
            round_to = 'second';
        end


        xticks(...
            dateshift(...
                linspace(start_time, stop_time, 20),...
            "start", round_to));

       
        grid on; grid minor;
    
    
        fig_h.YAxis(1).FontSize = 14;   fig_h.YAxis(2).FontSize = 14;
    
        fig_h.XAxis.FontSize = 14; 
    
    
        set(gca,'xlim', [start_time, stop_time], 'TickLength', [0 0],...
        'GridAlpha',0.4,'MinorGridAlpha',0.7, 'GridColor', 'k', 'MinorGridColor', 'k'); 
    
    
    
        filename = sprintf('%s-->%s.png', ...
                        datestr(dateshift(start_time, 'start', 'day'), 'YY-mm-DD'),...
                        datestr(dateshift(stop_time, 'start', 'day'), 'YY-mm-DD'));
    
        exportgraphics(gcf, [save_dir,filename]);
    end
end

  
% 
% ld0_bias0 = sprintf('%.0f, ' ,round(logspace(log10(11500), log10(512812), 4),-3))
% ld0_bias0 = ld0_bias0(1:end-2)
%mean(on_off_vec)
    %%
    %{
    if length(u_stream_sess) > 1
    
       
    
        starts = [];
        stops = [];
        
        
        for i = 1:length(u_stream_sess)
        
            i_stream_sess_vec = stream_sess_vec == u_stream_sess(i);
        
            if i == 1
                 starts = [starts; 1];
            else
                starts = [starts; find(diff(i_stream_sess_vec) == 1) + 1];  
            end
        
            if i == length(u_stream_sess)
        
                stops = [stops; length(stream_sess_vec)];
            else
                stops = [stops; find(diff(i_stream_sess_vec)== -1)];  
            end
            
            for j = 1:length(starts)
            
                x = [t_vec(starts(j)), t_vec(stops(j)),...
                     t_vec(stops(j)), t_vec(starts(j))];
            
                y = [0,0,50, 50];
            
            
                patch(x, y, colors(j,:),'FaceAlpha', .1,...
                        'HandleVisibility','off');
                hold on
            end
    
        end
    end









    %}