
 function [app_oi, ss_tbl_oi, g_chan_oi] = aDBS_toi(cfg, proc_app, app_ss_tbl, proc_g_chan)
       switch cfg.dates
           case 'DateRange'

                date_range  = datetime(cfg.date_range, 'TimeZone', 'America/Los_Angeles', 'InputFormat','dd-MMM-uuuu');
                
                i_entries   = find(ge(proc_app.time_INS, date_range(1)) & ...
                                 le(proc_app.time_INS, date_range(2)));
                app_oi      = proc_app(i_entries,:);
            
            
                i_entries   = find(ge(app_ss_tbl.timeStart, date_range(1)) & ...
                                 le(app_ss_tbl.timeStart, date_range(2)));
                ss_tbl_oi   = app_ss_tbl(i_entries,:);
            
            
                i_entries   = find(ge(proc_g_chan.time_INS, date_range(1)) & ...
                                 le(proc_g_chan.time_INS, date_range(2)));
            
                g_chan_oi   = proc_g_chan(i_entries, :);
    
    
           case 'AllTime'
    
                app_oi      = proc_app;
                ss_tbl_oi   = app_ss_tbl;
                g_chan_oi   = proc_g_chan;
    
       end
end