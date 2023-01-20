cfg             = [];
%cfg.dates       = 'AllTime';

cfg.dates       = 'DateRange';
cfg.date_range  = {'12-Dec-2022', '22-Jan-2023'};

proc_app        = INS_logs.RCS02R.app;
proc_group      = group_changes;
%proc_group      = INS_logs.RCS02R.group_changes;
%%
if strcmp(cfg.dates, 'DateRange') == 1

    date_range  = datetime(cfg.date_range, 'TimeZone', 'America/Los_Angeles');
    
    i_entries   = find(ge(proc_app.time, date_range(1)) & ...
                     le(proc_app.time, date_range(2)));
    
    app_oi   = proc_app(i_entries,:);


    i_entries   = find(ge(proc_group.time_INS, date_range(1)) & ...
                     le(proc_group.time_INS, date_range(2)));
    
    group_oi   = proc_group(i_entries,:);



end





