cfg             = [];

cfg.date_range  = {'12-Dec-2022', '15-Dec-2022'};

proc_app        = INS_logs.RCS02R.app;



date_range  = datetime(cfg.date_range, 'TimeZone', 'America/Los_Angeles');
i_entries   = find(ge(proc_app.time, date_range(1)) & ...
                 le(proc_app.time, date_range(2)));

entries   = proc_app(i_entries,:);

figure



rcs