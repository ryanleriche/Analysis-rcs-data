pt_id = 'RCS05';

%%

pt_meta          = pt_META.(pt_id);

blind_start       = pt_meta.dates(strcmp(pt_meta.desc, 'blinded_testing_start'));
blind_stop        = pt_meta.dates(strcmp(pt_meta.desc, 'blinded_testing_stop'));

%% from processed stimLog.json and AppLog.txt, focus on blinded testing indices

L_stim_log        = L_stim_log (i_stim_log, :);
%
L_app_log         = INS_logs_API_t_synced.RCS05L.app;
i_app_log         = find(ge(L_app_log.time_INS, blind_start) &...
                         le(L_app_log.time_INS, blind_stop));

L_app_log         = L_app_log(i_app_log, :);

%%
cfg               = [];
cfg.ephy_anal_dir = [dirs.rcs_pia, '/ephy_analysis/aDBS_offline_sessions/'];

cfg.dates         = 'DateRange';
cfg.date_range    = cellstr(datestr([blind_start; blind_stop]));

%%% state-current relationship (12 am - 12 pm)
cfg.plt_state_dur = 'sub_session_duration';

pt_sides           = {'RCS05L','RCS05R'};

for i = 1:length(pt_sides)

    aDBS_sum.(['blinded_', pt_sides{i}]) ...
        ...
        = plot_longitudinal_aDBS(...
        ...
    cfg,    pt_sides{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);
end


%% merge aDBS summary to stimLog



%%

pt_sides           = {'RCS05L', 'RCS05R'};

for i_hemi = 1:length(pt_sides)

    long_aDBS   = aDBS_sum.(['blinded_', pt_sides{i_hemi}]);
    stim_log    = stimLog.(pt_sides{i_hemi});

    i_stim_log  = find(ge(stim_log.time_stimLog, blind_start) &...
                             le(stim_log.time_stimLog, blind_stop));

    stim_log    = stim_log(i_stim_log, :);

    for i_long = 1  :  height(long_aDBS)
        
        if ~isempty(long_aDBS.sess_name{i_long})
            i_stimLog_sess = find(...
                                  ismember(stim_log.sess_name, ...
                                  long_aDBS.sess_name{i_long})...
                                  );
        
        
            stim_log.percentDutyCycle(i_stimLog_sess) = long_aDBS.avg_percent_on(i_long);
        end
    end


end
    par_db_aDBS_ss.(pt_sides{i})


tmp_par_db = par_db.RCS05L;



i_sess = find(strcmp(tmp_par_db.sess_name,'Session1677728486980' ))








