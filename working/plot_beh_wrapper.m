CONFIG_plot_beh;
%% import REDcap daily, weekly, and monthly surveys from stages 1,2 and 3
% as of Apr. 2023, only daily surveys are analysis-ready/organized

REDcap                  = RCS_redcap_painscores(rcs_API_token);

%% import RCS databases, and INS logs per pt side as structures
%{
*   saves RCS session summaries as databases (db) from constellation of
    .jsons saved during streaming sessions

*   saves INS logs (with AppLog.txt, and EventLog.txt changes) as INS_logs
    importantly contains Adaptive state changes/stim defintions and Group changes

*   from db created a parsed database (par_db) allowing for .xlsx
    exportable and human readable one line summaries per streaming session

    -> Note LAST stim/sense settings are returned for parsimonious report
       of aDBS settings between streaming sessions

*   from db create a stimLog containing every change in stim parameter
    during a streaming session ("misses" offline PTM intiated changed)

Inital run takes hours for running multiple pts w/ 1000s of streaming
sessions.
%}
cfg_rcs_db.load_EventLog      = true;

% option to load previous database for efficient processing
cfg_rcs_db.ignoreold_db                = false;
cfg_rcs_db.ignoreold_INS_logs          = false;
cfg_rcs_db.ignoreold_par_db            = true;    % <-- keep as true to avoid version issues

% specify patient hemispheres
%%% pts to update database from scratch locally:
pt_sides        = {'RCS02R','RCS05R', 'RCS05L','RCS04R','RCS04L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};


for i = 1  : length(pt_sides)
    %%% process RCS .jsons into searchable database
    [db.(pt_sides{i}), bs.(pt_sides{i})] ...
        ...
        = makeDatabaseRCS_Ryan(...
        ...
        cfg_rcs_db, pt_sides{i});

    %%% process INS logs .txts based on unique entries only
        % (INS logs have mostly repeating entries)
    INS_logs.(pt_sides{i})  ...
        ...
        = RCS_logs( ...
        ...
        cfg_rcs_db, pt_sides{i});


    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings
    [par_db.(pt_sides{i}), ss_var_oi.(pt_sides{i})] ...
        ...
        = makeParsedDatabaseRCS(...
        ...
        cfg_rcs_db, pt_sides{i}, db);

    %%% find nearest (yet, preceding) streaming session to INS log entry
        % accounts for INS to API time latency
    [par_db_aDBS_ss.(pt_sides{i}), INS_logs_API_t_synced.(pt_sides{i})] ...
    ...
        = align_INSLogs_to_API_time(...
    ...
    pt_sides{i}, INS_logs, par_db, ss_var_oi);

end
%% plot aDBS performance longitudinally
%%% specify which dates to return:
cfg_rcs_db.dates         = 'DateRange';
cfg_rcs_db.date_range    = {'28-Mar-2023'; '01-Jul-2023'};

%%% return every aDBS ever tried (takes much longer):
%cfg.dates        = 'AllTime';

%%% state-current relationship (12 am - 12 pm)
cfg_rcs_db.plt_state_dur = 'sub_session_duration';

%%% state-current relationship (from 1-2 am and 1-2 pm):
%cfg.plt_state_dur = 'two_chunks'; 

%%% plot aDBS performance over months
% w/ aligned INS logs, plot requested dates
for i = 1:length(pt_sides)

    aDBS_sum.(pt_sides{i}) ...
        ...
        = plot_longitudinal_aDBS(...
        ...
    cfg_rcs_db,    pt_sides{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);
end


%%%
%% behavioral plotting only
%%
% last N days for: 
cfg_rcap                     = [];

cfg_rcap.pt_id               = 'RCS05';
cfg_rcap.dates               = 'PreviousDays';
cfg_rcap.ndays               = 10;

%%% return every pain score:
%cfg.dates        = 'AllTime';

cfg_rcap.subplot             = true;
cfg_rcap.stim_parameter      = '';

    plot_timeline(cfg_rcap, REDcap, fluct_sum_stats);
    
    
cfg_rcap                     = [];

cfg_rcap.pt_id               = 'RCS06';
cfg_rcap.dates               = 'PreviousDays';
cfg_rcap.ndays               = 10;

%%% return every pain score:
%cfg.dates        = 'AllTime';

cfg_rcap.subplot             = true;
cfg_rcap.stim_parameter      = '';

    plot_timeline(cfg_rcap, REDcap, fluct_sum_stats);



cfg_rcap                     = [];

cfg_rcap.pt_id               = 'RCS07';
cfg_rcap.dates               = 'PreviousDays';
cfg_rcap.ndays               = 10;

%%% return every pain score:
%cfg.dates        = 'AllTime';

cfg_rcap.subplot             = true;
cfg_rcap.stim_parameter      = '';

    plot_timeline(cfg_rcap, REDcap, fluct_sum_stats);