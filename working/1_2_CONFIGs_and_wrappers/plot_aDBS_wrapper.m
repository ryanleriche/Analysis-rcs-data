%% see user_local_CONFIG.m script to specify directories and general set-up
user_local_CONFIG;

% load REDcap pain surveys and set-up directories
general_setup;

% see "dirs" structure for main directories used, and "rcs_db_cfg" for the
% configuration of the RC+S databasing

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

% option to load previous database for efficient processing
    rcs_db_cfg.ignoreold_db                = false;
    rcs_db_cfg.ignoreold_INS_logs          = false;
    rcs_db_cfg.ignoreold_par_db            = true;    % <-- keep as true to avoid version issues
    
    %rcs_db_cfg.ignoreold_aDBS_plts         = false;

% specify which dates to plot aDBS longitudinal plots:
    %rcs_db_cfg.dates        = 'AllTime'; %%% return every aDBS ever tried (takes much longer):
    rcs_db_cfg.dates         = 'DateRange';
    rcs_db_cfg.date_range    = {'01-Jul-2023'; '3-Jun-2027'};

% specify patient hemispheres
    pt_sides        = {'RCS02R', 'RCS04L','RCS04R', 'RCS05L', 'RCS05R',...
                       'RCS06L', 'RCS06R','RCS07L', 'RCS07R'};

for i = 1  : length(pt_sides)
    %%% process RCS .jsons into searchable database
    [db.(pt_sides{i}), bs.(pt_sides{i})] ...
        ...
        = makeDatabaseRCS_Ryan(...
        ...
        rcs_db_cfg, pt_sides{i});

    %%% process INS logs .txts based on unique entries only
        % (INS logs have mostly repeating entries)
    INS_logs.(pt_sides{i})  ...
        ...
        = RCS_logs( ...
        ...
        rcs_db_cfg, pt_sides{i});

    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings
    [par_db.(pt_sides{i}), ss_var_oi.(pt_sides{i})] ...
        ...
        = makeParsedDatabaseRCS(...
        ...
        rcs_db_cfg, pt_sides{i}, db);

    %%% find nearest (yet, preceding) streaming session to INS log entry
        % accounts for INS to API time latency
    [par_db_aDBS_ss.(pt_sides{i}), INS_logs_API_t_synced.(pt_sides{i})] ...
    ...
        = align_INSLogs_to_API_time(...
        ...
        pt_sides{i}, INS_logs, par_db, ss_var_oi);

    %%% plot aDBS performance over requested dates
        aDBS_sum.(pt_sides{i}) ...
        ...
        = plot_longitudinal_aDBS_2(...
        ...
        rcs_db_cfg,    pt_sides{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);

end