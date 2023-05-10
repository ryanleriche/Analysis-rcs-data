% call CONFIG_aDBS to load directories, and REDcap data 

CONFIG_aDBS;
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

sub_cfg                    = [];
sub_cfg.load_EventLog      = true;

% option to load previous database for efficient processing
sub_cfg.ignoreold_db                = false;
sub_cfg.ignoreold_INS_logs          = false;
sub_cfg.ignoreold_par_db            = true;    % <-- keep as true to avoid version issues

% specify patient hemispheres
%%% pts to update database from scratch locally:
pt_sides        = {'RCS02R','RCS05R', 'RCS05L','RCS04R','RCS04L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};

for i = 1  : length(pt_sides)
    %%% process RCS .jsons into searchable database
    [db.(pt_sides{i}), bs.(pt_sides{i})] ...
        ...
        = makeDatabaseRCS_Ryan(...
        ...
        sub_cfg, pt_sides{i});

    %%% process INS logs .txts based on unique entries only
        % (INS logs have mostly repeating entries)
    INS_logs.(pt_sides{i})  ...
        ...
        = RCS_logs( ...
        ...
        sub_cfg, pt_sides{i});

    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings
    [par_db.(pt_sides{i}), ss_var_oi.(pt_sides{i})] ...
        ...
        = makeParsedDatabaseRCS(...
        ...
        sub_cfg, pt_sides{i}, db);

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
sub_cfg.dates         = 'DateRange';
sub_cfg.date_range    = {'28-Mar-2023'; '01-Jul-2023'};
%cfg.dates        = 'AllTime'; %%% return every aDBS ever tried (takes much longer):

sub_cfg.plt_state_dur = 'sub_session_duration'; %%% state-current relationship (12 am - 12 pm)

%%% plot aDBS performance over requested dates
for i = 1:length(pt_sides)

    aDBS_sum.(pt_sides{i}) ...
        ...
        = plot_longitudinal_aDBS(...
        ...
    sub_cfg,    pt_sides{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);
end