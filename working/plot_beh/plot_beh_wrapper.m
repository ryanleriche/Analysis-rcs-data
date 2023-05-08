% load CONFIG_beh_wrapper

[dirs,    rcs_API_token,   pcs_API_token, ... -> input and output directories, and API tokens
 PFS,     PFS_sum_stats,...                   -> pain fluctuation study (PFS) data and summary statistics
 pt_META, stage_dates]...                     -> hard-coded pt meta data (RCS Stage dates pulled from patient iternaries, Google Drive folders, etc
...
    = CONFIG_beh_wrapper;

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

    -> Note LAST stim and sense setting is returned for parsimonious report
       of aDBS settings between streaming sessions

*   from db create a stimLog containing every change in stim parameter
    during a streaming session ("misses" offline PTM intiated changed)


(pg. 12 of the 4NR010 Research Lab Programmer Guide M979053A001)

"Multiple programs can be combined to define a group. Each group and its associated programs can be used
to provide a therapy for specific symptoms or specific patient activities. Pulse width, amplitude, amplitude
limits, and electrode polarity are programmed separately for each program within the group (ie, each
program within the group can have different values). Pulse width limits, rate, rate limits, SoftStart/Stop,
Cycling, and Active Recharge are programmed for each group (ie, each program within the group will have
the same values)."

Inital run can take hours if running multiple pts w/ 1000s of streaming
sessions.
%}

cfg                    = [];
cfg.load_EventLog      = true;

% option to load previous database for efficient processing
cfg.ignoreold_db                = false;
cfg.ignoreold_INS_logs          = true;
cfg.ignoreold_par_db            = true;

cfg.raw_dir                     = [dirs.rcs_pia, 'raw/'];
cfg.proc_dir                    = [dirs.rcs_pia, 'processed/'];

cfg.ephy_anal_dir               = [dirs.rcs_pia, '/ephy_analysis/aDBS_offline_sessions/'];

% specify patient hemispheres
%%% pts to update database from scratch locally:
pt_sides        = {'RCS02R','RCS05R', 'RCS05L','RCS04R','RCS04L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};

%%% pts for NEJM submission:
%pt_sides           = {'RCS02R', 'RCS05R', 'RCS05L'};

%pt_sides           = {'RCS04R','RCS04L'};

%pt_sides           = {'RCS02R','RCS05L','RCS05R'};

for i = 1  : length(pt_sides)
    %%% process RCS .jsons into searchable database
    [db.(pt_sides{i}), bs.(pt_sides{i})] ...
        ...
        = makeDatabaseRCS_Ryan(...
        ...
        cfg, pt_sides{i});

    %%% process INS logs .txts based on unique entries only
        % (INS logs have mostly repeating entries)
    INS_logs.(pt_sides{i})  ...
        ...
        = RCS_logs( ...
        ...
        cfg, pt_sides{i});


    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings
    [par_db.(pt_sides{i}), ss_var_oi.(pt_sides{i})] ...
        ...
        = makeParsedDatabaseRCS(...
        ...
        cfg, pt_sides{i}, db);

    %%% finally merge stim parameters to REDcap survey during said parameters
    [stimLog.(pt_sides{i}), REDcap.(pt_sides{i})]  ...
        ...
        = align_REDcap_to_stimLog(...
        ...
        cfg, pt_sides{i}, db, REDcap);
end
%%
%{
EventLog.txt          -> tracks Group changes and TherapyStatus offline

stimLog.json          -> tracks all Group changes during streaming sessions
                         (misses cycling and other settings)

DeviceSettings.json   -> shows comphrensive stimulation settings 

Start from stimLog.json and add comprehensive stimulation settings 
from nearest previous DeviceSettings.json (time_stimLog)

W/ API-latency corrected INS_log entries (built parsimoniously from
EventLog.txt files), find nearest previous stimLog entry and infer
comprensive settings

merge the comprehensive INS_log and stimLog entires into single table w/:

* time_aligned column
    - datetime vector in PST from the INS_log or stimLog
    from whenever *either* changed

* source
    -cell array listing Left or Right and if INS or API
        -(L_INS_time, L_API_time, R_INS_time, R_API_time)
%}

pt_sides           = {'RCS02R', 'RCS05R','RCS05L'};


for i = 1: length(pt_sides)
    %%% find nearest (yet, preceding) streaming session to INS log entry
    % --> accounts for INS to API time latency
    [par_db_aDBS_ss.(pt_sides{i}), INS_logs_API_t_synced.(pt_sides{i})] ...
    ...
        = align_INSLogs_to_API_time(...
    ...
    pt_sides{i}, INS_logs, par_db, ss_var_oi);
end


%%
%%% specify which dates to return:
cfg.dates         = 'DateRange';
cfg.date_range    = {'28-Mar-2023'; '01-Jul-2023'};

%%% return every aDBS ever tried (takes much longer):
%cfg.dates        = 'AllTime';

%%% state-current relationship (12 am - 12 pm)
cfg.plt_state_dur = 'sub_session_duration';

%%% state-current relationship (from 1-2 am and 1-2 pm):
%cfg.plt_state_dur = 'two_chunks'; 

%%% plot aDBS performance over months
% w/ aligned INS logs, plot requested dates
%pt_sides           = {'RCS02R','RCS05R','RCS05L'};

pt_sides           = {'RCS02R','RCS05R','RCS05L'};

for i = 1:length(pt_sides)

    aDBS_sum.(pt_sides{i}) ...
        ...
        = plot_longitudinal_aDBS(...
        ...
    cfg,    pt_sides{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);
end



%%
% last N days for: 
cfg                     = [];

cfg.pt_id               = 'RCS05';
cfg.dates               = 'PreviousDays';
cfg.ndays               = 10;

%%% return every pain score:
%cfg.dates        = 'AllTime';

cfg.subplot             = true;
cfg.stim_parameter      = '';

    plot_timeline(cfg, REDcap, fluct_sum_stats);
    
    
cfg                     = [];

cfg.pt_id               = 'RCS06';
cfg.dates               = 'PreviousDays';
cfg.ndays               = 10;

%%% return every pain score:
%cfg.dates        = 'AllTime';

cfg.subplot             = true;
cfg.stim_parameter      = '';

    plot_timeline(cfg, REDcap, fluct_sum_stats);



cfg                     = [];

cfg.pt_id               = 'RCS07';
cfg.dates               = 'PreviousDays';
cfg.ndays               = 10;

%%% return every pain score:
%cfg.dates        = 'AllTime';

cfg.subplot             = true;
cfg.stim_parameter      = '';

    plot_timeline(cfg, REDcap, fluct_sum_stats);

%% import RCS databases, and INS logs per pt side
%{

* saves RCS session summaries as databases (db) in struct per pt side
* badsessions w/n 'bs' struct

(pg. 12 of the 4NR010 Research Lab Programmer Guide M979053A001)

"Multiple programs can be combined to define a group. Each group and its associated programs can be used
to provide a therapy for specific symptoms or specific patient activities. Pulse width, amplitude, amplitude
limits, and electrode polarity are programmed separately for each program within the group (ie, each
program within the group can have different values). Pulse width limits, rate, rate limits, SoftStart/Stop,
Cycling, and Active Recharge are programmed for each group (ie, each program within the group will have
the same values)."

Inital run can take hours if running multiple pts w/ 1000s of streaming
sessions.
%}

cfg                    = [];
cfg.load_EventLog      = true;

% option to load previous database for efficient processing
cfg.ignoreold          = false;
cfg.raw_dir            = [pia_dir, 'raw/'];

% specify patient hemispheres
cfg.pt_sides         = {'RCS02R', 'RCS04R','RCS04L', 'RCS05L','RCS05R',...
                        'RCS06R','RCS06L', 'RCS07L', 'RCS07R'};

% pts not yet on closed-loop:
% cfg.pt_sides        = {'RCS04R','RCS04L', 'RCS06R','RCS06L'};

cfg.proc_dir           = [pia_dir, 'processed/databases/'];

    % process RCS .jsons into searchable database
    [db, bs] =  makeDatabaseRCS_Ryan(cfg);

%%% process INS logs .txts based on unique entries only
% (INS logs have mostly repeating entries)
% same cfg as 'makeDatabaseRCS()', but w/ different output folder

cfg.proc_dir  = [pia_dir, 'processed/INS_logs/'];

    INS_logs  = RCS_logs(cfg);

%% per streaming session, visualize INS to API latency and impedance
% (i.e., how long is the INS ahead OR behind internet time [generally behind])

% commented out as no need to frequently run:

%{
% prevents figures from popping up (see in cfg.proc_dir folder)
set(0,'DefaultFigureVisible','off')

%cfg             = [];
cfg.dates       = 'AllTime';
cfg.proc_dir    = [pia_dir, 'processed/aDBS_offline_sessions/'];

    plt_INS_lat_per_session(cfg, db);

%%% plot impedance over time (contacts to case during a "Lead Integrity Test")
cfg.proc_dir    = [pia_dir, 'processed/impedance_per_session/'];

    plt_impedance_per_session(cfg, db);

set(0,'DefaultFigureVisible','on')
%}
%% unpack all sense, LD, and stimulation settings as own variable in table
% --> allows for programmatic discernment of unique RC+S settings

cfg                  = [];
cfg.ignoreold        = false;
cfg.raw_dir          = [pia_dir, 'raw/'];
cfg.proc_dir         = [pia_dir, 'processed/parsed_databases/'];

cfg.pt_sides         = {'RCS02R', 'RCS04R','RCS04L', 'RCS05L','RCS05R',...
                        'RCS06R','RCS06L', 'RCS07L', 'RCS07R'};

    [par_db, ss_var_oi]    = makeParsedDatabaseRCS(cfg, db);

%% plot aDBS performance over months

cfg             = [];
cfg.pt_sides    = {'RCS02R', 'RCS05L', 'RCS05R', 'RCS07L', 'RCS07R'};
cfg.proc_dir    = [pia_dir, 'processed/aDBS_offline_sessions/'];

%%% specify which dates to return:
cfg.dates         = 'DateRange';
cfg.date_range    = {'01-Jan-2023'; '30-May-2023'};

%%% return every aDBS ever tried (takes much longer):
%cfg.dates        = 'AllTime';


%%% state-current relationship (12 am - 12 pm)
cfg.plt_state_dur = 'sub_session_duration';

%%% state-current relationship (from 1-2 am and 1-2 pm):
%cfg.plt_state_dur = 'two_chunks'; 


%%% find nearest (yet, preceding) streaming session to INS log entry
% --> accounts for INS to API time latency
   [app_SS_tbl, INS_logs_proc, INS_ss_merge_g_changes] ...
    ...
        = align_stream_sess_to_INSLogs(...
    ...
    cfg, INS_logs, par_db, ss_var_oi);
             
% w/ aligned INS logs, plot requested dates
   plot_longitudinal_aDBS(cfg, REDcap, INS_logs_proc, app_SS_tbl, INS_ss_merge_g_changes);