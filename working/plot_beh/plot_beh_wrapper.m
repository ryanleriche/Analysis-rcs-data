%% user-inputs
% where RCS files are saved from PIA server
pia_dir     = '/Users/Leriche/pia_server/datastore_spirit/human/rcs_chronic_pain/rcs_device_data/';

% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
github_dir      = '/Users/Leriche/Github/';

% where DropBox desktop is saved locally
dropbox_dir     = ['/Users/Leriche/Dropbox (UCSF Department of Neurological Surgery)/',...
                   'SUBNETS Dropbox/Chronic Pain - Activa and Summit 2.0'];

% application programming interface (API) token which is essentially a
% password to access REDcap remotely, and is unique per researcher per
% study (e.g., Ryan has a unique token for the RCS and PCS studies)

rcs_API_token   = '95FDE91411C10BF91FD77328169F7E1B';
pcs_API_token   = 'DB65F8CB50CFED9CA5A250EFD30F10DB';

% pulls/organizes arms from REDcap (go into fxn to add new arms)
cd([github_dir, 'Analysis-rcs-data/working']);         

addpath(genpath([github_dir, 'Analysis-rcs-data/']));
addpath(genpath([github_dir, 'rcs-simulation/']));

REDcap                  = RCS_redcap_painscores(rcs_API_token);

% stage dates and, home/clinic visits for RCS pts 1-7 w/ brief descriptions
[visits, stage_dates]   = make_visit_dates;

% per RCS pt, organize pain fluctuation study (pain prior to stage 0)
%{
If this fails b/c of a time out issue go into MATLAB's 'webwrite()'fxn
and below line ~126   --> 'options = validateRequestMethod(options);' 
add                   --> 'options.Timeout   = 20;'
%}  

fluct             = RCS_redcap_painscores(rcs_API_token, pcs_API_token, {'FLUCT'});

% pain fluctuation study descriptive stats provides basis for clinical
% outcomes of interest
cfg               = [];
cfg.dates         = 'AllTime';
pts               = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};
for i = 1:length(pts)
    fluct_sum_stats.(pts{i})  = calc_sum_stats(cfg, fluct.(pts{i}));
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