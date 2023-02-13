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


% need to further distill by pts initals
% fluct                 = RCS_redcap_painscores(rcs_API_token, pcs_API_token, {'FLUCT'});

%%
% last 7 days for: 
cfg                     = [];

cfg.pt_id               = 'RCS05';
cfg.dates               = 'PreviousDays';
cfg.ndays               = 10;

cfg.subplot             = true;
cfg.stim_parameter      = '';

    plot_timeline(cfg, REDcap);
    
    
cfg                     = [];

cfg.pt_id               = 'RCS06';
cfg.dates               = 'PreviousDays';
cfg.ndays               = 10;

cfg.subplot             = true;
cfg.stim_parameter      = '';

    plot_timeline(cfg, REDcap);
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

%}

cfg                    = [];
cfg.load_EventLog      = true;
cfg.ignoreold          = false;
cfg.raw_dir            = [pia_dir, 'raw/'];

% specify which patient's INS
pt_sides               = {'RCS02R','RCS05L','RCS05R','RCS07L', 'RCS07R'};


for i =  1:length(pt_sides)

    % process RCS .jsons into searchable database
    cfg.pt_id_side                     = pt_sides{i};
    cfg.proc_dir           = [pia_dir, 'processed/databases/'];

    [db.(pt_sides{i}), bs.(pt_sides{i})]...
        ...
        =  makeDatabaseRCS_Ryan(cfg);

    % now INS logs--inital run can take 15-30 minutes
    cfg.proc_dir           = [pia_dir, 'processed/INS_logs/'];

    INS_logs.(pt_sides{i})...
        ...
        = RCS_logs(cfg);

end

%% unpack all sense, LD, and stimulation settings as own variable in table
% --> allows for programmatic discernment of unique RC+S settings

pt_sides        = {'RCS02R','RCS05L','RCS05R','RCS07L', 'RCS07R'};


cfg                    = [];
cfg.ignoreold          = true;
cfg.raw_dir            = [pia_dir, 'raw/'];
cfg.proc_dir           = [pia_dir, 'processed/parsed_databases/'];


for i= 1:length(pt_sides)

    cfg.pt_id_side = pt_sides{i};


    [par_db.(pt_sides{i}), ss_var_oi.(pt_sides{i})]...
    ...
        = makeParsedDatabaseRCS(...
    ...
    cfg, db);


end

%% find nearest (yet, preceding) streaming session to INS log entry
%%% --> accounts for INS to API time latnecy

for i=  1: length(pt_sides)


   [app_SS_tbl.(pt_sides{i}), INS_logs_proc.(pt_sides{i}), INS_ss_merge_g_changes.(pt_sides{i})] ...
    ...
        = align_stream_sess_to_INSLogs(...
    ...
    INS_logs.(pt_sides{i}), par_db.(pt_sides{i}), ss_var_oi.(pt_sides{i}));
end

%% plot aDBS performance in-real time
cfg             = [];
cfg.dates       = 'AllTime';
cfg.save_dir    = [github_dir, 'Analysis-rcs-data/working/plot_ephy/aDBS_offline_sessions/'];

cfg.dates          = 'DateRange';
cfg.date_range     = {'31-Dec-2022'; '30-May-2023'};

% to avoid figures popping up in MATLAB (see cfg.save_dir for the aDBS longitudinal plots)
set(0,'DefaultFigureVisible','off')

for i=  1 : length(pt_sides)
        cfg.pt_id_side = pt_sides{i};

    plot_longitudinal_aDBS(cfg, REDcap, INS_logs_proc, app_SS_tbl, INS_ss_merge_g_changes);

end

%% per streaming session, visualize INS to API latnecy 
% (i.e., how long is the INS ahead OR behind internet time)
% commented out since only needs to be checked every month or so



%{
%pt_sides               = {'RCS02R','RCS04L','RCS04R','RCS05L','RCS05R','RCS06L','RCS06R','RCS07L','RCS07R'};

pt_sides               = {'RCS02R','RCS05L','RCS05R','RCS07L', 'RCS07R'};

%pt_sides               = {'RCS02R','RCS05L'};
%pt_sides               = {'RCS02R'};

cfg             = [];
cfg.dates       = 'AllTime';
cfg.save_dir    = [github_dir, 'Analysis-rcs-data/working/plot_ephy/aDBS_offline_sessions/'];

for i= 1:length(pt_sides)

    cfg.pt_id_side = pt_sides{i};

    plt_INS_lat_per_session(cfg, db.(pt_sides{i}));
end
%}

%%
%
%
%
%
%
%
%
%%
%% RCS02 --> use INS logs to generate more accurate stim groups
cfg                    = [];
cfg.rootdir            = pia_raw_dir;
cfg.pt_id              = 'RCS02R';
       
INS_logs.RCS02R        = RCS_logs(cfg);

% add stim params to REDcap based off of EventLog.txt, and DeviceSettings.json files
cfg                    = [];
cfg.stage_dates        = stage_dates{2};
cfg.pt_id              = 'RCS02';

INSLog_w_redcap.RCS02R  = get_INSLog_stim_params(cfg,...
                                     INS_logs.RCS02R.group_changes,...
                                     db.RCS02R,...
                                     REDcap.RCS02,...
                                     visits.RCS02...
                                    );

% unilateral implant AND used INS logs to capture PTM intiated group changes
[wrt_stim_REDcap.RCS02, stimGroups.RCS02] ...
    ...
    = make_stim_groups(...
    ...
'RCS02', [], INSLog_w_redcap.RCS02R, visits.RCS02);


cfg                          = [];
cfg.pt_id                    ='RCS02';
cfg.min_n_reports            = 5;

    plted_stim_groups.RCS02  = plot_stim_groups(cfg, stimGroups.RCS02);

%% generate box plots of pain metrics wrt stim parameters
%{
w/ REDcap reports aligned to stimLog.json files
    * group by contacts
    * w/n contacts further group by unique amp-PW-freq-cycling parameters
    - all this is done automatically--even to saving of the pain metrics by
      respective stim groups as .pngs



Next Steps: 

DEFINE baseline pain from pre-Stage 0 fluct data

* evaluate 50% and 30% decrease 
    (see Farrar et al., 2011 for 30% justification based off of patient global
    impression of change (PGIC) of much improved -> very much improved 
    corresponding to 30% decrease across varying pain etiologies, 
    placebo vs pregabalin, ages, btwn, sexes, etc.)


stats framework:

see for multivariate approach (predict many pain metrics)
https://www.mathworks.com/help/stats/specify-the-response-and-design-matrices.html

* RCS02 stage 3:
    * phase 1: ol- blinded testing (visits, and days off)
    * phase 2: cl-parameter testing (for our purposes)
    * phase 3: cl vs ol testing



* run Kolmogorov–Smirnov test (see's if data are normal--explore pain
  metric distribtion more generally)

* run Kruskal-Walis test for the pain metrics across each of the contact-pairs 

    * based on N pain reports, then run Kruskal-Walis test w/n
      amp-freq-PW-cycling space

* (?) follow-up with a Wilconxin signed rank test for signifigant groups

%}

% RCS04, RCS05, RCS06, and RCS07 can be handled together
pts = {'RCS04', 'RCS05', 'RCS06', 'RCS07'};

for i = length(pts)

    [wrt_stim_REDcap.(pts{i}), stimGroups.(pts{i})] ...
    ...
    = make_stim_groups(...
    ...
    pts{i}, REDcap.([pts{i}, 'L']), REDcap.([pts{i}, 'R']), visits.(pts{i}));


    cfg                         = [];
    cfg.min_n_reports           = 5;
    cfg.pt_id                   = pts{i};


    plted_stim_groups.(pts{i})  ...
        ...
        = plot_stim_groups(...
        ...
    cfg, stimGroups.(pts{i}));

end


% RCS04 --> stim groups after starting buprenorphine (after July 2022 home visit)

i_epoch          = find(ge(REDcap.RCS04.time, visits.RCS04.dates(11) + duration('24:00:00')));

[~, stimGroups.RCS04_postJul22] ...
    ...
    = make_stim_groups(...
    ...
 'RCS04', REDcap.RCS04L(i_epoch : end,:), REDcap.RCS04R(i_epoch: end, :), visits.RCS04);


cfg                          = [];
cfg.pt_id                    ='RCS04 (Jul 13th to Now)';
cfg.min_n_reports            = 2;

    plted_stim_groups.RCS04_postJul22 = plot_stim_groups(cfg, stimGroups.RCS04_postJul22);

%% distributions of pain metrics and relationship BETWEEN pain metrics
pts = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};

cfg                     = [];
cfg.dates               = 'AllTime';


for i = 1:length(pts)
    cfg.pt_id  = pts{i};       

    plot_hist(cfg, REDcap);          plot_versus(cfg, REDcap);
end

% RCS02: explore NRS and VAS "mismatch" 
cfg.pt_id               = 'RCS02-mismatch';
cfg.stage_dates         = stage_dates{2}; % starts at Stage 1

      plot_timeline(cfg, REDcap);

%{
Takeaways:
    * VAS 50s dominate, but is otherwise bimodal
    * for most pts, MPQ affective subscore is uninformative
    * NRS has normal(ish) distribution
%}

%% Behavioral clustering: identify low/high pain states + reduce dimensionality
%{
* filter-out 50s (a neccesary alteration to identify natural pain subspaces)

* z-score metrics to themselves for easier comparison of magnitude btwn
    metrics

* cluster based off local density while prioritizing distance
    btwn high density points (Rodriguez & Laio, 2014, Science)

* determine first PC to have summary predictive metric for neural analyses

    * assess/finalize stability of inputs to PC and outputs

%}

pain_space      = [];

cfg             = [];     
cfg.dates       = 'AllTime';
cfg.pca         = false;

for i = 1 : length(pts) - 1

    cfg.pt_id  = pts{i};

    [pain_space.(pts{i})] = plot_pain_space(cfg, REDcap);
end

