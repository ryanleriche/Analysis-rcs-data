%% user-inputs

% where RCS files are saved from PIA server
pia_raw_dir               = '/Users/Leriche/pia_server/datastore_spirit/human/rcs_chronic_pain/rcs_device_data/raw/';

% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
github_dir            = '/Users/Leriche/Github/Analysis-rcs-data/';

% application programming interface (API) token which is essentially a
% password to access REDcap remotely, and is unique per researcher per
% study (e.g., Ryan has a unique token for the RCS and PCS studies)

API_token             = '95FDE91411C10BF91FD77328169F7E1B';


% pulls/organizes arms from REDcap (go into fxn to add new arms)
cd([github_dir, 'working']);         addpath(genpath(github_dir));

REDcap                 = RCS_redcap_painscores(API_token);

%% Stage start dates, home, and clinic visits for RCS pts 1-6
stage_dates             = {{''}, {'08-Sep-2020'; '31-Jan-2021'; '31-May-2022'},... % RCS02
                           {''}, {'13-May-2021'; '12-Jul-2021'},... % RCS04
                           {'21-Jul-2021';'07-Sep-2021'},... % RCS05
                           {'07-Jun-2022','18-Aug-2022'}}; % RCS06


visits              = struct;
visits.RCS05        = table;
visits.RCS05.dates  = datetime(...
   {'21-Jul-2021';   '22-Jul-2021';    '23-Jul-2021';...
    '07-Sep-2021';   '08-Sep-2021';    '09-Sep-2021';...
    '01-Dec-2021';   '02-Dec-2021';...
    '10-Aug-2022';   '11-Aug-2022';...
    '05-Sep-2022'},...
    ...
    'TimeZone', 'America/Los_Angeles');

visits.RCS05.desc   = ...
    {'s1_implant';      's1_inpatient_d1';   's1_inpatient_d2';...
    's2_inclinic_d1';   's2_inclinic_d2' ;   's2_inclinic_d3';...
    's2_home_d1'    ;   's2_home_d2';...
    's2_inclinic_d1';   's2_inclinic_d2';...
    's3_start'};




visits.RCS02        = table;
visits.RCS02.dates  = datetime(...
   {'08-Sep-2020'; '09-Sep-2020'; '10-Sep-2020'; '11-Sep-2020';...
    '13-Oct-2020';...
    '18-Oct-2020'; ...
    '02-Nov-2021'; '09-Nov-2021';...
    '01-Feb-2021'; '02-Feb-2021';  '03-Feb-2021';...
    '13-Apr-2021'; '14-Apr-2021';...
    '27-Sep-2021'; '28-Sep-2021';...
    '31-May-2022';...
    '28-Jun-2022'; '29-Jun-2022'}, ...
    ...
    'TimeZone', 'America/Los_Angeles');

visits.RCS02.desc   = ...
    {'s1_implant';     's1_inpatient';    's1_inclinic_d1';  's1_inclinic_d2';...
     'remove_hardware_inpatient';...
     's2_start';
     'washout_testing'; 'washout_testing';
     's2_inclinic_d1'; 's2_inclinic_d2' ; 's2_inclinic_d3';...
     's2_home_d1'    ; 's2_home_d2';...
     's2_inclinic_d1'; 's2_inclinic_d2' ;...
     's3_start';...
     's3_inclinic_d1'; 's3_inclinic_d2'};


visits.RCS04        = table;
visits.RCS04.dates  = datetime(...
   {'13-May-2021'; '14-May-2021'; '15-May-2021'; '16-May-2021';...   
    '12-Jul-2021'; '13-Jul-2021'; '14-Jul-2021';...
    '17-Aug-2021'; '18-Aug-2021'; '19-Aug-2021';...
    '13-Jul-2022'}, ...
    ...
    'TimeZone', 'America/Los_Angeles');

visits.RCS04.desc   = ...
    {'s1_implant';      's1_inpatient_d1';  's1_inpatient_d2'; 's1_inpatient_d3'; ...
     's2_inclinic_d1';  's2_inclinic_d2' ;  's2_inclinic_d3';...
     's2_home_d1'    ;  's2_home_d2';       's2_home_d2';...
     's2_home_d1'    ...
     };



visits.RCS06        = table;
visits.RCS06.dates  = datetime(...
   {'01-Mar-2022';...   
    '07-Jun-2021'}, ...
    ...
    'TimeZone', 'America/Los_Angeles');

visits.RCS06.desc   = ...
    {'s1_implant';...
     's2'};


%% import RCS Databases per pt side

% save RCS session summaries as DataBases w/n 'db' structure and
% badsessions w/n 'bs' structure 

cfg                    = [];
cfg.load_EventLog      = false;
cfg.ignoreold          = false;




% RCS04
[db.RCS04L, bs.RCS04L] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS04L', cfg);

[db.RCS04R, bs.RCS04R] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS04R', cfg);

% RCS05
[db.RCS05L, bs.RCS05L] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS05L',cfg);

[db.RCS05R, bs.RCS05R] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS05R',cfg);

% RCS02
[db.RCS02R, bs.RCS02R] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS02R',cfg);


%% RCS04 plot daily metrics
% Note leave 'cfg.stim_parameter' blank as alignment lacks INS log data
% right now

% Specify, cfg before calling functions--see below for examples.
cfg                     = [];
cfg.pt_id               = 'RCS04';

cfg.dates               = 'AllTime';
cfg.stage_dates         = stage_dates{4}; % starts at Stage 1
cfg.subplot             = false;

cfg.stim_parameter      = '';
   
    plot_timeline(cfg, REDcap.RCS04);

%% last 7 days for all pts
cfg                     = [];
cfg.pt_id               = 'RCS04';
cfg.stage_dates         = stage_dates{4}; % starts at Stage 1
cfg.subplot             = true;

cfg.stim_parameter      = '';

cfg.dates               = 'PreviousDays';
cfg.ndays               = 7;
cfg.subplot             = true;

    plot_timeline(cfg, REDcap.RCS04);

cfg.pt_id               = 'RCS05';
cfg.stage_dates         = stage_dates{5}; % starts at Stage 1

        plot_timeline(cfg, REDcap.RCS05);


cfg.pt_id               = 'RCS02';
cfg.stage_dates         = stage_dates{2}; % starts at Stage 1

      plot_timeline(cfg, REDcap.RCS02);


cfg.pt_id               = 'RCS06';
cfg.stage_dates         = stage_dates{6}; % starts at Stage 1

      plot_timeline(cfg, REDcap.RCS06);
      

% visually inspect pain metric distributions
cfg             = [];
cfg.pt_id       = 'RCS04';
cfg.dates       = 'AllTime';

    plot_hist(cfg, REDcap.RCS04);

cfg            = [];
cfg.dates      = 'AllTime';
% SEE Table for easy access to common summary statistics
RCS04_sum_stats      = calc_sum_stats(cfg, REDcap.RCS04);


%% ************************************************************************
% clustering behavioral metrics to identify low/high pain states + reduce dimensionality


% build intuition of beh distribution via histograms
cfg             = [];
cfg.dates       = 'AllTime';
cfg.pt_id       = 'RCS02';       

plot_hist(cfg, REDcap.RCS02);

% saveas(gcf, [cd, '/plot_beh/figs/RCS02/', 'RCS02_beh_hist.png']);

cfg.pt_id       = 'RCS04';       plot_hist(cfg, REDcap.RCS04);

% saveas(gcf, [cd, '/plot_beh/figs/RCS04/', 'RCS04_beh_hist.png']);


cfg.pt_id       = 'RCS05';       plot_hist(cfg, REDcap.RCS05);

% saveas(gcf, [cd, '/plot_beh/figs/RCS05/', 'RCS05_beh_hist.png']);


%{
Conclusion:
 VAS 50s dominate, but is otherwise bimodal, MPQ affective is uninformative, 
 and NRS has normal(ish) distribution

Next Steps:

* filter-out 50s (a neccesary alteration to identify natural pain subspaces)

* z-score metrics to themselves for easier comparison of magnitude btwn
    metrics

* cluster based off local density while prioritizing distance
    btwn high density points (Rodriguez & Laio, 2014, Science)

* determine first PC to have summary predictive metric for neural analyses

    * assess/finalize stability of inputs to PC and outputs

%}
