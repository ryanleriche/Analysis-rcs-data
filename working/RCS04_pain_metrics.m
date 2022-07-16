%% user-inputs
% the patient's ID and specific side as outputted by RCS
PATIENTIDside         = 'RCS04L';

% where RCS files are saved from PIA server
rootdir               = '/Users/Leriche/pia_server/datastore_spirit/human/rcs_chronic_pain/rcs_device_data/raw';

% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
github_dir            = '/Users/Leriche/Github/';

% application programming interface (API) token which is essentially a
% password to access REDcap remotely, and is unique per researcher per
% study (e.g., Ryan has a unique token for the RCS and PCS studies)

API_token             = '95FDE91411C10BF91FD77328169F7E1B';


% pulls/organizes arms from REDcap (go into fxn to add new arms)
cd(github_dir);         addpath(genpath(github_dir));

pt_pain                 = RCS_redcap_painscores(API_token);

%% import RCS files
patientroot_dir = fullfile(rootdir,char(regexp(PATIENTIDside,...
                        '\w*\d\d','match'))); %match the PATIENTID up to 2 digits: ie RCS02

% scbs_dir        = [patientroot_dir,...
%                     '/SummitData/SummitContinuousBilateralStreaming/', PATIENTIDside];

[RCSdatabase_out, badsessions] = ...
    makeDataBaseRCSdata(patientroot_dir, PATIENTIDside, 'ignoreold');


%% RCS04 plot daily metrics

% Specify, cfg before calling functions--see below for examples.
cfg             = [];
cfg.pt_id       = 'RCS04';
cfg.dates       = 'AllTime';
cfg.stage_dates = {'13-May-2021', '12-Jul-2021'}; % starts at Stage 1
   
    plot_timeline(cfg, pt_pain.RCS04);


cfg.dates       = 'DateRange';
cfg.date_range  = {'01-Jun-2022'; '30-Jun-2022'};

    plot_timeline(cfg, pt_pain.RCS04) ;


cfg.dates       = 'DateRange';
cfg.date_range  = {'13-Jul-2022'; '16-Jul-2022'};
    plot_timeline(cfg, pt_pain.RCS04) ;


% visually inspect pain metric distributions
cfg             = [];
cfg.pt_id       = 'RCS04';
cfg.dates       = 'AllTime';

    plot_hist(cfg, pt_pain.RCS04);

cfg            = [];
cfg.dates      = 'AllTime';
% SEE Table for easy access to common summary statistics
RCS04_sum_stats      = calc_sum_stats(cfg, pt_pain.RCS04);

%%
% pain metric A "versus" pain metric B --> see how metrics visually covary
cfg             = [];
cfg.pt_id       = 'RCS04';
cfg.dates       = 'AllTime';

% go into this fxn to change the pain metrics on each axis

% RBL reccomends leaving NRS as the colorbar though, as it is a categorical
% metric--visually creating planes which muddles visualization
    plot_versus(cfg, pt_pain.RCS04);

%%

cfg             = [];
cfg.pt_id       = 'RCS05';
cfg.dates       = 'AllTime';

   plot_versus(cfg, pt_pain.RCS05);


%% RCS04 Stim Plan

% 07/13/22; Home Programs
%{
LCaud
    open-loop
        A:     9+11-    1 mA     125 Hz     300 mcs     60s/20s
        B:     C+10-    1 mA     125 Hz     300 mcs     60s/20s  
        C:     C+9-     1 mA     125 Hz     300 mcs     60s/20s

    closed-loop
        D:     C+10-    0-1 mA   100 Hz     300 mcs     aDBS

RThal
    open-loop
        A:     9+11-    1 mA     150 Hz     200 mcs     60s/20s
        B:     C+11-    3 mA     100 Hz     200 mcs     60s/20s

    closed-loop
        D:     C+10-    0-1 mA   100 Hz     300 mcs     aDBS

RACC
    open-loop
        C:     0+3-     2 mA     100 Hz     300 mcs     60s/20s

%}

% Stim Comparisons
%{

07/15/22 - 07/29/22

    RACC
        open-loop
            C:     0+3-     2 mA     100 Hz     300 mcs     60s/20s

    *w/ option to decrease amplitude

    If lowering the amplitude 0.1-0.5 mA does not help after 5-6 hours
    -->
        Lower amplitude further (min of 0.5 mA)?
        Keep waiting?
        Turn off stim? 

07/20/22, Stim Sweep
    
    IF pain was controlled at given amplitude (given expected compliance of
    2 mA max w/ option to decrease amplitude).
    -->
        sweep
        [80, 0, 90, 0, 100, 0, 110, 0, 100, 0, 90, 100, 0, 110]

    Else
    -->
        sweep
        [0, 2, 0, 3, 0, 2.5, 0, 2.75, 0, 3, 0, 2.5, 0, 2, 0, 2.75]


    What other sweep should we try?
        IJ mentioned









___________________________________________________________________________
bilateral v. unilateral stim

cumalative effect of stim

    regular stim-sweeps

    total energy delieverd as fxn of pain

%}

% Control Comparisons
%{

%}





















