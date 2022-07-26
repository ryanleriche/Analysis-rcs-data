%% user-inputs

% where RCS files are saved from PIA server
pia_raw_dir               = '/Users/Leriche/pia_server/datastore_spirit/human/rcs_chronic_pain/rcs_device_data/raw';

% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
github_dir            = '/Users/Leriche/Github/Analysis-rcs-data/';

% application programming interface (API) token which is essentially a
% password to access REDcap remotely, and is unique per researcher per
% study (e.g., Ryan has a unique token for the RCS and PCS studies)

API_token             = '95FDE91411C10BF91FD77328169F7E1B';


% pulls/organizes arms from REDcap (go into fxn to add new arms)
cd([github_dir, 'working']);         addpath(genpath(github_dir));

pt_pain                 = RCS_redcap_painscores(API_token);

% Stage start dates I-III for RCS pts 1-6
stage_dates             = [{''}, {'08-Sep-2020', '31-Jan-2021', '31-May-2022'},... % RCS02
                           {''}, {'13-May-2021', '12-Jul-2021'},... % RCS04
                           {'21-Jul-2021','07-Sep-2021'},... % RCS05
                           {'07-Jun-2022'}]; % RCS06

%% import RCS DataBses per pt side

% save RCS session summaries as DataBases w/n 'db' structure and
% badsessions w/n 'bs' structure 

% RCS02

[db.RCS02R, bs.RCS02R] = ...
    makeDataBaseRCSdata(pia_raw_dir, 'RCS02R');


% RCS04

[db.RCS04L, bs.RCS04L] = ...
    makeDbRCS(pia_raw_dir, 'RCS04L','ignoreold');



[db.RCS04L, bs.RCS04L] = ...
     makeDataBaseRCSdata(pia_raw_dir, 'RCS04L');



[db.RCS04R, bs.RCS04R] = ...
    makeDataBaseRCSdata(pia_raw_dir, 'RCS04R');


% RCS05
[db.RCS05L, bs.RCS05L] = ...
    makeDataBaseRCSdata(pia_raw_dir, 'RCS05L','ignoreold');

[db.RCS05R, bs.RCS05R] = ...
    makeDataBaseRCSdata(pia_raw_dir, 'RCS05R');


%% exploring how to best align "nearest" REDcap report to pain metrics

nearest_beh_to_sess = interp1(pt_pain.RCS04.time, pt_pain.RCS04.time,...
        db_RCSXXX.time, 'nearest');

i_wn_30_min = le(abs(nearest_beh_to_sess - db_RCSXXX.time),...
                 '00:30:00');


RCSXX_paintable = db_RCSXXX(i_wn_30_min,:);





%% RCS04 plot daily metrics

% Specify, cfg before calling functions--see below for examples.
cfg                     = [];
cfg.pt_id               = 'RCS04';

cfg.dates               = 'AllTime';
cfg.stage_dates         = {'13-May-2021', '19-Jul-2021'}; % starts at Stage 1
cfg.subplot             = true;

cfg.stim_parameter      = '';
   
    plot_timeline(cfg, pt_pain.RCS04, db.RCS04L);


cfg.dates       = 'DateRange';
cfg.date_range  = {'01-Jun-2022'; '30-Jun-2022'};

    plot_timeline(cfg, pt_pain.RCS04, db.RCS04L) ;


cfg.dates       = 'DateRange';
cfg.date_range  = {'13-Jul-2022'; '19-Jul-2022'};
cfg.subplot     = true;
    plot_timeline(cfg, pt_pain.RCS04) ;

    % add color brewer and remove edges of patches

% visually inspect pain metric distributions
cfg             = [];
cfg.pt_id       = 'RCS04';
cfg.dates       = 'AllTime';

    plot_hist(cfg, pt_pain.RCS04);

cfg            = [];
cfg.dates      = 'AllTime';
% SEE Table for easy access to common summary statistics
RCS04_sum_stats      = calc_sum_stats(cfg, pt_pain.RCS04);

% pain metric A "versus" pain metric B --> see how metrics visually covary
cfg             = [];
cfg.pt_id       = 'RCS04';
cfg.dates       = 'AllTime';

% go into this fxn to change the pain metrics on each axis

% RBL reccomends leaving NRS as the colorbar though, as it is a categorical
% metric--visually creating planes which muddles visualization
    plot_versus(cfg, pt_pain.RCS04);


%% 07/21/22 Lab Meeting 

% 'plot_timeline' demo
cfg                     = [];
cfg.pt_id               = 'RCS04';

cfg.dates               = 'AllTime';
cfg.stage_dates         = {'13-May-2021', '12-Jul-2021'}; % starts at Stage 1

cfg.dates               = 'PreviousDays';
cfg.date_range          = {'13-Jul-2022'; '19-Jul-2022'};
cfg.ndays               = 7;
cfg.subplot             = true;

cfg.stim_parameter      = '';
   
    plot_timeline(cfg, pt_pain.RCS04, db.RCS04L, db.RCS04R);




% 'plot_hist' demo

cfg.pt_id       = 'RCS02';       plot_hist(cfg, pt_pain.RCS02);

cfg.pt_id       = 'RCS04';       plot_hist(cfg, pt_pain.RCS04);

cfg.pt_id       = 'RCS05';       plot_hist(cfg, pt_pain.RCS05);

% 'plot_versus' demo

cfg.pt_id       = 'RCS02';       plot_versus(cfg, pt_pain.RCS02);

cfg.pt_id       = 'RCS04';       plot_versus(cfg, pt_pain.RCS04);

cfg.pt_id       = 'RCS05';       plot_versus(cfg, pt_pain.RCS05);

%%

% clean-up fields for easy visualization
%{
    if height(stimLogSettings) > 1
        db_out(d).stimparams = stimLogSettings;

    else

                            
        stim_para = strsplit(char(db_out(d).stimparams),',');
        
        db_out(d).contacts = stim_para{1};
        db_out(d).amp      = str2double(stim_para{2}(2:end -2));
        db_out(d).PW       = str2double(stim_para{3}(2:end -2));
        db_out(d).freq     = str2double(stim_para{4}(2:end -2));


    end

%}




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

n_runs          = 5;            t_off           = 20;
n_iter          = 2;            t_on            = 60;


t_in_sec        = n_iter * n_runs * (t_off+ t_on);
duty_cycle      = t_on / (t_off + t_on);

t_in_min        = t_in_sec / 60

% logrhythmically spaced frequencies
freq_start       = 10;
freq_stop        = 150;

freq_exp_start   = log10(freq_start);
freq_exp_stop    = log10(freq_stop);

log_freqs        = floor(logspace(freq_exp_start, freq_exp_stop, n_runs))

lin_freqs        = floor(linspace(freq_start, freq_stop, n_runs))

amp              = repmat(2, 1, n_runs)
pw               = repmat(300, 1, n_runs)

% element-wise multiplication based on percent of time stim-sweep is on
TEED             = amp.^2 .* lin_freqs .* pw * (duty_cycle)

% need impendance to find mean total electrical energy delievered (TEED) per sec

TEED_now         = 2.^2 * 100 * 300 * 60/(60+20)

TEED_stim_sweep_percentage   = TEED ./ TEED_now





% Stim Comparisons
%{

07/15/22 - 07/29/22

    RACC
        open-loop
            C:     0+3-     2 mA     100 Hz     300 mcs     60s/20s

    *w/ option to decrease amplitude

    Q: If lowering the amplitude 0.1-0.5 mA does not help after 5-6 hours
    -->
        Lower amplitude further (min of 0.5 mA)?
        Keep waiting?
        Turn off stim?

    Q: Min N of days per contact-pair?

        Q: Within given contact-pair min N of days at given
           amplitude, pulse width, frequency?


Stim Sweeps as therapy (originally IJ's idea)


    
   analgesic effect of:

        *cumulative regular stim sweeps (daily, weekly, etc.)
    
            LCaudB:     C+10-    1 mA     125 Hz     300 mcs     60s/20s 


        "meditation-like" calming by focusing on internal state

            Q: Try double-blinded sham stim sweeps?

        

07/22/22, Stim Sweep
    
    (Program right C for reference)
    RACC
            C:     0+3-     2 mA     100 Hz     300 mcs     60s/20s
    -->
        freq:
        [80, 0, 90, 0, 100, 0, 110, 0, 100, 0, 90, 100, 0, 110]

        PW:
        [0, 200, 0, 250, 0, 300, 0, 200, 0, 300, 0, 250, 0]


        * 2 w/ PW & freq values shuffled while maintaining zeros btwn.








___________________________________________________________________________
bilateral v. unilateral stim

cumalative effect of stim

    regular stim-sweeps

    total energy delieverd as fxn of pain

%}

% Control Comparisons
%{

%}