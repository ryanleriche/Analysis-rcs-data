%% user-inputs

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


%% align RCSdatabase stim parameters to REDcap pain metrics

% import RCS files

% the patient's ID and specific side as outputted by RCS
PATIENTIDside         = 'RCS04L';

patientroot_dir = fullfile(rootdir,char(regexp(PATIENTIDside,...
                        '\w*\d\d','match'))); %match the PATIENTID up to 2 digits: ie RCS02

[db.RCS04L, db.RCS04L_badsessions] = ...
    makeDataBaseRCSdata(patientroot_dir, PATIENTIDside);

PATIENTIDside         = 'RCS04R';
[db.RCS04R, db.RCS04R_badsessions] = ...
    makeDataBaseRCSdata(patientroot_dir, PATIENTIDside);



%% parse through and split stim parameters into their own columns

db.RCS04R = parse_db(db.RCS04R);

db.RCS04L = parse_db(db.RCS04L);



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
cfg.stage_dates         = {'13-May-2021', '12-Jul-2021'}; % starts at Stage 1
cfg.subplot             = false;

cfg.stim_parameter      = 'contacts';
   
    plot_timeline(cfg, pt_pain.RCS04, db_RCSXXX);


cfg.dates       = 'DateRange';
cfg.date_range  = {'01-Jun-2022'; '30-Jun-2022'};

    plot_timeline(cfg, pt_pain.RCS04) ;


cfg.dates       = 'DateRange';
cfg.date_range  = {'13-Jul-2022'; '18-Jul-2022'};
cfg.subplot     = true;
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
cfg.pt_id               = 'RCS04L';

cfg.dates               = 'AllTime';
cfg.stage_dates         = {'13-May-2021', '12-Jul-2021'}; % starts at Stage 1
cfg.subplot             = false;

cfg.stim_parameter      = 'contacts';
   
    plot_timeline(cfg, pt_pain.RCS04, db.RCS04L);

cfg.pt_id               = 'RCS04R';

    plot_timeline(cfg, pt_pain.RCS04, db.RCS04R);



% 'plot_hist' demo

   plot_hist(cfg, pt_pain.RCS02);


cfg.pt_id       = 'RCS04';       plot_hist(cfg, pt_pain.RCS04);
cfg.pt_id       = 'RCS05';       plot_hist(cfg, pt_pain.RCS05);

% 'plot_versus' demo

cfg.pt_id       = 'RCS02';       plot_versus(cfg, pt_pain.RCS02);

cfg.pt_id       = 'RCS04';       plot_versus(cfg, pt_pain.RCS04);

cfg.pt_id       = 'RCS05';       plot_versus(cfg, pt_pain.RCS05);




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

%% local functions
function db_RCSXXX = parse_db(db_RCSXXX)

    for i_sess = 1 : height(db_RCSXXX)
    
        if length(db_RCSXXX.stimparams{i_sess}) <= 1
    
            i_para = strsplit(char(db_RCSXXX.stimparams{i_sess}),',');
        
            if ~isempty(i_para{1})
        
                db_RCSXXX.contacts{i_sess}    = i_para{1};
                db_RCSXXX.amp(i_sess)         = str2double(i_para{2}(2:end -2));
                db_RCSXXX.PW(i_sess)          = str2double(i_para{3}(2:end -2));
                db_RCSXXX.freq(i_sess)        = str2double(i_para{4}(2:end -2));
        
            else % without stim parameters 
                db_RCSXXX.contacts{i_sess}    = '';
                db_RCSXXX.amp(i_sess)         = NaN;
                db_RCSXXX.PW(i_sess)          = NaN;
                db_RCSXXX.freq(i_sess)        = NaN;
    
            end
    
        else
            db_RCSXXX.contacts{i_sess}    = 'MANY';
            db_RCSXXX.amp(i_sess)         = NaN;
            db_RCSXXX.PW(i_sess)          = NaN;
            db_RCSXXX.freq(i_sess)        = NaN;
        end
    end
end