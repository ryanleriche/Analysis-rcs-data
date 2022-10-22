%% user-inputs

% where RCS files are saved from PIA server
pia_raw_dir               = '/Users/Leriche/pia_server/datastore_spirit/human/rcs_chronic_pain/rcs_device_data/raw/';

% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
github_dir            = '/Users/Leriche/Github/Analysis-rcs-data/';

% application programming interface (API) token which is essentially a
% password to access REDcap remotely, and is unique per researcher per
% study (e.g., Ryan has a unique token for the RCS and PCS studies)

rcs_API_token             = '95FDE91411C10BF91FD77328169F7E1B';
pcs_API_token             = 'DB65F8CB50CFED9CA5A250EFD30F10DB';

% pulls/organizes arms from REDcap (go into fxn to add new arms)
cd([github_dir, 'working']);         addpath(genpath(github_dir));

REDcap                 = RCS_redcap_painscores(rcs_API_token);

% fluct = RCS_redcap_painscores(rcs_API_token, pcs_API_token, {'FLUCT'});

%% Stage start dates, home, and clinic visits for RCS pts 1-6
stage_dates             = {{''}, {'08-Sep-2020'; '31-Jan-2021'; '31-May-2022'},... % RCS02
                           {''}, {'13-May-2021'; '12-Jul-2021'},... % RCS04
                           {'21-Jul-2021';'07-Sep-2021'},... % RCS05
                           {'07-Jun-2022','18-Aug-2022'}}; % RCS06



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




visits.RCS06        = table;
visits.RCS06.dates  = datetime(...
   {'07-Jun-2022'; '08-Jun-2022'; '09-Jun-2022';...
    '17-Aug-2022'; '18-Aug-2022'}, ...
    ...
    'TimeZone', 'America/Los_Angeles');

visits.RCS06.desc   = ...
    {'s1_implant';  's1_inpatient_d1';  's1_inpatient_d2';...
     's2_inclinic_d1';  's2_inclinic_d2'};

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


% RCS06
[db.RCS06L, bs.RCS06L] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS06L',cfg);

[db.RCS06R, bs.RCS06R] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS06R',cfg);


%% import RCS logs (state changes from INS)
cfg = [];
cfg.pull_adpt_logs      = false;
cfg.pull_event_logs     = true;
cfg.ld_detection_events = false;
cfg.pull_recharge_sess  = false;

cfg.pull_mirror_logs    = false;

%{
X see how INS log changes compare to db_RCSXX.timeStarts
    - 'timeStarts' are in 'PacketGenTime' (API estimate of packet was generated on
       the INS as in Sellers et al., 2021 in Front. Hum. Neurosci.)

    - INS logs are in 'timestamp' which is implemented in the INS firmware
      which has 1s resolution (Sellers et al., 2021)


determine when w/n streaming session REDcap survey occured

    using RCS data start and stop times

    go into 'read_adaptive_txt_log' and further deeper fxns for all needed
    fields

    
    


%}

[textlog.RCS04L]         = RCS_logs(pia_raw_dir,'RCS04L', cfg);

[textlog.RCS04R]         = RCS_logs(pia_raw_dir,'RCS04R', cfg);


[textlog.RCS05R]         = RCS_logs(pia_raw_dir,'RCS05R', cfg);

%% concatenate StimLog.json outputs from RCS database

cfg                             = [];
cfg.stage_dates                 = stage_dates{2};
cfg.pt_id                       = 'RCS02';


% RCS02
[stimLog_w_redcap.RCS02R] = ...
    align_REDcap_to_stimLog(cfg, db.RCS02R, REDcap.RCS02);


% RCS04
cfg.stage_dates = stage_dates{4};
cfg.pt_id       = 'RCS04';

[stimLog_w_redcap.RCS04L] = ...
    align_REDcap_to_stimLog(cfg, db.RCS04L, REDcap.RCS04);

[stimLog_w_redcap.RCS04R] = ...
    align_REDcap_to_stimLog(cfg, db.RCS04R, REDcap.RCS04);


% RCS05
cfg.stage_dates = stage_dates{5};
cfg.pt_id       = 'RCS05';

[stimLog_w_redcap.RCS05L] = ...
    align_REDcap_to_stimLog(cfg, db.RCS05L, REDcap.RCS05);

[stimLog_w_redcap.RCS05R] = ...
    align_REDcap_to_stimLog(cfg, db.RCS05R, REDcap.RCS05);

% RCS06
cfg.stage_dates = stage_dates{6};
cfg.pt_id       = 'RCS06';

[stimLog_w_redcap.RCS06L] = ...
    align_REDcap_to_stimLog(cfg, db.RCS06L, REDcap.RCS06);

[stimLog_w_redcap.RCS06R] = ...
    align_REDcap_to_stimLog(cfg, db.RCS06R, REDcap.RCS06);
%% stim parameters grouping based off of 'DeviceSettings.json' from 'StimSettingsOut' var in db.RCSXXX
%{
%{
* creates db_beh.RCSXX with timestamps, contacts, amp, PW, freq, cycling,
  group, and laterality (unilateral vs bilateral stim)

    -manually inspect laterality to handle edge cases

* align nearest REDcap report as row in db_beh.RCSXX


for RCS05 and RCS02 especially look at pain reports within session 5 - 15 minutes
   

%}


% cfg                 = [];
% cfg.pt_id           = 'RCS04';
% cfg.pia_raw_dir     = pia_raw_dir;
% 
% cfg.plot_sess_dur    = false;
% 

cfg                 = [];
cfg.pt_id           = 'RCS04';
cfg.pia_raw_dir     = pia_raw_dir;
cfg.plot_sess_dur   = false;

   db_beh.RCS04 = db_sort_beh(cfg, db.RCS04L, db.RCS04R, REDcap.RCS04);


cfg.pt_id = 'RCS05';
    db_beh.RCS05 = db_sort_beh(cfg, db.RCS05L, db.RCS05R, REDcap.RCS05, visits.RCS05);


cfg.pt_id = 'RCS02';
    db_beh.RCS02 = db_sort_beh(cfg,[], db.RCS02R, REDcap.RCS02);


cfg.pt_id = 'RCS06';
    [db_beh.RCS06] = db_sort_beh(cfg,db.RCS06L, db.RCS06R, REDcap.RCS06, visits.RCS06  );

%}


%% validate StimLog.json's w/ db.RCSXXX 'stimSettingsOut' (from DeviceSettings.json)

%{
for the REDcap reports reported in/near a steaming session compare the stim
parameters btwn. that and the StimLog.json



%}
% code (commented as used for sanity check btwn StimLog.json and DeviceSettings.json)
%{
i_stimLog_redcaps = cell2mat(stimLog.i_redcap);

i_emp = cellfun(@isempty, db_beh.RCS05.stimContacts);

i_L = cellfun(@(x) strcmp('L', x(1)), db_beh.RCS05.stimContacts(~i_emp));

db_beh_RCS05L = db_beh.RCS05(~i_emp,:);

db_beh_RCS05L = db_beh_RCS05L(i_L,:);


for i = 1 : height(db_beh_RCS05L)

   % per db_beh, return index in stimLog if ANY of the REDcap indices are shared
   in_stimLog_at = find(cellfun(@(x) any(ismember(db_beh_RCS05L.i_redcap_near{i}, x)), stimLog.i_redcap));

    if  any(~strcmp(db_beh_RCS05L.stimGroup(i), stimLog.activeGroup(in_stimLog_at)))||...
        any(db_beh_RCS05L.stimfreq(i)  ~=  stimLog.stimfreq(in_stimLog_at)) ||...
        any(db_beh_RCS05L.stimPW(i)  ~=  stimLog.stimPW(in_stimLog_at))
    
        disp([num2str(i), ' index of db_beh_RCS05L does not match ', num2str(in_stimLog_at'),...
            ' index of StimLog.json (for "iredcap_near")'])

    end

     % per db_beh, return index in stimLog if ANY of the REDcap indices are shared
   in_stimLog_at = find(cellfun(@(x) any(ismember(db_beh_RCS05L.i_redcap{i}, x)), stimLog.i_redcap));

    if  any(~strcmp(db_beh_RCS05L.stimGroup(i), stimLog.activeGroup(in_stimLog_at)))||...
        any(db_beh_RCS05L.stimfreq(i)  ~=  stimLog.stimfreq(in_stimLog_at)) ||...
        any(db_beh_RCS05L.stimPW(i)  ~=  stimLog.stimPW(in_stimLog_at))
    
        disp([num2str(i), ' index of db_beh_RCS05L does not match ', num2str(in_stimLog_at'),...
            ' index of StimLog.json (for "iredcap")'])

    end

end


all_db_beh_i_redcap = vertcat(db_beh_RCS05L.i_redcap{:});

if length(all_db_beh_i_redcap) ~= length(unique(all_db_beh_i_redcap)) 
    
    disp('****The same REDcap report(s) "assigned" to different stim settings*****')

end

per_assigned = length(unique(all_db_beh_i_redcap)) ./ height(redcap) * 100;

if  per_assigned ~= 1

    disp(['Only ', num2str(per_assigned),'% of REDcap report(s) assigned to stim settings'])

end
%}
%% generate box plots of pain metrics wrt stim parameters
%{

input(s)
cfg

align_REDcap_to_stimLog.RCSXXL
align_REDcap_to_stimLog.RCSXXR

REDcap.fluct
___________________________________________________________________________

* w/ aligned REDcap reports, first visualize the parameter space per each
   "major contact pair" (best analgesic contact pairs as assessed clinically--
    heuristically those contacts w/ a 10+ pain reports)

* visualized explored amp-PW-freq space per contact
    * quantify N reports nested by freq, amp, and PW

* box plots of NRS, VAS, MPQ, and PC1 (for RCS05 only) based of
  parsimonious amp-PW-freq space *per* contact


___________________________________________________________________________

output(s)
beh_ss


DEFINE baseline pain from pre-Stage 0 fluct data

* evaluate 50% and 30% decrease 
    (see Farrar et al., 2011 for 30% justification based off of patient global
    impression of change (PGIC) of much improved -> very much improved 
    corresponding to 30% decrease across varying pain etiologies, 
    placebo vs pregabalin, ages, btwn, sexes, etc.)


RCS04 stats framework:

see for multivariate approach (predict many pain metrics)
https://www.mathworks.com/help/stats/specify-the-response-and-design-matrices.html


* run Kolmogorov–Smirnov test (see's if data are normal--explore pain
  metric distribtion more generally)

* run Kruskal-Walis test for the pain metrics across each of the contact-pairs 
  (grouping all freq, amp, and PW parameters w/n a contact--lack statistical 
  power for freq/PW/amp)

* for LCaud, RThal, and bilateral stim contacts look run run a Kruskal-Walis 
  test across amplitude w/n contacts 

* follow-up with a Wilconxin signed rank test (?) for signifigant groups

%}

% RCS02
[wrt_stim_REDcap.RCS02, stimGroups.RCS02] ...
    ...
    = make_stim_groups(...
    ...
'RCS02', [], stimLog_w_redcap.RCS02R, REDcap.RCS02, visits.RCS02);


cfg                          = [];
cfg.pt_id                    ='RCS02';
cfg.min_n_reports            = 5;

    plted_stim_groups.RCS02  = plot_stim_groups(cfg, stimGroups.RCS02);



% RCS04
[wrt_stim_REDcap.RCS04, stimGroups.RCS04] ...
    ...
    = make_stim_groups(...
    ...
 'RCS04', stimLog_w_redcap.RCS04L, stimLog_w_redcap.RCS04R, REDcap.RCS04, visits.RCS04);


cfg                          = [];
cfg.pt_id                    ='RCS04';
cfg.min_n_reports            = 5;

    plted_stim_groups.RCS04  = plot_stim_groups(cfg, stimGroups.RCS04);



% RCS05
[wrt_stim_REDcap.RCS05, stimGroups.RCS05] ...
    ...
    = make_stim_groups(...
    ...
'RCS05', stimLog_w_redcap.RCS05L, stimLog_w_redcap.RCS05R, REDcap.RCS05, visits.RCS05);

cfg                          = [];
cfg.pt_id                    ='RCS05';
cfg.min_n_reports            = 5;

    plted_stim_groups.RCS05  = plot_stim_groups(cfg, stimGroups.RCS05);


% RCS06
[wrt_stim_REDcap.RCS06, stimGroups.RCS06] ...
    ...
    = make_stim_groups(...
    ...
'RCS06', stimLog_w_redcap.RCS06L, stimLog_w_redcap.RCS06R, REDcap.RCS06, visits.RCS06);

cfg                          = [];
cfg.pt_id                    ='RCS06';
cfg.min_n_reports            = 5;

    plted_stim_groups.RCS06  = plot_stim_groups(cfg, stimGroups.RCS06);



%% RCS04 plot daily metrics

% Specify, cfg before calling functions--see below for examples.
cfg                     = [];
cfg.pt_id               = 'RCS04';

cfg.dates               = 'AllTime';
cfg.stage_dates         = stage_dates{4}; % starts at Stage 1
cfg.subplot             = false;

cfg.stim_parameters     = {'contacts'};
   
    plot_timeline(cfg, REDcap.RCS04);

%% last 7 days for all pts
cfg                     = [];
cfg.pt_id               = 'RCS04';
cfg.stage_dates         = stage_dates{4}; % starts at Stage 1
cfg.subplot             = true;

cfg.stim_parameter      = '';

cfg.dates               = 'PreviousDays';
cfg.ndays               = 10;
cfg.subplot             = true;

    plot_timeline(cfg, REDcap.RCS04);

cfg.pt_id               = 'RCS05';
cfg.stage_dates         = stage_dates{5}; % starts at Stage 1

        plot_timeline(cfg, REDcap.RCS05, db_beh.RCS05);


cfg.pt_id               = 'RCS02';
cfg.stage_dates         = stage_dates{2}; % starts at Stage 1

      plot_timeline(cfg, REDcap.RCS02, db_beh.RCS05);


cfg.ndays               = 21;
cfg.pt_id               = 'RCS06';
cfg.stage_dates         = stage_dates{6}; % starts at Stage 1

      plot_timeline(cfg, REDcap.RCS06, db_beh.RCS06);

%% inspect Stage I
cfg                     = [];
cfg.pt_id               = 'RCS04';
cfg.stage_dates         = stage_dates{4}; % starts at Stage 1
cfg.subplot             = true;

cfg.stim_parameter      = '';

cfg.dates               = 'DateRange';
cfg.date_range          = stage_dates{4}(1:2);
cfg.subplot             = true;

    plot_timeline(cfg, REDcap.RCS04, db_beh.RCS05);

% RCS05
cfg                     = [];
cfg.pt_id               = 'RCS05';
cfg.stage_dates         = stage_dates{5}; % starts at Stage 1
cfg.subplot             = true;

cfg.stim_parameter      = '';

cfg.dates               = 'DateRange';
cfg.date_range          = stage_dates{5}(1:2);
cfg.subplot             = true;

    plot_timeline(cfg, REDcap.RCS05, db_beh.RCS05);

% RCS02
cfg                     = [];
cfg.pt_id               = 'RCS02';
cfg.stage_dates         = stage_dates{2}; % starts at Stage 1
cfg.subplot             = true;

cfg.stim_parameter      = '';

cfg.dates               = 'DateRange';
cfg.date_range          = stage_dates{2}(1:2);
cfg.subplot             = true;

    plot_timeline(cfg, REDcap.RCS02, db_beh.RCS05);

% RCS06
cfg                     = [];
cfg.pt_id               = 'RCS06';
cfg.stage_dates         = stage_dates{6}; % starts at Stage 1
cfg.subplot             = true;

cfg.stim_parameter      = '';

cfg.dates               = 'DateRange';
cfg.date_range          = stage_dates{6}(1:2);
cfg.subplot             = true;

    plot_timeline(cfg, REDcap.RCS06, db_beh.RCS05);
%%
cfg                     = [];
cfg.pt_id               = 'RCS04';
cfg.stage_dates         = stage_dates{4}; % starts at Stage 1
cfg.subplot             = true;

cfg.stim_parameter      = 'all';

cfg.dates               = 'DateRange';
cfg.date_range               = {'14-Jul-2022'; '6-Sep-2022'};
cfg.subplot             = false;

    plot_timeline(cfg,...
        wrt_stim_REDcap.RCS04);


%% visually inspect pain metric distributions

cfg             = [];
cfg.pt_id       = 'RCS04';
cfg.dates       = 'AllTime';

    plot_hist(cfg, REDcap.RCS04);


cfg            = [];
cfg.dates      = 'AllTime';
% SEE Table for easy access to common summary statistics
RCS04_sum_stats      = calc_sum_stats(cfg, REDcap.RCS04);

% pain metric A "versus" pain metric B --> see how metrics visually covary
cfg             = [];
cfg.pt_id       = 'RCS04';
cfg.dates       = 'AllTime';

% go into this fxn to change the pain metrics on each axis

% RBL reccomends leaving NRS as the colorbar though, as it is a categorical
% metric--visually creating planes which muddles visualization
    plot_versus(cfg, REDcap.RCS04);



%% ************************************************************************
% clustering behavioral metrics to identify low/high pain states + reduce dimensionality


% build intuition of beh distribution via histograms
cfg             = [];
cfg.dates               = 'DateRange';
cfg.date_range          = stage_dates{4}(1:2);
cfg.pt_id       = 'RCS02';       

plot_hist(cfg, REDcap.RCS02, db_beh.RCS05);

saveas(gcf, [cd, '/plot_beh/figs/RCS02/', 'RCS02_beh_hist.png']);

cfg.pt_id       = 'RCS04';       plot_hist(cfg, REDcap.RCS04);

saveas(gcf, [cd, '/plot_beh/figs/RCS04/', 'RCS04_beh_hist.png']);

cfg.date_range          = stage_dates{5}(1:2);
cfg.pt_id       = 'RCS05';       plot_hist(cfg, REDcap.RCS05, db_beh.RCS05);

saveas(gcf, [cd, '/plot_beh/figs/RCS05/', 'RCS05_beh_hist.png']);


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

%%
cfg             = [];
cfg.pt_id       = 'RCS02';       
cfg.dates       = 'AllTime';
cfg.pca         = false;

[beh_fig.RCS02_filt, beh_fig.RCS02_z_space, beh_fig.RCS02_dec, beh_fig.RCS02_cl, ...
    beh_anl.z_RCS02] ...
    =...
    plot_versus(cfg, REDcap.RCS02);


% rinse & repeat for RCS04
cfg.pt_id       = 'RCS04'; 

[beh_fig.RCS04_filt, beh_fig.RCS04_z_space, beh_fig.RCS04_dec, beh_fig.RCS04_cl,...
    beh_anl.z_RCS04] ...
    =...
    plot_versus(cfg, REDcap.RCS04);



% rinse & repeat for RCS05
cfg.pt_id       = 'RCS05'; 

[beh_fig.RCS05_filt, beh_fig.RCS05_z_space, beh_fig.RCS05_dec, beh_fig.RCS05_cl, ...
    beh_anl.z_RCS05] ...
    =...
    plot_versus(cfg, REDcap.RCS05);

% rinse & repeat for RCS06
cfg.pt_id       = 'RCS06'; 

[beh_fig.RCS06_filt, beh_fig.RCS06_z_space, beh_fig.RCS06_dec, beh_fig.RCS06_cl, ...
    beh_anl.z_RCS06] ...
    =...
    plot_versus(cfg, REDcap.RCS06);

% save([cd,'/plot_beh/beh_anl'],'beh_anl');

save([cd,'/plot_beh/beh_anl_no_pca'],'beh_anl');


fig_names = fieldnames(beh_fig);
for i = 1 :  length(fig_names)

    savefig(beh_fig.(fig_names{i}),[cd,'/plot_beh/figs/', fig_names{i},'.fig'] )

    saveas(beh_fig.(fig_names{i}),...
        [cd,'/plot_beh/figs/', fig_names{i}, '.png']);

end


beh_anl_out     = fieldnames(beh_anl);

n_pts           = length(beh_anl_out);

for i = 1 :  n_pts

    writetable(beh_anl.(beh_anl_out{i}),...
        ['/Users/Leriche/Dropbox (UCSF Department of Neurological Surgery)/UFlorida_UCSF_RCS_collab/Pain Reports/beh_clustered/',...
             beh_anl_out{i}(3:end)], 'FileType','spreadsheet')

end


%%

%{

"ask question and it's opposite"

Chronback's coefficient alpha

%}


%% organize Streaming Notes, clinic dates, etc


%{
RCS06
    08/17/22 clinic visit

    instructed on Streaming Notes to use HIS timezone 
    
    sliders on pt tablet "jumps"--needs to zoom in to get VAS right where he
    wants it

    CTM occasionally drops during sessions--probably the right side--not
    100%

    for MPQ "sickening" field has been answered--last 3-7 days Crohn's Flare
%}












openfig([cd,'/plot_beh/figs/RCS02_z_space.fig'])

openfig([cd,'/plot_beh/figs/RCS04_z_space.fig'])

openfig([cd,'/plot_beh/figs/RCS05_z_space.fig'])

% visualize distribution of PC 1--wait for now--could be useful for
% distributed closed-loop

%{

figure('Units', 'Inches', 'Position', [0, 0, 15, 10])

sgtitle('Dist of PC 1 per patient')


subplot(311)
title('RCS02');
hold on;
histogram(beh_anl.RCS02_pc.score(:,1)) 

set(gca,'fontSize',14, 'TickLength', [0 0]); 
        grid on;    grid MINOR;      box off


subplot(312)
title('RCS04')
hold on
histogram(beh_anl.RCS04_pc.score(:,1)) 

set(gca,'fontSize',14, 'TickLength', [0 0]); 
        grid on;    grid MINOR;      box off

subplot(313)
title('RCS05');
hold on;
histogram(beh_anl.RCS05_pc.score(:,1)) 

xlabel('PC 1')

 set(gca,'fontSize',14, 'TickLength', [0 0]); 
        grid on;    grid MINOR;      box off
%}

%% animate pain space visualization
%{

v = RecordRotation(h);

%%



playRecording(v);
%%
function V2 = RecordRotation(h)

n = 250;
V = zeros(n,2);
T = timer('period',0.05,'executionmode','fixedrate',...
    'TimerFcn',@captureAzEl,'TasksToExecute',n);
rotate3d on
drawnow;
start(T);
drawnow;
wait(T);
V2 = V;

      function captureAzEl(src,~)
          cnt = get(src,'TasksExecuted');
          [V(cnt,1), V(cnt,2)] = view;  
          drawnow;
      end

end

%%

function playRecording(V)
frame_period = 0.01;
n = length(V);
T = timer('period',frame_period,'executionmode','fixedrate',...
    'TimerFcn',@playRecording,'TasksToExecute',n);
rotate3d on
drawnow;
start(T);
drawnow;
wait(T);

      function playRecording(src,~)
          cnt = get(src,'TasksExecuted');
          view(V(cnt,1), V(cnt,2));  
          drawnow;        
      end
  end
%}


%% RCS04 Stim Plan

% 07/13/22; Home Programs
%{
LCaud
    open-loop
        A:     9+11-    1 mA     125 Hz     300 mcs     60s/20s
        B:     C+10-    2 mA     125 Hz     300 mcs     60s/20s  
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
        C:     0+3-     2 mA     100 Hz     300 mcs     60s/20sk

%}

% experiment w/ stim sweeps parameters
%{
n_runs          = 5;            t_off           = 20;
n_iter          = 2;            t_on            = 60;


t_in_sec        = n_iter * n_runs * (t_off+ t_on);
duty_cycle      = t_on / (t_off + t_on);

t_in_min        = t_in_sec / 60

% logrhythmically spaced frequencies
freq_start       = 10;
freq_stop        = 175;

freq_exp_start   = log10(freq_start);
freq_exp_stop    = log10(freq_stop);

log_freqs        = floor(logspace(freq_exp_start, freq_exp_stop, n_runs))

lin_freqs        = floor(linspace(freq_start, freq_stop, n_runs))

% amp              = repmat(2, 1, n_runs)
pw               = repmat(300, 1, n_runs)

% element-wise multiplication based on percent of time stim-sweep is on
% TEED             = amp.^2 .* lin_freqs .* pw * (duty_cycle)

% need impendance to find mean total electrical energy delievered (TEED) per sec

TEED_now         = 2.^2 * 100 * 300 * 60/(60+20)

amp              = sqrt(TEED_now ./ (lin_freqs.* pw .* duty_cycle))


%}



n_runs          = 5;            t_off           = 20;
n_iter          = 2;            t_on            = 60;


t_in_sec        = n_iter * n_runs * (t_off+ t_on);
duty_cycle      = t_on / (t_off + t_on);

t_in_min        = t_in_sec / 60




%%
cont_I  = repmat(["c+1-", "c+2-", "1+2-", "c+1-", "c+2-", "1+2-"; ...
                         "1", "1", "1", "2", "2", "2"], 1, 4);


shuffled_cont_I  = cont_I(:,randperm(length(cont_I)))';


cont_II  = repmat(["0+3-", "0+2-", "0+1-", "0+3-", "0+2-", "0+1-"; ...
                         "1", "1", "1", "2", "2", "2"], 1, 4);


shuffled_cont_II  = cont_II(:,randperm(length(cont_II)))';


cont_III  = repmat(["X", "Y", "X", "Y"; "0", "0", "2", "2"], 1, 4);

shuffled_cont_III  = cont_III(:,randperm(length(cont_III)))';




length(shuffled_cont_I) + length(shuffled_cont_II) * 40* 2 / 60

%%
amps             = repmat([1 2], 1, 3);
shuffled_amps    = amps(:,randperm(length(amps)))';

shuffled_amps    = reshape([shuffled_amps' ; zeros(size(shuffled_amps'))],[],1);



4*16*40 / 60




((60 * 6) + (20 *5)) * 6 /60








