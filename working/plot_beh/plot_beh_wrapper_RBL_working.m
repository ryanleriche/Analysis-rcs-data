%% load CONFIG_beh_wrapper

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
cfg.ignoreold_INS_logs          = false;
cfg.ignoreold_par_db            = false;

cfg.raw_dir                     = [dirs.rcs_pia, 'raw/'];
cfg.proc_dir                    = [dirs.rcs_pia, 'processed/'];

cfg.ephy_anal_dir               = [dirs.rcs_pia, '/ephy_analysis/aDBS_offline_sessions/'];

% specify patient hemispheres
%%% pts to update database from scratch locally:
%pt_sides        = {'RCS02R','RCS05R', 'RCS05L','RCS04R','RCS04L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};

%%% pts for NEJM submission:
%pt_sides           = {'RCS02R', 'RCS05R', 'RCS05L'};

%pt_sides           = {'RCS04R','RCS04L'};

pt_sides           = {'RCS02R','RCS05L','RCS05R'};

for i = 1: length(pt_sides)
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

%pt_sides           = {'RCS02R','RCS05R','RCS05L'};

pt_sides           = {'RCS05L','RCS05R'};


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
cfg.date_range    = {'28-Feb-2023'; '27-Apr-2023'};

%%% return every aDBS ever tried (takes much longer):
%cfg.dates        = 'AllTime';

%%% state-current relationship (12 am - 12 pm)
cfg.plt_state_dur = 'sub_session_duration';

%%% state-current relationship (from 1-2 am and 1-2 pm):
%cfg.plt_state_dur = 'two_chunks'; 

%%% plot aDBS performance over months
% w/ aligned INS logs, plot requested dates
%pt_sides           = {'RCS02R','RCS05R','RCS05L'};


pts                = {'RCS05'};

for i = 1:length(pts)

    DBS_sum.(pts{i}) ...
        ...
        = plot_timeline_DBS(...
        ...
    cfg,    pts{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);
end











%%
pt_sides           = {'RCS05L'};

for i = 1:length(pt_sides)

    aDBS_sum.(pt_sides{i}) ...
        ...
        = plot_longitudinal_aDBS(...
        ...
    cfg,    pt_sides{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);
end


%% unilateral implant AND used INS logs to capture PTM intiated group changes
%{
[REDcap_INSLog.RCS02R, proc_group_changes.RCS02R]  ...
    ...
    = explicit_INS_stim_params(...
    ...
cfg, INS_logs.RCS02R.group_changes, db.RCS02R, REDcap.RCS02, visits.RCS02);



%% add stim params to REDcap based off of EventLog.txt, and DeviceSettings.json files
cfg                     = [];
cfg.stage_dates         = stage_dates{2};
cfg.pt_id               = 'RCS02';

INSLog_w_redcap.RCS02R  = get_INSLog_stim_params(cfg,...
                                     INS_logs.RCS02R.group_changes,...
                                     db.RCS02R,...
                                     REDcap.RCS02,...
                                     visits.RCS02...
                                    );

[wrt_stim_REDcap.RCS02, stimGroups.RCS02] ...
    ...
    = make_stim_groups(...
    ...
'RCS02', [], INSLog_w_redcap.RCS02R, visits.RCS02);


cfg                          = [];
cfg.pt_id                    ='RCS02';
cfg.min_n_reports            = 5;

    plted_stim_groups.RCS02  = plot_stim_groups(cfg, stimGroups.RCS02);

%}
%% generate box plots of pain metrics wrt stim parameters
%{

input(s)

___________________________________________________________________________

* w/ aligned REDcap reports, first visualize the parameter space per each
   "major contact pair" (best analgesic contact pairs as assessed clinically--
    heuristically those contacts w/ a 10+ pain reports)

* box plots of NRS, VAS, MPQ, and PC1 (for RCS05 only) based of
  parsimonious amp-PW-freq space *per* contact


___________________________________________________________________________

output(s)
wrt_stim_REDcap.RCS0X | combined left and right stim params PER REDcap survey
stimGroups.RCS0X      | surveys grrouped per major contact pair, stim OFF, @ 0mA, etc


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

% RCS04, RCS05, RCS06, and RCS07 can be handled together
%pts = {'RCS04', 'RCS05', 'RCS06', 'RCS07'};

pts = {'RCS02', 'RCS05'};

%pts = {'RCS02','RCS04', 'RCS05','RCS07'};

for i = 1:length(pts)

    switch pts{i}
        case 'RCS02'
            [wrt_stim_REDcap.(pts{i}), stimGroups.(pts{i})] ...
            ...
            = make_stim_groups(...
            ...
            pts{i}, [], REDcap.([pts{i}, 'R']), pt_META.(pts{i}));

        otherwise
            [wrt_stim_REDcap.(pts{i}), stimGroups.(pts{i})] ...
            ...
            = make_stim_groups(...
            ...
            pts{i}, REDcap.([pts{i}, 'L']), REDcap.([pts{i}, 'R']), pt_META.(pts{i}));
    end
    % add pain fluctuation study as "stim_group"
    stimGroups.(pts{i}).('Pre-trial baseline') = {PFS.(pts{i})};
end
%% return stim boxplots based on specified grouping of stim parameters
cfg                          = [];
cfg.min_n_reports            = 10;
cfg.min_n_reports_subspace   = 2;

% general folder for 'by_contacts'
cfg.proc_dir          = [dirs.rcs_pia, 'processed/pain_per_DBS_parameters/'];
% sub-folder for specified stim params
cfg.proc_subdir       = 'by_freq_amp_pw_cyc';                             

% pain metrics to plot (stim boxplots in seperate folders)
cfg.plt_metrics       = {'painVAS', 'mayoNRS',  'MPQtotal', 'unpleasantVAS'};


cfg.exclude_VAS_50s.decision = true;
cfg.exclude_VAS_50s.plot     = true;


% see 'stimGroups' for variable names (do NOT include side--done automatically)
% to find all combinations of said stim parameters, and
% plot as seperate boxplots w/ clear labels
cfg.seperate_by               = {'rateInHz', 'pulseWidthInMicroseconds','ampInMilliamps',...  
                                 'cycleOnInSecs', 'cycleOffInSecs'};
cfg.include_cl                = true;

cfg.pt_lbls.RCS02 = 'Patient 1';
cfg.pt_lbls.RCS05 = 'Patient 2';

% option to include the closed-loop settings (of a given contact) as
% furthest right boxplot (treats closed-loop setting as one group)
cfg.include_cl        = true;               

set(0,'DefaultFigureVisible','off');            close all

%pts = {'RCS02', 'RCS04','RCS05', 'RCS06','RCS07'};

pts = {'RCS02', 'RCS05'};

for i = 1:length(pts)
    plted_stim_groups.(pts{i})  = plot_stim_groups(cfg,pts{i}, stimGroups);
end
%%
%%% repeat but keep duty cycle and frequency the same (group amplitude and
%%% pulse width together--returned as mean amplitude and pulse width for
%%% given frequency and duty cycle)
cfg                          = [];
cfg.min_n_reports            = 10;
cfg.min_n_reports_subspace   = 2;


cfg.proc_dir          = [dirs.rcs_pia, 'processed/pain_per_DBS_parameters/'];
cfg.proc_subdir       = 'same_dutyCycle_and_freq';

%**note** plots are ranked by mean of FIRST pain metric specified below
cfg.plt_metrics       = {'mayoNRS', 'painVAS', 'MPQtotal', 'unpleasantVAS'};

cfg.exclude_VAS_50s.decision = true;
cfg.exclude_VAS_50s.plot     = true;


cfg.seperate_by              = {'percentDutyCycle', 'rateInHz'};
cfg.include_cl               = true;

cfg.pt_lbls.RCS02 = 'Patient 1';
cfg.pt_lbls.RCS05 = 'Patient 2';

set(0,'DefaultFigureVisible','off');         close all

pts = {'RCS02','RCS05'};

for i = 1:length(pts)
    plted_stim_groups.(pts{i})  = plot_stim_groups(cfg,pts{i}, stimGroups);
end
%% parsimonious open and closed-loop programs for 2023 NEJM submission
cfg                          = [];
cfg.min_n_reports            = 30;
cfg.min_n_reports_subspace   = 5;

cfg.proc_dir          = [dirs.rcs_pia, 'beh_analysis/pain_per_DBS_parameters/NEJM_2023_submission/'];
cfg.proc_subdir       = 'same_dutyCycle_and_freq';

%**note** plots are ranked by mean of FIRST pain metric specified below
cfg.plt_metrics       = {'mayoNRS', 'painVAS', 'MPQtotal', 'unpleasantVAS'};

cfg.exclude_VAS_50s.decision = true;
cfg.exclude_VAS_50s.plot     = false;

cfg.seperate_by               = {'percentDutyCycle', 'rateInHz'};
cfg.include_cl                = true;

cfg.pt_lbls.RCS02 = 'Patient 1';
cfg.pt_lbls.RCS05 = 'Patient 2';

set(0,'DefaultFigureVisible','off');         close all

pts = {'RCS02','RCS05'};

for i = 1:length(pts)
    
    grps  = stimGroups.(pts{i}).Properties.VariableNames;

    switch pts{i}
        case 'RCS02'
            rmv_grps = {'s3', 'visits','s2','RACC c+2-', 'RACC c+1-'};

        case 'RCS05'
            rmv_grps = {'s3', 'visits','s2','LCaud 0+2-', 'RThal c+3-'};

    end
    i_grp = ~contains(grps, rmv_grps);

    stimGroups_parsi.(pts{i})   = stimGroups.(pts{i})(:, i_grp);

    %
    [plted_stim_groups.(pts{i}),...
     plted_stimGroups_parsi_subspace.(pts{i})]...
        ...
        = plot_stim_groups(cfg, pts{i}, stimGroups_parsi);
end


%%
%%%
%% assess blinded testing
%{
04/15/23

Stage 3 design

    * 6 days on either sham, ol-stim, or cl-stim
    * 1 day washout
    * repeat for N weeks


Determining sham, ol-stim, versus cl-stim


sham:      both sides set to 0 mA Groups (A, B, or C)

ol-stim:   one side at >0mA and other side at 0 mA (Note both sides should still be ON)

cl-stim:   one side at Group D On

washout:   both sides OFF


%}

%% for RCS02's Stage 3, have own boxplots delineating open-loop programs versus sham
% as well as the "day-by-day" showing each stim parameters over given sets
% of days
cfg                   = [];
cfg.proc_dir          = [dirs.rcs_pia, 'beh_analysis/pain_per_DBS_parameters/NEJM_2023_submission/'];
cfg.proc_subdir       = 'Stage 3 (ol-DBS_versus_sham)';


cfg.plt_metrics       = {'mayoNRS','painVAS',  'MPQtotal', 'unpleasantVAS'};
cfg.pt_id             = 'RCS02';

cfg.pt_meta           = pt_META;
cfg.N_days            = duration('3:00:00:00');

cfg.seperate_by       = {'R_stim_groups'};

close all
set(0,'DefaultFigureVisible','off')
s3.RCS02 = struct;

[s3.RCS02.groups, s3.RCS02.sham_vs_stim_stats]...
    ...
    = plot_RCS02_stage3_v2(...
    ...
cfg, stimGroups.RCS02.s3{1}, stimLog.RCS02R);
%% seperate out RCS05's blinded testing

cfg                   = [];
cfg.proc_dir          = [dirs.rcs_pia, 'beh_analysis/pain_per_DBS_parameters/NEJM_2023_submission/'];
cfg.proc_subdir       = '/blinded testing/';


cfg.plt_metrics       = {'mayoNRS','painVAS',  'MPQtotal', 'unpleasantVAS'};
pt_id                 = 'RCS05';

cfg.pt_meta           = pt_META;


s2_blind.RCS05        = plot_RCS05_blinded(cfg, pt_id, stimGroups, INS_logs_API_t_synced);




%% RCS04 --> stim groups after starting buprenorphine (after July 2022 home visit)
cfg                   = [];
cfg.pt_id             = 'RCS04_July22_Apr23';


i_epoch          = find(ge(REDcap.RCS04.time, pt_META.RCS04.dates(end) + duration('24:00:00')));

[~, stimGroups.(cfg.pt_id)] ...
    ...
    = make_stim_groups(...
    ...
 'RCS04', REDcap.RCS04L(i_epoch : end,:), REDcap.RCS04R(i_epoch: end, :), pt_META.RCS04);

 stimGroups.(cfg.pt_id).('Pre-trial baseline') = {fluct.('RCS04')};
%%
%cfg                          = [];
cfg.min_n_reports            = 3;

cfg.min_n_reports_subspace   = 2;

cfg.proc_dir          = [dirs.rcs_pia, '/beh_analysis/pain_per_DBS_parameters/'];
cfg.proc_subdir       = 'same_dutyCycle_and_freq';

cfg.plt_metrics       = {'painVAS', 'mayoNRS',  'MPQtotal', 'unpleasantVAS'};

cfg.seperate_by       = {'percentDutyCycle', 'rateInHz'};

cfg.include_cl        = false;

cfg.exclude_VAS_50s.decision   = false;
cfg.exclude_VAS_50s.plot       = true;

cfg.pt_lbls.RCS04_July22_Apr23 = 'RCS04_July22_Apr23';

set(0,'DefaultFigureVisible','off');         close all

    plted_stim_groups.(cfg.pt_id )= plot_stim_groups(cfg,'RCS04_July22_Apr23', stimGroups);
%%
%%
%%

%%



%%
%% distributions of pain metrics and relationship BETWEEN pain metrics
pts = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};

cfg                     = [];
cfg.dates               = 'AllTime';

for i = 4%1:length(pts)
    cfg.pt_id  = pts{i};       

    %plot_hist(cfg, REDcap);          
    plot_versus(cfg, REDcap);
end

%% RCS02: explore NRS and VAS "mismatch" 
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
pts = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};

pain_space      = [];

cfg             = [];     
cfg.dates       = 'AllTime';
cfg.pca         = false;
cfg.plt_VAS     = true;
cfg.VAS_only    = false;

cfg.save_xlsx   = false;

cfg.CBDP_method = 'top_two';

cfg.source_dir  = ['/Users/Leriche/',...
                   'Dropbox (UCSF Department of Neurological Surgery)/',...
                   'UFlorida_UCSF_RCS_collab/Pain Reports/beh_clustered/'];

cfg.fig_dir       = ['/Users/Leriche/', ...
                    'Dropbox (UCSF Department of Neurological Surgery)/', ...
                    'UFlorida_UCSF_RCS_collab/beh_analysis/figs/beh_only/stages1_2_3'];

set(0,'DefaultFigureVisible','off')

for i =  1:length(pts)

    cfg.pt_id  = pts{i};

    [pain_space.(pts{i})] = plot_pain_space(cfg, REDcap);
end

set(0,'DefaultFigureVisible','on')
%%
%
%% ******************** in-progress versions below ******************** %%
%
%%
% last N days for: 
close all
cfg                     = [];

cfg.pt_id               = 'RCS04';
cfg.dates               = 'PreviousDays';
cfg.ndays               = 10;

cfg.subplot             = true;
cfg.sum_stat_txt        = true;
cfg.stim_parameter      = '';

    plot_timeline(cfg, REDcap, PFS_sum_stats);
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

% visualize distribution of PC 1--wait for now--could be useful for distributed closed-loop


%% parse db to human-readable format
i                   = cellfun(@(x) length(x) == 1, db.RCS02R.duration);
db_RCS02R           = db.RCS02R(i, :);

[~, i_u] = unique(db_RCS02R.sess_name);

db_RCS02R = db_RCS02R(i_u, :);

db_RCS02R.timeStart = cellfun(@(x) x, db_RCS02R.timeStart);
db_RCS02R.timeStop  = cellfun(@(x) x, db_RCS02R.timeStop);
db_RCS02R.duration  = cellfun(@(x) x, db_RCS02R.duration);

% take sessions of useful yet managable duration
i_sess              = ge(db_RCS02R.duration , duration('00:07:00')) &...
                      le(db_RCS02R.duration , duration('01:00:00'));

db_RCS02R           = db_RCS02R(i_sess, :);

i_stim              = cellfun(@(x) ~isempty(x), db_RCS02R.stimSettingsOut);

db_RCS02R.activeGroup(i_stim)  ...
    = cellfun(@(x) x.activeGroup, db_RCS02R.stimSettingsOut(i_stim));


db_RCS02R  = sortrows(db_RCS02R, 'timeStart', 'descend');

% try to simulate certain sessions
sess_oi    = [2, 15, 18:19];
db_RCS02R  = db_RCS02R(sess_oi,:);

nickname = {'i'; 'ii';'iii';'iv'};
db_RCS02R.sess_name = cellfun(@(x,y) [x,'_',y], db_RCS02R.sess_name, nickname, 'UniformOutput', false);

db_RCS02R.per_TD_lost = nan(height(db_RCS02R),1);

for i_sess = 1 : height(db_RCS02R)

    data_dir = db_RCS02R.path{i_sess};
    
    % data using Analysis-rcs-data
    
    [unifiedDerivedTimes,...
        timeDomainData, ~, ~,...
        ~, ~, ~, ...
        PowerData, ~, ~,...
        FFTData, ~, ~,...
        ...
        AdaptiveData, ~, ~, timeDomainSettings, powerSettings,...
        ~, eventLogTable, metaData, stimSettingsOut, stimMetaData,...
        stimLogSettings, DetectorSettings, AdaptiveStimSettings, ...
        AdaptiveEmbeddedRuns_StimSettings, ~] ...
        ...
        = ProcessRCS(data_dir, 3);
    
    dataStreams         = {timeDomainData, PowerData, AdaptiveData, FFTData};

    comb_dt = createCombinedTable(dataStreams, unifiedDerivedTimes, metaData);
    
    [comb_dt_chunks, per_TD_lost]    = chunks_and_gaps(comb_dt);

    db_RCS02R.per_TD_lost(i_sess)    = per_TD_lost;
    db_RCS02R.comb_dt_chunks(i_sess) = comb_dt_chunks;
end
%% simulate LD activity 
set(0,'DefaultFigureVisible','off')
for i_sess = 1 : height(db_RCS02R)

    sim_tbl     = td_to_fft_pb(i_sess, db_RCS02R);
end
set(0,'DefaultFigureVisible','on')
%%
% pull out all sessions w/ FFT channel streamed

fftSettings  = vertcat(db_RCS02R.fftSettings{:});
fftSettings  = struct2table(fftSettings.fftConfig);

% length of 1 implies that it is a channel name rather than 'Disabled'
i_fft_stream = find(cellfun(@(x)  length(x.fftConfig.fftStreamChannel) == 1, ...
                        db_RCS02R.fftSettings));

i_sess       = i_fft_stream(1);

td_to_fft_pb(i_sess, db_RCS02R);

%%

% need aDBS session to see actual LD outputs as validation
% --> use most recent one

i_sess                 = find(strcmp(db_RCS02R.activeGroup, 'D'));
i_sess                 = i_sess(2);



fprintf(['aDBS %s | starting at %s | duration %s (HH:MM:SS.sss)', newline], ...
    db_RCS02R.sess_name{i_sess}, ...
    db_RCS02R.timeStart(i_sess),...
    db_RCS02R.duration(i_sess));

%% animate pain space visualization

%v = RecordRotation();

%
%{
playRecording(v);
%%
function V2 = RecordRotation()

n = 300;
V = zeros(n,2);
T = timer('period',0.001,'executionmode','fixedrate',...
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



% n_runs          = 5;            t_off           = 20;
% n_iter          = 2;            t_on            = 60;
% 
% 
% t_in_sec        = n_iter * n_runs * (t_off+ t_on);
% duty_cycle      = t_on / (t_off + t_on);
% 
% t_in_min        = t_in_sec / 60;
% 


%% old code:
% (previously used DeviceSettings.json alone to get stim params)
% --> now likely redundant and less accurate than
% 'align_REDcap_to_stimLog()' fxn

    % stim parameters grouping based off of 'DeviceSettings.json' from 'StimSettingsOut' var in db.RCSXXX
    %{
    %{
    * creates db_beh.RCSXX with timestamps, contacts, amp, PW, freq, cycling,
      group, and laterality (unilateral vs bilateral stim)
    
        -manually inspect laterality to handle edge cases
    
    * align nearest REDcap report as row in db_beh.RCSXX
    
    
    for RCS05 and RCS02 especially look at pain reports within session 5 - 15 minutes
       
    
    %}
    
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


    % validate StimLog.json's w/ db.RCSXXX 'stimSettingsOut' (from DeviceSettings.json)

    %  used for sanity check btwn StimLog.json and DeviceSettings.json)
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
