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
cfg.ignoreold_par_db            = true;

cfg.raw_dir                     = [dirs.rcs_pia, 'raw/'];
cfg.proc_dir                    = [dirs.rcs_pia, 'processed/'];

cfg.ephy_anal_dir               = [dirs.rcs_pia, '/ephy_analysis/'];

% specify patient hemispheres
%%% pts to update database from scratch locally:
%pt_sides        = {'RCS02R','RCS05R', 'RCS05L','RCS04R','RCS04L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};

%%% pts for NEJM submission:
%pt_sides           = {'RCS02R', 'RCS05R', 'RCS05L'};

%pt_sides           = {'RCS04R','RCS04L'};

pt_sides           = {'RCS02R','RCS05L','RCS05R','RCS07L', 'RCS07R'};

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
cfg.dates         = 'AllTime';
cfg.date_range    = {'03-Mar-2023', '01-Jul-2023'};

%%% return every aDBS ever tried (takes much longer):
%cfg.dates        = 'AllTime';

%%% state-current relationship (12 am - 12 pm)
cfg.plt_state_dur = 'sub_session_duration';
%%% state-current relationship (from 1-2 am and 1-2 pm):
%cfg.plt_state_dur = 'two_chunks'; 

%%% plot aDBS performance over months
% w/ aligned INS logs, plot requested dates
%pt_sides           = {'RCS02R','RCS05R','RCS05L'};
pts                = {'RCS02'};

for i = 1:length(pts)

    DBS_sum.(pts{i}) ...
        ...
        = plot_timeline_DBS(...
        ...
    cfg,    pts{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);
   
end



%% plotting adaptive DBS only
pt_sides           = {'RCS02R', 'RCS05R','RCS05L', 'RCS07R','RCS07L'};
cfg.dates         = 'DateRange';
cfg.date_range    = {'03-Mar-2023', '01-Jul-2023'};


for i = 1:length(pt_sides)

    aDBS_sum.(pt_sides{i}) ...
        ...
        = plot_longitudinal_aDBS(...
        ...
    cfg,    pt_sides{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);
end
%%
%% plot duty cycle versus pain
adbs_sum        = aDBS_sum.(pt_sides{i});
r_cap_tbl       = REDcap.(pts{i});

t_aDBS_sett = adbs_sum.timeStop_INS_log - adbs_sum.timeStart_INS_log ;



adbs_sum(le(t_aDBS_sett, duration('72:00:00')), :) = [];

pain_metrics = {'mayoNRS', 'painVAS', 'MPQtotal'};

med_vars = cellfun(@(x) sprintf('median_%s',x), pain_metrics, 'UniformOutput',false);
iqr_vars  = cellfun(@(x) sprintf('iqr_%s',x), pain_metrics, 'UniformOutput',false);
n_vars    = cellfun(@(x) sprintf('N_%s',x), pain_metrics, 'UniformOutput',false);


for i_aDBS = 1 : height(adbs_sum)

    i_rcap = find(...
                    isbetween(r_cap_tbl.time,...
                            adbs_sum.timeStart_INS_log(i_aDBS), ...
                            adbs_sum.timeStop_INS_log(i_aDBS)));

    adbs_sum{i_aDBS, med_vars} ...
        = median(r_cap_tbl{i_rcap, pain_metrics}, 'omitnan');


    adbs_sum{i_aDBS, iqr_vars} ...
        = iqr(r_cap_tbl{i_rcap, pain_metrics});
    
    adbs_sum{i_aDBS, n_vars} ...
        = length(r_cap_tbl{i_rcap, pain_metrics});


end

figure
errorbar(adbs_sum.avg_percent_on,...
         adbs_sum.median_mayoNRS, ...
         adbs_sum.iqr_mayoNRS/2, 'o', 'LineWidth',2);


ylim([0, 10]); xlim([0, 100]);

ylabel('mayoNRS');    xlabel('percent duty cycle');     

title(['median w/ IQR pain',newline,'for given aDBS duty cycle'])

set(gca, 'FontSize', 12)


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
cfg.proc_dir                 = [dirs.rcs_pia, 'processed/pain_per_DBS_parameters/'];

% sub-folder for specified stim params
cfg.proc_subdir       = 'by_freq_amp_pw_cyc';                             

%**note** plots are ranked by mean of FIRST pain metric specified below
cfg.plt_metrics       = {'mayoNRS', 'painVAS', 'MPQtotal', 'unpleasantVAS'};

cfg.exclude_VAS_50s.decision = true;
cfg.exclude_VAS_50s.plot     = true;

% see 'stimGroups' for variable names (do NOT include side--done automatically)
% to find all combinations of said stim parameters, and
% plot as seperate boxplots w/ clear labels
cfg.seperate_by               = {'rateInHz', 'pulseWidthInMicroseconds','ampInMilliamps',...  
                                 'cycleOnInSecs', 'cycleOffInSecs'};
% option to include the closed-loop settings (of a given contact) as
% furthest right boxplot (treats closed-loop setting as one group)
cfg.include_cl                = true;

cfg.pt_lbls.RCS02 = 'Patient 1';
cfg.pt_lbls.RCS05 = 'Patient 2';

pts = {'RCS02', 'RCS05'};

for i = 1:length(pts)
    plted_stim_groups.(pts{i})  = plot_stim_groups(cfg, dirs, pts{i}, stimGroups);
end
%%% repeat but keep duty cycle and frequency the same (group amplitude and
% pulse width together--returned as mean amplitude and pulse width for
% given frequency and duty cycle)

cfg.proc_subdir       = 'same_dutyCycle_and_freq';
cfg.seperate_by       = {'percentDutyCycle', 'rateInHz'};

for i = 1:length(pts)
    plted_stim_groups.(pts{i})  = plot_stim_groups(cfg,pts{i}, stimGroups);
end
%% parsimonious open and closed-loop programs for 2023 NEJM submission
cfg                          = [];
cfg.min_n_reports            = 30;
cfg.min_n_reports_subspace   = 5;

cfg.proc_dir                 = [dirs.rcs_pia, 'beh_analysis/pain_per_DBS_parameters/NEJM_2023_submission/'];
cfg.proc_subdir              = 'same_dutyCycle_and_freq';

%**note** plots are ranked by mean of FIRST pain metric specified below
cfg.plt_metrics              = {'mayoNRS', 'painVAS', 'MPQtotal', 'unpleasantVAS'};

cfg.exclude_VAS_50s.decision = true;
cfg.exclude_VAS_50s.plot     = false;

cfg.seperate_by              = {'percentDutyCycle', 'rateInHz'};
cfg.include_cl               = true;

cfg.pt_lbls.RCS02            = 'Patient 1';
cfg.pt_lbls.RCS05            = 'Patient 2';

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

    plted_stim_groups.(cfg.pt_id )= plot_stim_groups(cfg,'RCS04_July22_Apr23', stimGroups);




%
%
%
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

cfg.save_xlsx   = true;

cfg.CBDP_method = 'manual'; %cfg.CBDP_method = 'top_two';

cfg.clus_fp_plots = false;

cfg.source_dir  = ['/Users/Leriche/',...
                   'Dropbox (UCSF Department of Neurological Surgery)/',...
                   'UFlorida_UCSF_RCS_collab/Pain Reports/beh_clustered/'];

cfg.fig_dir       = ['/Users/Leriche/', ...
                    'Dropbox (UCSF Department of Neurological Surgery)/', ...
                    'UFlorida_UCSF_RCS_collab/beh_analysis/figs/beh_only/stages1_2_3'];

pts = {'RCS05'};

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

