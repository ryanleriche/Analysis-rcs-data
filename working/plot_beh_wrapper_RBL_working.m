% see CONFIG_plot_beh.m script to specify directories, pull REDcap, 
% and generate subdirectories
CONFIG_plot_beh;

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
%}

cfg_rcs_db.load_EventLog      = true;

% option to load previous database for efficient processing
cfg_rcs_db.ignoreold_db                = false;
cfg_rcs_db.ignoreold_INS_logs          = false;
cfg_rcs_db.ignoreold_par_db            = true;

% specify patient hemispheres
%%% pts to update database from scratch locally:
%pt_sides        = {'RCS02R','RCS05R', 'RCS05L','RCS04R','RCS04L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};
pt_sides           = {'RCS02R', 'RCS07R','RCS04L'};

pt_sides           = {'RCS05L'};
for i = 1  : length(pt_sides)
    %%% process RCS .jsons into searchable database
    [db.(pt_sides{i}), bs.(pt_sides{i})] ...
        ...
        = makeDatabaseRCS_Ryan(...
        ...
        cfg_rcs_db, pt_sides{i});

    %%% process INS logs .txts based on unique entries only
        % (INS logs have mostly repeating entries)
    INS_logs.(pt_sides{i})  ...
        ...
        = RCS_logs( ...
        ...
        cfg_rcs_db, pt_sides{i});


    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings
    [par_db.(pt_sides{i}), ss_var_oi.(pt_sides{i})] ...
        ...
        = makeParsedDatabaseRCS(...
        ...
        cfg_rcs_db, pt_sides{i}, db);

end
%%
for i = 1: length(pt_sides)
    %%% find nearest (yet, preceding) streaming session to INS log entry
    % --> accounts for INS to API time latency
    [par_db_aDBS_ss.(pt_sides{i}), INS_logs_API_t_synced.(pt_sides{i})] ...
    ...
        = align_INSLogs_to_API_time(...
    ...
    pt_sides{i}, INS_logs, par_db, ss_var_oi);
end

%%% plotting adaptive DBS only
cfg_rcs_db.dates         = 'DateRange';
cfg_rcs_db.date_range    = {'03-Mar-2023', '01-Jul-2023'};
cfg_rcs_db.plt_state_dur = 'sub_session_duration';

for i = 1:length(pt_sides)

    aDBS_sum.(pt_sides{i}) ...
        ...
        = plot_longitudinal_aDBS(...
        ...
    cfg_rcs_db,    pt_sides{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);

    
end

 
%%
%%% specify which dates to return:
cfg_rcs_db.dates         = 'DateRange';
cfg_rcs_db.date_range    = {'20-Apr-2023', '01-Jul-2023'};

%%% return every aDBS ever tried (takes much longer):
%cfg.dates        = 'AllTime';

%%% state-current relationship (12 am - 12 pm)
cfg_rcs_db.plt_state_dur = 'sub_session_duration';

%%% state-current relationship (from 1-2 am and 1-2 pm):
%cfg.plt_state_dur = 'two_chunks'; 

%%% plot aDBS performance over months
% w/ aligned INS logs, plot requested dates
%pt_sides           = {'RCS02R','RCS05R','RCS05L'};
pts                = {'RCS02'};

for i = 1:length(pts)

    DBS_sum.(pts{i}) ...
        ...
        = plot_longitudinal_all_DBS(...
        ...
    cfg_rcs_db,    pts{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);
   
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


pt_sides           = {'RCS02R','RCS05L','RCS05R','RCS07L', 'RCS07R'};

for i = 1  : length(pt_sides)
    %%% finally merge stim parameters to REDcap survey during said parameters
    [stimLog.(pt_sides{i}), REDcap.(pt_sides{i})]  ...
        ...
        = align_REDcap_to_stimLog(...
        ...
        cfg_rcs_db, pt_sides{i}, db, REDcap);
end

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
cfg_rcs_db                          = [];
cfg_rcs_db.min_n_reports            = 10;
cfg_rcs_db.min_n_reports_subspace   = 2;
cfg_rcs_db.proc_dir                 = [dirs.rcs_pia, 'processed/pain_per_DBS_parameters/'];

% sub-folder for specified stim params
cfg_rcs_db.proc_subdir       = 'by_freq_amp_pw_cyc';                             

%**note** plots are ranked by mean of FIRST pain metric specified below
cfg_rcs_db.plt_metrics       = {'mayoNRS', 'painVAS', 'MPQtotal', 'unpleasantVAS'};

cfg_rcs_db.exclude_VAS_50s.decision = true;
cfg_rcs_db.exclude_VAS_50s.plot     = true;

% see 'stimGroups' for variable names (do NOT include side--done automatically)
% to find all combinations of said stim parameters, and
% plot as seperate boxplots w/ clear labels
cfg_rcs_db.seperate_by               = {'rateInHz', 'pulseWidthInMicroseconds','ampInMilliamps',...  
                                 'cycleOnInSecs', 'cycleOffInSecs'};
% option to include the closed-loop settings (of a given contact) as
% furthest right boxplot (treats closed-loop setting as one group)
cfg_rcs_db.include_cl                = true;

cfg_rcs_db.pt_lbls.RCS02 = 'Patient 1';
cfg_rcs_db.pt_lbls.RCS05 = 'Patient 2';

pts = {'RCS02', 'RCS05'};

for i = 1:length(pts)
    plted_stim_groups.(pts{i})  = plot_stim_groups(cfg_rcs_db, dirs, pts{i}, stimGroups);
end
%%% repeat but keep duty cycle and frequency the same (group amplitude and
% pulse width together--returned as mean amplitude and pulse width for
% given frequency and duty cycle)

cfg_rcs_db.proc_subdir       = 'same_dutyCycle_and_freq';
cfg_rcs_db.seperate_by       = {'percentDutyCycle', 'rateInHz'};

for i = 1:length(pts)
    plted_stim_groups.(pts{i})  = plot_stim_groups(cfg_rcs_db,pts{i}, stimGroups);
end
%% parsimonious open and closed-loop programs for 2023 NEJM submission
cfg_rcs_db                          = [];
cfg_rcs_db.min_n_reports            = 30;
cfg_rcs_db.min_n_reports_subspace   = 5;

cfg_rcs_db.proc_dir                 = [dirs.rcs_pia, 'beh_analysis/pain_per_DBS_parameters/NEJM_2023_submission/'];
cfg_rcs_db.proc_subdir              = 'same_dutyCycle_and_freq';

%**note** plots are ranked by mean of FIRST pain metric specified below
cfg_rcs_db.plt_metrics              = {'mayoNRS', 'painVAS', 'MPQtotal', 'unpleasantVAS'};

cfg_rcs_db.exclude_VAS_50s.decision = true;
cfg_rcs_db.exclude_VAS_50s.plot     = false;

cfg_rcs_db.seperate_by              = {'percentDutyCycle', 'rateInHz'};
cfg_rcs_db.include_cl               = true;

cfg_rcs_db.pt_lbls.RCS02            = 'Patient 1';
cfg_rcs_db.pt_lbls.RCS05            = 'Patient 2';

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
        = plot_stim_groups(cfg_rcs_db, pts{i}, stimGroups_parsi);
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
cfg_rcs_db                   = [];
cfg_rcs_db.proc_dir          = [dirs.rcs_pia, 'beh_analysis/pain_per_DBS_parameters/NEJM_2023_submission/'];
cfg_rcs_db.proc_subdir       = 'Stage 3 (ol-DBS_versus_sham)';


cfg_rcs_db.plt_metrics       = {'mayoNRS','painVAS',  'MPQtotal', 'unpleasantVAS'};
cfg_rcs_db.pt_id             = 'RCS02';

cfg_rcs_db.pt_meta           = pt_META;
cfg_rcs_db.N_days            = duration('3:00:00:00');

cfg_rcs_db.seperate_by       = {'R_stim_groups'};

close all
set(0,'DefaultFigureVisible','off')
s3.RCS02 = struct;

[s3.RCS02.groups, s3.RCS02.sham_vs_stim_stats]...
    ...
    = plot_RCS02_stage3_v2(...
    ...
cfg_rcs_db, stimGroups.RCS02.s3{1}, stimLog.RCS02R);
%% seperate out RCS05's blinded testing

cfg_rcs_db                   = [];
cfg_rcs_db.proc_dir          = [dirs.rcs_pia, 'beh_analysis/pain_per_DBS_parameters/NEJM_2023_submission/'];
cfg_rcs_db.proc_subdir       = '/blinded testing/';


cfg_rcs_db.plt_metrics       = {'mayoNRS','painVAS',  'MPQtotal', 'unpleasantVAS'};
pt_id                 = 'RCS05';

cfg_rcs_db.pt_meta           = pt_META;


s2_blind.RCS05        = plot_RCS05_blinded(cfg_rcs_db, pt_id, stimGroups, INS_logs_API_t_synced);

%% RCS04 --> stim groups after starting buprenorphine (after July 2022 home visit)
cfg_rcs_db                   = [];
cfg_rcs_db.pt_id             = 'RCS04_July22_Apr23';


i_epoch          = find(ge(REDcap.RCS04.time, pt_META.RCS04.dates(end) + duration('24:00:00')));

[~, stimGroups.(cfg_rcs_db.pt_id)] ...
    ...
    = make_stim_groups(...
    ...
 'RCS04', REDcap.RCS04L(i_epoch : end,:), REDcap.RCS04R(i_epoch: end, :), pt_META.RCS04);

 stimGroups.(cfg_rcs_db.pt_id).('Pre-trial baseline') = {fluct.('RCS04')};
%%
%cfg                          = [];
cfg_rcs_db.min_n_reports            = 3;

cfg_rcs_db.min_n_reports_subspace   = 2;

cfg_rcs_db.proc_dir          = [dirs.rcs_pia, '/beh_analysis/pain_per_DBS_parameters/'];
cfg_rcs_db.proc_subdir       = 'same_dutyCycle_and_freq';

cfg_rcs_db.plt_metrics       = {'painVAS', 'mayoNRS',  'MPQtotal', 'unpleasantVAS'};

cfg_rcs_db.seperate_by       = {'percentDutyCycle', 'rateInHz'};

cfg_rcs_db.include_cl        = false;

cfg_rcs_db.exclude_VAS_50s.decision   = false;
cfg_rcs_db.exclude_VAS_50s.plot       = true;

cfg_rcs_db.pt_lbls.RCS04_July22_Apr23 = 'RCS04_July22_Apr23';

    plted_stim_groups.(cfg_rcs_db.pt_id )= plot_stim_groups(cfg_rcs_db,'RCS04_July22_Apr23', stimGroups);




%
%
%
%% distributions of pain metrics and relationship BETWEEN pain metrics
pts = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};

cfg_rcs_db                     = [];
cfg_rcs_db.dates               = 'AllTime';

for i = 4%1:length(pts)
    cfg_rcs_db.pt_id  = pts{i};       

    %plot_hist(cfg, REDcap);          
    plot_versus(cfg_rcs_db, REDcap);
end

%% RCS02: explore NRS and VAS "mismatch" 
cfg_rcs_db.pt_id               = 'RCS02-mismatch';
cfg_rcs_db.stage_dates         = stage_dates{2}; % starts at Stage 1

      plot_timeline(cfg_rcs_db, REDcap);

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

cfg_rcs_db             = [];     
cfg_rcs_db.dates       = 'AllTime';
cfg_rcs_db.pca         = false;
cfg_rcs_db.plt_VAS     = true;
cfg_rcs_db.VAS_only    = false;

cfg_rcs_db.save_xlsx   = true;

cfg_rcs_db.CBDP_method = 'manual'; %cfg.CBDP_method = 'top_two';

cfg_rcs_db.clus_fp_plots = false;

cfg_rcs_db.source_dir  = ['/Users/Leriche/',...
                   'Dropbox (UCSF Department of Neurological Surgery)/',...
                   'UFlorida_UCSF_RCS_collab/Pain Reports/beh_clustered/'];

cfg_rcs_db.fig_dir       = ['/Users/Leriche/', ...
                    'Dropbox (UCSF Department of Neurological Surgery)/', ...
                    'UFlorida_UCSF_RCS_collab/beh_analysis/figs/beh_only/stages1_2_3'];

pts = {'RCS05'};

for i =  1:length(pts)

    cfg_rcs_db.pt_id  = pts{i};

    [pain_space.(pts{i})] = plot_pain_space(cfg_rcs_db, REDcap);
end


%%
%
%% ******************** in-progress versions below ******************** %%
%
%%
% last N days for: 
close all; set(0,'DefaultFigureVisible','on')
cfg_rcap                     = [];

cfg_rcap.pt_id               = 'RCS04';
cfg_rcap.dates               = 'PreviousDays';
cfg_rcap.ndays               = 7;

cfg_rcap.stage_dates         = stage_dates;

cfg_rcap.subplot             = true;
cfg_rcap.sum_stat_txt        = true;
cfg_rcap.stim_parameter      = '';

    plot_timeline(cfg_rcap, REDcap, PFS_sum_stats);