%% see user_local_CONFIG.m script to specify directories and general set-up
user_local_CONFIG;

% load REDcap pain surveys and set-up directories
general_setup;

% see "dirs" struct for main directories used, and "rcs_db_cfg" for the
% configuration of the RC+S databasing

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

% option to load previous database for efficient processing
rcs_db_cfg.ignoreold_db                = false;
rcs_db_cfg.ignoreold_INS_logs          = false;
rcs_db_cfg.ignoreold_par_db            = false;

% specify patient hemispheres
%%% pts to update database from scratch locally:
%pt_sides        = {'RCS02R','RCS05R', 'RCS05L','RCS04R','RCS04L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};
%pt_sides           = {'RCS02R', 'RCS07R','RCS04L'};

pt_sides           = {'RCS05R','RCS05L'};
for i = 1  : length(pt_sides)
    %%% process RCS .jsons into searchable database
    [db.(pt_sides{i}), bs.(pt_sides{i})] ...
        ...
        = makeDatabaseRCS_Ryan(...
        ...
        rcs_db_cfg, pt_sides{i});

    %%% process INS logs .txts based on unique entries only
        % (INS logs have mostly repeating entries)
    INS_logs.(pt_sides{i})  ...
        ...
        = RCS_logs( ...
        ...
        rcs_db_cfg, pt_sides{i});


    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings
    [par_db.(pt_sides{i}), ss_var_oi.(pt_sides{i})] ...
        ...
        = makeParsedDatabaseRCS(...
        ...
        rcs_db_cfg, pt_sides{i}, db);


end
%%
%%% specify which dates to return:
rcs_db_cfg.dates         = 'DateRange';
rcs_db_cfg.date_range    = {'20-Apr-2023', '01-Jul-2023'};

%%% return every aDBS ever tried (takes much longer):
%cfg.dates        = 'AllTime';

%%% state-current relationship (12 am - 12 pm)
rcs_db_cfg.plt_state_dur = 'sub_session_duration';

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
    rcs_db_cfg,    pts{i},    REDcap,     INS_logs_API_t_synced,      par_db_aDBS_ss);
   
end

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
%pts = {'RCS02','RCS04', 'RCS05','RCS07'};

pts = {'RCS05'};
for i = 1  : length(pt_sides)
    %%% finally merge stim parameters to REDcap survey during said parameters
    [stimLog.(pt_sides{i}), REDcap.(pt_sides{i})]  ...
        ...
        = align_REDcap_to_stimLog(...
        ...
        rcs_db_cfg, pt_sides{i}, db, REDcap);
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
            pts{i}, REDcap.([pts{i}, 'L']), REDcap.([pts{i}, 'R']), pt_meta.(pts{i}));
    end
    % add pain fluctuation study as "stim_group"
    stimGroups.(pts{i}).('Pre-trial baseline') = {PFS.(pts{i})};
end
%% return stim boxplots based on specified grouping of stim parameters
stim_cfg                          = [];
stim_cfg.min_n_reports            = 10;
stim_cfg.min_n_reports_subspace   = 2;
stim_cfg.proc_dir                 = [dirs.rcs_pia, 'beh_analysis/pain_per_DBS_parameters/'];

% sub-folder for specified stim params
stim_cfg.proc_subdir       = 'by_freq_amp_pw_cyc';                             

%**note** plots are ranked by mean of FIRST pain metric specified below
stim_cfg.plt_metrics       = {'mayoNRS', 'painVAS', 'MPQtotal', 'unpleasantVAS'};

stim_cfg.exclude_VAS_50s.decision = true;
stim_cfg.exclude_VAS_50s.plot     = true;

% see 'stimGroups' for variable names (do NOT include side--done automatically)
% to find all combinations of said stim parameters, and
% plot as seperate boxplots w/ clear labels
stim_cfg.seperate_by               = {'rateInHz', 'pulseWidthInMicroseconds','ampInMilliamps',...  
                                 'cycleOnInSecs', 'cycleOffInSecs'};
% option to include the closed-loop settings (of a given contact) as
% furthest right boxplot (treats closed-loop setting as one group)
stim_cfg.include_cl                = true;

stim_cfg.pt_lbls.RCS02 = 'Patient 1';
stim_cfg.pt_lbls.RCS05 = 'RCS05';
stim_cfg.pt_lbls.RCS06 = 'Patient RCS06';

pts = {'RCS05'};

for i = 1:length(pts)
    plted_stim_groups.(pts{i})  = plot_stim_groups(stim_cfg, pts{i}, stimGroups);
end
%%% repeat but keep duty cycle and frequency the same (group amplitude and
% pulse width together--returned as mean amplitude and pulse width for
% given frequency and duty cycle)

stim_cfg.proc_subdir       = 'same_cycling_and_freq';
stim_cfg.seperate_by       = {'cycleOnInSecs', 'cycleOffInSecs', 'rateInHz'};

for i = 1:length(pts)
    plted_stim_groups.(pts{i})  = plot_stim_groups(stim_cfg, pts{i}, stimGroups);
end


stim_cfg.proc_subdir       = 'same_dutyCycle_and_freq';
stim_cfg.seperate_by       = {'percentDutyCycle', 'rateInHz'};

for i = 1:length(pts)
    plted_stim_groups.(pts{i})  = plot_stim_groups(stim_cfg, pts{i}, stimGroups);
end
%% parsimonious open and closed-loop programs for 2023 NEJM submission
stim_cfg                          = [];
stim_cfg.min_n_reports            = 30;
stim_cfg.min_n_reports_subspace   = 5;

stim_cfg.proc_dir                 = [dirs.rcs_pia, 'beh_analysis/pain_per_DBS_parameters/NEJM_2023_submission/'];
stim_cfg.proc_subdir              = 'same_dutyCycle_and_freq';

%**note** plots are ranked by mean of FIRST pain metric specified below
stim_cfg.plt_metrics              = {'mayoNRS', 'painVAS', 'MPQtotal', 'unpleasantVAS'};

stim_cfg.exclude_VAS_50s.decision = true;
stim_cfg.exclude_VAS_50s.plot     = false;

stim_cfg.seperate_by              = {'percentDutyCycle', 'rateInHz'};
stim_cfg.include_cl               = true;

stim_cfg.pt_lbls.RCS02            = 'Patient 1';
stim_cfg.pt_lbls.RCS05            = 'Patient 2';

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
        = plot_stim_groups(stim_cfg, pts{i}, stimGroups_parsi);
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
stim_cfg                   = [];
stim_cfg.proc_dir          = [dirs.rcs_pia, 'beh_analysis/pain_per_DBS_parameters/NEJM_2023_submission/'];
stim_cfg.proc_subdir       = 'Stage 3 (ol-DBS_versus_sham)';


stim_cfg.plt_metrics       = {'mayoNRS','painVAS',  'MPQtotal', 'unpleasantVAS'};
stim_cfg.pt_id             = 'RCS02';

stim_cfg.pt_meta           = pt_META;
stim_cfg.N_days            = duration('3:00:00:00');

stim_cfg.seperate_by       = {'R_stim_groups'};

close all
set(0,'DefaultFigureVisible','off')
s3.RCS02 = struct;

[s3.RCS02.groups, s3.RCS02.sham_vs_stim_stats]...
    ...
    = plot_RCS02_stage3_v2(...
    ...
stim_cfg, stimGroups.RCS02.s3{1}, stimLog.RCS02R);
%% seperate out RCS05's blinded testing

stim_cfg                   = [];
stim_cfg.proc_dir          = [dirs.rcs_pia, 'beh_analysis/pain_per_DBS_parameters/NEJM_2023_submission/'];
stim_cfg.proc_subdir       = '/blinded testing/';


stim_cfg.plt_metrics       = {'mayoNRS','painVAS',  'MPQtotal', 'unpleasantVAS'};
pt_id                 = 'RCS05';

stim_cfg.pt_meta           = pt_META;


s2_blind.RCS05        = plot_RCS05_blinded(stim_cfg, pt_id, stimGroups, INS_logs_API_t_synced);

%% RCS04 --> stim groups after starting buprenorphine (after July 2022 home visit)
stim_cfg                   = [];
stim_cfg.pt_id             = 'RCS04_July22_Apr23';


i_epoch          = find(ge(REDcap.RCS04.time, pt_META.RCS04.dates(end) + duration('24:00:00')));

[~, stimGroups.(stim_cfg.pt_id)] ...
    ...
    = make_stim_groups(...
    ...
 'RCS04', REDcap.RCS04L(i_epoch : end,:), REDcap.RCS04R(i_epoch: end, :), pt_META.RCS04);

 stimGroups.(stim_cfg.pt_id).('Pre-trial baseline') = {fluct.('RCS04')};
%%
%cfg                          = [];
stim_cfg.min_n_reports            = 3;

stim_cfg.min_n_reports_subspace   = 2;

stim_cfg.proc_dir          = [dirs.rcs_pia, '/beh_analysis/pain_per_DBS_parameters/'];
stim_cfg.proc_subdir       = 'same_dutyCycle_and_freq';

stim_cfg.plt_metrics       = {'painVAS', 'mayoNRS',  'MPQtotal', 'unpleasantVAS'};

stim_cfg.seperate_by       = {'percentDutyCycle', 'rateInHz'};

stim_cfg.include_cl        = false;

stim_cfg.exclude_VAS_50s.decision   = false;
stim_cfg.exclude_VAS_50s.plot       = true;

stim_cfg.pt_lbls.RCS04_July22_Apr23 = 'RCS04_July22_Apr23';

    plted_stim_groups.(stim_cfg.pt_id )= plot_stim_groups(stim_cfg,'RCS04_July22_Apr23', stimGroups);




%
%
%
%% distributions of pain metrics and relationship BETWEEN pain metrics
% pts = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};
% 
% cfg                     = [];
% cfg.dates               = 'AllTime';
% 
% for i = 4%1:length(pts)
%     cfg.pt_id  = pts{i};       
% 
%     %plot_hist(cfg, REDcap);          
%     plot_versus(cfg, REDcap);
% end

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
% cfg is the configuration for plotting and clustering w/n 'plot_psyphy_space' fxn
cfg                 = [];

% specify 'manual' for GUI to appear and select clustering
% OR
% specify 'top_2' clusters to return 2 clusters w/ highest density*distance
cfg.CBDP_method     = 'top_2';

% where all raw, z-scored, clustered, fingerprints, and .xlsx of metrics are saved
cfg.proc            = fullfile(dirs.dropbox,'DATA ANALYSIS',...
                             '(Ryan) Stage1_2_3_group_level_behavioral_analysis',...
                             'REDcap_FitBit_only/','psychophysio_fingerprinting/');

cfg.proc_xlsx       = fullfile(rcs_db_cfg.proc_dir, "REDcap/");
% 'true' to cluster based on N principal components needed to describe 95% of variaance
% 'false' to cluster on z-score metrics themselves
cfg.pca             = false;

% name of subdirectoy and pt ids (used folder creation/saving)
cfg.proc_subdir    = 'psy_only';
cfg.save_fig       = true;

% no need to print figures if they're all saved out (can delete these for
% troubleshooting purposes)
close all
set(0,'DefaultFigureVisible','off')

for i =  1  :  length(pts)

    [r_cap] = preproc_prisim_rcap(pts{i}, cfg, REDcap.(pts{i}));

    cfg.var_oi             = r_cap.vars_oi;
    cfg.view_scatter3      = {'mayoNRS', 'painVAS', 'unpleasantVAS'};
    cfg.color_var          = 'MPQtotal';

    switch pts{i}
        case 'RCS05'
            cfg.CBDP_method     = 'top_3';
            
        otherwise
            cfg.CBDP_method     = 'top_2';
    end

    plot_psyphy_space(pts{i}, cfg, r_cap.tbl)
end

%%
%
%% ******************** in-progress versions below ******************** %%
%
%%
% last N days for: 
set(0,'DefaultFigureVisible','on')
cfg_rcap                     = [];

cfg_rcap.pt_id               = 'RCS04';
cfg_rcap.dates               = 'PreviousDays';
cfg_rcap.ndays               = 7;

cfg_rcap.stage_dates         = stage_dates;

cfg_rcap.subplot             = true;
cfg_rcap.sum_stat_txt        = true;
cfg_rcap.stim_parameter      = '';

    plot_timeline(cfg_rcap, REDcap, PFS_sum_stats);

%% stage 3 visualization based on dates from Note-taking Google sheet
s3_week = table;
s3_week.starts = {'20-Jul-2023';...
                  '27-Jul-2023';...
                  '03-Aug-2023';
                  '10-Aug-2023';...
                  '17-Aug-2023'...
                  ...
                  ...
                  };

s3_week.ends   = {'27-Jul-2023';...
                  '03-Aug-2023';...
                  '10-Aug-2023';...
                  '17-Aug-2023';...
                  '24-Aug-2023'...
                  ...
                  ...
                  };


set(0,'DefaultFigureVisible','on'); close all

for i_w = 1 : height(s3_week)
    cfg_rcap                     = [];
    
    cfg_rcap.pt_id               = 'RCS02';
    cfg_rcap.dates               = 'DateRange';
    cfg_rcap.date_range          = [s3_week.starts(i_w), s3_week.ends(i_w)];
    
    cfg_rcap.stage_dates         = stage_dates;
    
    cfg_rcap.subplot             = true;
    cfg_rcap.sum_stat_txt        = true;
    cfg_rcap.stim_parameter      = '';
    
    plot_timeline_s3(cfg_rcap, REDcap, PFS_sum_stats);

end
%% 
cfg_rcap                     = [];

cfg_rcap.pt_id               = 'RCS02';
cfg_rcap.dates               = 'PreviousDays';
cfg_rcap.ndays               = 6;

cfg_rcap.stage_dates         = stage_dates;

cfg_rcap.subplot             = true;
cfg_rcap.sum_stat_txt        = true;
cfg_rcap.stim_parameter      = '';

    plot_timeline_s3(cfg_rcap, REDcap, PFS_sum_stats);








%%
cfg_rcap                     = [];

cfg_rcap.pt_id               = 'RCS02';
cfg_rcap.dates               = 'PreviousDays';
cfg_rcap.ndays               = 4;

cfg_rcap.stage_dates         = stage_dates;

cfg_rcap.subplot             = true;
cfg_rcap.sum_stat_txt        = true;
cfg_rcap.stim_parameter      = '';

    plot_timeline(cfg_rcap, REDcap, PFS_sum_stats);
    