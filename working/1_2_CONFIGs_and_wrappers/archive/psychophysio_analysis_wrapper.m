%% psychophysiological fingerprinting of neuropsychiatric states
%{
Psychological variables:
    assessment of states (pain, depression, anxiety, etc) from
    REDcap surveys

Physiological variables:
    step counts, hours of sleep, mean heart rate, and heart rate
    variability from FitBit

__________________________________________________________________________
High-level processing overview:

    * filter-out fencesitting of VAS 50s--a neccesary alteration to identify natural pain subspaces
    
    * z-score metrics to themselves for easier comparison of magnitude btwn metrics
        -> could consider rescaling by [min, max]
    
    * cluster based on density peaks (CBDP) (Rodriguez & Laio, 2014, Science)
        -> currently, manual selection or "top two" clusters

    * plot mean+-std of metrics w/n clusters (the "fingerprints" per se)

__________________________________________________________________________

Ryan Leriche, May 2023

%}
%%
% see CONFIG_psychophysio_analysis.m script to specify directories, pull REDcap, 
% and generate subdirectories
psychophysio_analysis_CONFIG;

%%
% cfg is the configuration for plotting and clustering w/n 'plot_psyphy_space' fxn
cfg                 = [];

% specify 'manual' for GUI to appear and select clustering
% OR
% specify 'top_two' clusters to return two clusters w/ highest density*distance
cfg.CBDP_method     = 'manual';

% where all raw, z-scored, clustered, fingerprints, and .xlsx of metrics are saved
cfg.proc            = fullfile(dirs.working,'psychophysio_fingerprinting/');

% 'true' to cluster based on N principal components needed to describe 95% of variaance
% 'false' to cluster on z-score metrics themselves
cfg.pca             = false;

% name of subdirectoy and pt ids (used folder creation/saving)
cfg.proc_subdir    = 'psy_only';
rcs_pts            = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};

% no need to print figures if they're all saved out (can delete these for
% troubleshooting purposes)
close all
set(0,'DefaultFigureVisible','off')

for i =  1  :  length(rcs_pts)

    [r_cap] = preproc_prisim_rcap(rcs_pts{i}, cfg, PRISM_r_cap.(rcs_pts{i}));
 
    cfg.var_oi             = r_cap.vars_oi;

    cfg.view_scatter3      = {'mayoNRS', 'painVAS', 'unpleasantVAS'};
    cfg.color_var          = 'MPQtotal';

        plot_psyphy_space(rcs_pts{i}, cfg, r_cap.tbl)

end

%% iterate w/ Presidio pts (DBS for depression trial)
% specify clustering decision method (uses rest of cfg from above)
cfg.CBDP_method     =  'top_two';

% specifying variable names based on visual inspection of files
psy_vars            = {'vas_depression','vas_anxiety', 'vas_energy','hamd2_score','sss'};
phy_vars            = {'TotalMinutesAsleep','TotalTimeInBed', 'RestingHeartRate', 'Steps'};

%%% import Presidio pt data to merged REDcap + FitBit structures
raw_dir    = fullfile(dirs.working, 'presidio_data/behavioral_data');
PR_rcap_fb = import_presidio_csvs(raw_dir);
pt_ids     = fieldnames(PR_rcap_fb);


%%% per file (unique pt and time window for FitBit data)
% one file for now, until standard format is established across Prisim and Presidio pts
i = 3;
pr_beh = PR_rcap_fb.(pt_ids{i});

%%% cluster by psychological reporting of neuropsychiatric states
cfg.proc_subdir    = 'psy_only';
cfg.var_oi         = psy_vars;
cfg.pca            = true;
cfg.view_scatter3  = {'vas_depression','vas_anxiety', 'vas_energy'};
cfg.color_var      =  'hamd2_score';

    plot_psyphy_space(pt_ids{i}, cfg, pr_beh)


%%% cluster by physiological reporting of neuropsychiatric states
cfg.proc_subdir    = 'phy_only';
cfg.var_oi         = phy_vars;
cfg.CBDP_method    = 'manual';

cfg.view_scatter3  =  {'TotalMinutesAsleep','Steps', 'RestingHeartRate'};
cfg.color_var      =  'TotalTimeInBed';

      plot_psyphy_space(pt_ids{i}, cfg, pr_beh)

%%% seeing psychological "energy" clusters--further breakdown by sleeping
%%% physiology
cfg.proc_subdir    = 'psy_phy_arousal';

cfg.var_oi         = {'TotalMinutesAsleep','TotalTimeInBed','vas_energy', 'sss'};

cfg.view_scatter3  =  {'TotalMinutesAsleep','vas_energy','sss'};
cfg.color_var      =  'sss';

      plot_psyphy_space(pt_ids{i}, cfg, pr_beh)





