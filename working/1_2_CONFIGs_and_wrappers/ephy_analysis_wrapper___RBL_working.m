%% see user_local_CONFIG.m script to specify directories and general set-up
user_local_CONFIG;

% load REDcap pain surveys and set-up directories
general_setup;

% see "dirs" structure for main directories used, and "rcs_db_cfg" for the
% configuration of the RC+S databasing

%% import RCS databases, and INS logs per pt side
% option to load previous database for efficient processing
    rcs_db_cfg.ignoreold_db                = false;
    rcs_db_cfg.ignoreold_INS_logs          = false;
    rcs_db_cfg.ignoreold_par_db            = false;
    rcs_db_cfg.ignoreold_td_parsing        = false;


    %rcs_db_cfg.pp_RCS_TD_subset            = 'stage1_only'; 
    rcs_db_cfg.pp_RCS_TD_subset           = 'network_mapping';
    rcs_db_cfg.pp_fft                      = '30s_pre_survey'; 
    % specify patient hemispheres
    %%% pts to update database from scratch locally:
    %pt_sides               = {'RCS02R','RCS04R','RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};
    
    pt_sides               = {'RCS02R','RCS05R', 'RCS05L'};

%% main loop (per pt hemisphere)
for i = 1 : length(pt_sides)
    %%% process RCS .jsons into searchable database 
    [db.(pt_sides{i}), bs.(pt_sides{i})] ...
        ...
        = makeDatabaseRCS_Ryan(...
        ...
        rcs_db_cfg , pt_sides{i});

    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings
    [par_db.(pt_sides{i}), ~]...
        ...
        = makeParsedDatabaseRCS(...
        ...
        rcs_db_cfg , pt_sides{i}, db);

    %%% find stim settings during each REDcap survey
    [stimLog.(pt_sides{i}), REDcap.(pt_sides{i})]  ...
        ...
        = align_REDcap_to_stimLog(...
        ...
        rcs_db_cfg , pt_sides{i}, db, REDcap);

    %%% plot API to INS latency
    % (i.e., how long is the INS ahead OR behind internet time [generally behind])
%         plt_INS_lat_per_session(rcs_db_cfg , pt_sides{i}, db);
%     
%     %%% plot impedance over time (contacts to case during a "Lead Integrity Test")
%         plt_impedance_per_session(rcs_db_cfg , pt_sides{i}, db);


end
%
%
%
%% pre-process RCS time-domain data
for i = 1 : length(pt_sides)
    %%% save time-domain LFPs (from RC+S streaming sessions) 
        % across pt hemispheres as clearly labelled .mat
        % first entry allows subsetting of streaming sessions based on
        % criteria within fxn
    par_db_out.(pt_sides{i})  ...
        ...
        = pp_RCS_ss_TD( ...
        ...
        rcs_db_cfg, pt_sides{i}, dirs.rcs_preproc,...
        par_db, stimLog);
end
%
%
%
%%% load RCS TD and format as FieldTrip structure
% allows different time windows pre-survey

for i = 1 : length(pt_sides)
    %%% save RCS time-domain data wrt REDcap survey
    ft_form_TD.(pt_sides{i}) ...
       = rcs_TD_peri_survey_2(rcs_db_cfg, dirs.rcs_preproc, pt_sides{i});
end


%% from FieldTrip, generate power spectra peri-survey
%%%% save Hanning window and multitaper multiplication (MTM) for comparison
% 07/17/23 --> using MTM
    % visually saw less variance in spectra, but
    % decision seems arbitrary, some cases where MTM "spreads" RC+S artifact
    % noise, but >40Hz is not analyzed anyways by FOOOF
%%%%%%

% first w/ 30s pre-survey
rcs_db_cfg.min_pre_rcap   = duration('0:00:30');
rcs_db_cfg.cont_chunk_dur = duration('0:00:10');
rcs_db_cfg.pp_fft          = '30s_pre_survey';


pt_sides               = {'RCS02R','RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};


for i = 1 : length(pt_sides)
    %%% save FFT per channel across sessions
    % summary spectrograms over sessions (i.e., days to months worth of
    % spectra saved in cfg_rcs.proc_dir subfolders)
        fft_taper_comparison_2(rcs_db_cfg, dirs, pt_sides{i}, ft_form_TD);
end

for i = 1 : length(pt_sides)
    % (1) remove noisy sessions by visual inspection
    % (2) interpolate over RC+S spectra artifacts in frequency-domain
        pp_rcs_spectra_2(rcs_db_cfg, dirs, pt_sides{i})
end

%% single channel spectrogram, (2023 BRAIN submission)

rcs_db_cfg.pwr_limits      = [-10, -5];
rcs_db_cfg.freq            = [1, 40];


rcs_db_cfg.y_lbl_txt       = 'days since implant';
rcs_db_cfg.rcs_ch_oi       = 'Ch2';
rcs_db_cfg.c_str           = 'Power (db)';



rcs_db_cfg.save_dir       = fullfile(dirs.rcs_pia, 'ephy_analysis/staged_spectra_stability/');

rcs_db_cfg.fig_name        = 'example_s1_spectrum';

plt_daily_spec(rcs_db_cfg, dirs, 'RCS06R');



%%

%%% then w/ 5 min pre-survey
%{
rcs_db_cfg.min_pre_rcap   = duration('0:05:00');
rcs_db_cfg.cont_chunk_dur = duration('0:00:10');
rcs_db_cfg.pp_fft          = '5m_pre_survey';

for i = 1 : length(pt_sides)

    %%% save FFT per channel across sessions
    % summary spectrograms over sessions (i.e., days to months worth of
    % spectra saved in cfg_rcs.proc_dir subfolders)
    fft_taper_comparison_2(rcs_db_cfg, dirs, pt_sides{i}, ft_form_TD);

        % (1) remove noisy sessions by visual inspection
    % (2) interpolate over RC+S spectra artifacts in frequency-domain
        pp_rcs_spectra_2(rcs_db_cfg, dirs, pt_sides{i})

end
%}
%% merge RCS and NK spectra 
%{
From FOOOFed RCS stage 1 spectra, identfy band peaks + width while plotting
all sensing hemisphere is single figure.

%}

%pt_sides    = {'RCS02R','RCS04R','RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};
close all
set(0,'DefaultFigureVisible','on')

pt_sides    = {'RCS02R', 'RCS04L', 'RCS04R','RCS05L', 'RCS05R', 'RCS06L','RCS06R','RCS07L', 'RCS07R'};
pts         = {'RCS02', 'RCS04', 'RCS05', 'RCS06','RCS07'};          


[rcs_30s.sense_chan, ...
    rcs_30s.fft_bins_inHz]  = import_fooof(...
                                            fullfile(dirs.rcs_preproc, 'spectra_per_sess/'),...
                                             pt_sides, ...
                                            '/stage1_only (pp work up)/30s_pre_survey/');

[nk.sense_chan, ...
    nk.fft_bins_inHz]       = import_fooof(...
                                            fullfile(dirs.nk_preproc, 'fooof_specs (30s_prior_to_survey, updated channels)'),...
                                            pts,...
                                            '');
%%
[nk, rcs_30s] = plt_nk_rcs_spectra(dirs, pt_sides, ...
                                  nk, rcs_30s, ...
                                  'ALL_hemispheres_nk_rcs (fixed, 30s_pre_survey)');
%%
peak_sum_tbl = find_nearest_FOOOFd_peak(pt_sides, nk, rcs_30s);

plot_primary_peaks_2(dirs, peak_sum_tbl, nk, rcs_30s, ...
    'ALL_hemispheres_nk_rcs (fixed, 30s_pre_survey)---major peaks')

%%
s0_s1_path = fullfile(dirs.dropbox,...
        '/MANUSCRIPTS/2023_Staged_Neuropsych_DBS_Feasibility/',...
        'staged_neuropsych_DBS_feasibility_data/',...
         'Prasad_electrode_locations_stage0_and_stage1_ver4.mat');


[peak_sum_tbl, s0_s1_dist_tbl] ...
    ...
    = import_interstage_distance(...
    ...
s0_s1_path, peak_sum_tbl);

peak_sum_tbl ...
    ...
    = plot_interstage_distance_to_spectra(...
    ...
dirs, peak_sum_tbl,s0_s1_dist_tbl, 'spec_power_to_dist_error_correlation.png');

%%



%%
pt_sides    = {'RCS02R', 'RCS04R', 'RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};
pts         = {'RCS02', 'RCS04', 'RCS05', 'RCS06','RCS07'};          


[rcs_5m.sense_chan, rcs_5m.fft_bins_inHz]   = import_fooof(...
                                            fullfile(dirs.rcs_preproc, 'spectra_per_sess/'),...
                                             pt_sides, ...
                                            '/stage1_only (pp work up)/5m_pre_survey/');

plt_nk_rcs_spectra(dirs, pt_sides, nk, rcs_5m, 'ALL_hemispheres_nk_rcs (fixed, 5m_pre_survey)')

%% compare peak frequencies between NK and RC+S spectra
%{

visualize NK and RCS peaks and variance

% find distance between peaks


% use NK peaks from mean_fooof_params --> look for nearest RC+S peaks

    * 
    %


%}

% peak_freq_tbl_5m = find_nearest_FOOOFd_peak(nk, rcs_5m);
% 
% plot_primary_peaks(peak_freq_tbl_5m, nk, rcs_5m)
% 
%%