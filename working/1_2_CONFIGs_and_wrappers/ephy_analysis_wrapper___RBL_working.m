% see CONFIG_ephy_analysis.m script to specify directories, pull REDcap, 
% and generate subdirectories
ephy_analysis_CONFIG;

%% import REDcap daily, weekly, and monthly surveys from stages 1,2 and 3
% as of Apr. 2023, only daily surveys are analysis-ready/organized
REDcap                = RCS_redcap_painscores(rcs_API_token);

%% import RCS databases, and INS logs per pt side
% option to load previous database for efficient processing
sub_cfg.ignoreold_db                = false;
sub_cfg.ignoreold_INS_logs          = false;
sub_cfg.ignoreold_par_db            = true;
sub_cfg.ignoreold_td_parsing        = false;


sub_cfg.pp_RCS_TD_subset            = 'stage1_only'; 
%sub_cfg.pp_RCS_TD_subset           = 'network_mapping';
sub_cfg.pp_fft                      = '30s_pre_survey'; 

% specify patient hemispheres
%%% pts to update database from scratch locally:
pt_sides               = {'RCS02R','RCS04R','RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};

%% main loop (per pt hemisphere)
for i = 1 : length(pt_sides)
    %%% process RCS .jsons into searchable database 
    [db.(pt_sides{i}), bs.(pt_sides{i})] ...
        ...
        = makeDatabaseRCS_Ryan(...
        ...
        sub_cfg , pt_sides{i});

    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings
    [par_db.(pt_sides{i}), ~]...
        ...
        = makeParsedDatabaseRCS(...
        ...
        sub_cfg , pt_sides{i}, db);

    %%% find stim settings during each REDcap survey
    [stimLog.(pt_sides{i}), REDcap.(pt_sides{i})]  ...
        ...
        = align_REDcap_to_stimLog(...
        ...
        sub_cfg , pt_sides{i}, db, REDcap);

    %%% plot API to INS latency
    % (i.e., how long is the INS ahead OR behind internet time [generally behind])
%         plt_INS_lat_per_session(sub_cfg , pt_sides{i}, db);
%     
%     %%% plot impedance over time (contacts to case during a "Lead Integrity Test")
%         plt_impedance_per_session(sub_cfg , pt_sides{i}, db);


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
        sub_cfg, pt_sides{i}, dirs.rcs_preproc,...
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
       = rcs_TD_peri_survey_2(sub_cfg, dirs.rcs_preproc, pt_sides{i});
end


%% from FieldTrip, generate power spectra peri-survey
%%%% save Hanning window and multitaper multiplication (MTM) for comparison
% 07/17/23 --> using MTM
    % visually saw less variance in spectra, but
    % decision seems arbitrary, some cases where MTM "spreads" RC+S artifact
    % noise, but >40Hz is not analyzed anyways by FOOOF
%%%%%%

% first w/ 30s pre-survey
sub_cfg.min_pre_rcap   = duration('0:00:30');
sub_cfg.cont_chunk_dur = duration('0:00:10');
sub_cfg.pp_fft          = '30s_pre_survey';


pt_sides               = {'RCS02R','RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};


for i = 1 : length(pt_sides)
    %%% save FFT per channel across sessions
    % summary spectrograms over sessions (i.e., days to months worth of
    % spectra saved in cfg_rcs.proc_dir subfolders)
        fft_taper_comparison_2(sub_cfg, dirs, pt_sides{i}, ft_form_TD);
end

for i = 1 : length(pt_sides)
    % (1) remove noisy sessions by visual inspection
    % (2) interpolate over RC+S spectra artifacts in frequency-domain
        pp_rcs_spectra_2(sub_cfg, dirs, pt_sides{i})
end

%%% then w/ 5 min pre-survey
%{
sub_cfg.min_pre_rcap   = duration('0:05:00');
sub_cfg.cont_chunk_dur = duration('0:00:10');
sub_cfg.pp_fft          = '5m_pre_survey';

for i = 1 : length(pt_sides)

    %%% save FFT per channel across sessions
    % summary spectrograms over sessions (i.e., days to months worth of
    % spectra saved in cfg_rcs.proc_dir subfolders)
    fft_taper_comparison_2(sub_cfg, dirs, pt_sides{i}, ft_form_TD);

        % (1) remove noisy sessions by visual inspection
    % (2) interpolate over RC+S spectra artifacts in frequency-domain
        pp_rcs_spectra_2(sub_cfg, dirs, pt_sides{i})

end
%}
%% merge RCS and NK spectra 
%{
From FOOOFed RCS stage 1 spectra, identfy band peaks + width while plotting
all sensing hemisphere is single figure.

%}

%pt_sides    = {'RCS02R','RCS04R','RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};
close all
pt_sides    = {'RCS02R', 'RCS04L', 'RCS04R','RCS05L', 'RCS05R', 'RCS06L','RCS06R','RCS07L', 'RCS07R'};
pts         = {'RCS02', 'RCS04', 'RCS05', 'RCS06','RCS07'};          


[rcs_30s.sense_chan, rcs_30s.fft_bins_inHz]   = import_fooof(...
                                            fullfile(dirs.rcs_preproc, 'spectra_per_sess/'),...
                                             pt_sides, ...
                                            '/stage1_only (pp work up)/30s_pre_survey/');

[nk.sense_chan, nk.fft_bins_inHz]     = import_fooof(...
                                            fullfile(dirs.nk_preproc, 'fooof_specs (30s_prior_to_survey, updated channels)'),...
                                            pts,...
                                            '');
%%
[nk, rcs_30] = plt_nk_rcs_spectra(dirs, pt_sides, ...
                                  nk, rcs_30s, ...
                                  'ALL_hemispheres_nk_rcs (fixed, 30s_pre_survey)');
%%
peak_sum_tbl = find_nearest_FOOOFd_peak(nk, rcs_30s);

for i = 1 : length(pt_sides)

    rcs_chs = rcs_30.sense_chan.Properties.RowNames;
    i_rcs   = find(contains(rcs_chs, pt_sides{i}));
    
    nk_chs  = nk.sense_chan.Properties.RowNames;
    i_nk    = find(contains(nk_chs, pt_sides{i}));

    peak_sum_tbl.abs_spec_diff_in_db(i) = mean(...
                        abs(...
                            nk.sense_chan{i_rcs,"mean_periodic_spec"}{1} ...
                            - ...
                            rcs_30.sense_chan{i_nk,"mean_periodic_spec"}{1}...
                            )...
                        );
end

peak_sum_tbl.Properties.RowNames = pt_sides;

peak_sum_tbl = sortrows(peak_sum_tbl, 'Row');

% plot_primary_peaks(dirs, peak_freq_tbl, nk, rcs_30s, ...
%     'ALL_hemispheres_nk_rcs (fixed, 30s_pre_survey)---single band comparison')

plot_primary_peaks_2(dirs, peak_sum_tbl, nk, rcs_30s, ...
    'ALL_hemispheres_nk_rcs (fixed, 30s_pre_survey)---major peaks')

%%
s0_s1_path = fullfile(dirs.dropbox,...
        '/MANUSCRIPTS/2023_Staged_Neuropsych_DBS_Feasibility/',...
        'staged_neuropsych_DBS_feasibility_data/',...
         'Prasad_electrode_locations_stage0_and_stage1_ver4.mat');


tmp            = load(s0_s1_path);
s0_s1_dist_tbl = tmp.data;              clear tmp;


s0_s1_dist_tbl(...
                cellfun(@isempty, s0_s1_dist_tbl.SlicerClosestLabel),...
                :) = [];


s0_s1_dist_tbl.Subject = upper(s0_s1_dist_tbl.Subject);
%%%
hemi_lbl_mat = split(...
                    replace(...   
                        replace(...
                                peak_sum_tbl.labels_1, {' ', '-'}, '_'), ...
                        {'RACCa', 'CM'}, {'RACC', 'CMPF'}),...
                    '_');

first_cont = cellfun(@(x, y) sprintf('%s%s', x, y), ...
            hemi_lbl_mat(:,2), hemi_lbl_mat(:,3), ...
            'UniformOutput', false);

second_cont = cellfun(@(x, y) sprintf('%s%s', x, y), ...
            hemi_lbl_mat(:,2), hemi_lbl_mat(:,4), ...
            'UniformOutput', false);

peak_sum_tbl.labels_1_1 = first_cont;
peak_sum_tbl.labels_1_2 = second_cont;


RCS_lbl_rmv = {'RIFG', 'RSFG','RSGC', 'LSGC','LCaud','+'};
TW_lbl_add  = {'RMFG', 'RMFG', 'RACC', 'LACC','LCN',''};
hemi_lbl_mat = split(...
                    replace(...   
                        replace(...
                                peak_sum_tbl.labels_2, {' ', '-'}, '_'), ...
                                RCS_lbl_rmv, TW_lbl_add),...
                    '_');

% translate from 0:3, 8:11 (i.e., RC+S labels) to Tom Wozny's 1:4 labels
hemi_lbl_mat(:, 3:4) = replace(hemi_lbl_mat(:, 3:4), ...
                               compose('%g', [8:11, 0:3]),...
                               compose('%g', [1:4, 1:4]));

first_cont = cellfun(@(x, y) sprintf('%s%s', x, y), ...
            hemi_lbl_mat(:,2), hemi_lbl_mat(:,3), ...
            'UniformOutput', false);

second_cont = cellfun(@(x, y) sprintf('%s%s', x, y), ...
            hemi_lbl_mat(:,2), hemi_lbl_mat(:,4), ...
            'UniformOutput', false);

peak_sum_tbl.labels_2_1 = first_cont;
peak_sum_tbl.labels_2_2 = second_cont;


%%

for i = 1 : length(pt_sides)

    pt_oi     = peak_sum_tbl.labels_1{i}(1:5);
    i_pt_oi   = strcmpi(s0_s1_dist_tbl.Subject, pt_oi);


    hemi_dist = s0_s1_dist_tbl(i_pt_oi, :);

    i_label   = contains(hemi_dist.Label,...
                    {peak_sum_tbl.labels_2_1{i}, peak_sum_tbl.labels_2_2{i}});

    peak_sum_tbl.avg_SliderEuclidainMin(i) ...
        ...
        = mean(hemi_dist(i_label, :).SlicerEuclidianMin);

    peak_sum_tbl.hemi_dist_tbl{i} = hemi_dist(i_label, :);


    peak_sum_tbl.bf01_max(i)...
        ...
        =max(peak_sum_tbl.bf01{i});


    peak_sum_tbl.bf01_mean(i)...
        ...
        =mean(peak_sum_tbl.bf01{i});

end

% close all
% set(0,'DefaultFigureVisible','on')

figure('Units','inches','Position',[5, 5, 5, 5]);

colors = brewermap(length(pt_sides), 'Set1');

for i = 1: length(pt_sides)
scatter(peak_sum_tbl.avg_SliderEuclidainMin(i), ...
        peak_sum_tbl.abs_spec_diff_in_db(i), ...
        200, colors(i, :),'filled');
hold on

end

legend(pt_sides, 'Location','northoutside', 'NumColumns',3);
ylabel('absolute spectra difference (db)')
xlabel('SliderEuclidainMin between bipolar pair')

set(gca, 'FontSize', 14)

i_keep = ~strcmp(peak_sum_tbl.Row, 'RCS05R');

[rho, p] = corr(peak_sum_tbl.abs_spec_diff_in_db(i_keep),...
                 peak_sum_tbl.avg_SliderEuclidainMin(i_keep), 'type','Pearson');


title(sprintf("Pearson's correlation: %.2f \rho; %.1e \p-value", rho, p))


%%













%%
pt_sides    = {'RCS02R', 'RCS04R', 'RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};
pts         = {'RCS02', 'RCS04', 'RCS05', 'RCS06','RCS07'};          


[rcs_5m.sense_chan, rcs_5m.fft_bins_inHz]   = import_fooof(...
                                            fullfile(dirs.rcs_preproc, 'spectra_per_sess/'),...
                                             pt_sides, ...
                                            '/stage1_only (pp work up)/5m_pre_survey/');

plt_nk_rcs_spectra(dirs, pt_sides, nk, rcs_5m, 'ALL_hemispheres_nk_rcs (fixed, 5m_pre_survey)')
%%
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
%
%
%
%
%
%

hemi_chan_lbls = hemi_sense_chan.Properties.RowNames;



i_cnt = 1;

%i_sense = 1:height(hemi_sense_chan);

i_sense = [1, 2, 4, 7:2:13, 14, 17];
for j  = i_sense
    nexttile(i_cnt)

    plt_spectrum = hemi_sense_chan{j,"pwr_spectra_aperiodic_rmv"}{1};

    [lineout, ~] = stdshade(plt_spectrum, 0.2, 'k', fft_bins_inHz);

    input = {lineout.YData, fft_bins_inHz, 'MinPeakProminence',0.1,'Annotate','extents'};

    [a, b, c, d] ...
     ...
        =  findpeaks(input{:});

     hemi_sense_chan{j, 'ss_avg_pks_pwr'} = {a};
     hemi_sense_chan{j, 'ss_avg_pks_freq'} = {b};

     hemi_sense_chan{j, 'ss_avg_half_prom_freq_width'} = {c};
     hemi_sense_chan{j, 'ss_avg_hprom_pwr'}            = {d};

    
    %%% nice built-in plotting call
    findpeaks(input{:}); hold on;

    ax = gca;
    h = findobj(ax, 'tag', 'HalfProminenceWidth');  h.LineWidth = 2;  h.Color = 'g';
    h = findobj(ax, 'tag', 'Prominence');           h.LineWidth = 2;
    
    h = findobj(ax, 'tag', 'Peak');                 h.MarkerSize = 10; 

  
    %ylim([-0.1, 1])
    ylabel('log_{10}(mV^{2}/Hz)');      xlabel('Frequency (Hz)');
    title(hemi_chan_lbls(j), 'Interpreter','none');

 
    [lineout, ~] = stdshade(plt_spectrum, 0.2, 'k', fft_bins_inHz);
    lineout.LineStyle = '-';              lineout.LineWidth = 2.25;
    ylim([-0.1, .8])

    i_cnt = i_cnt +1;
    if j ==1

        ax = gca;

        ax.Legend.Position    = [0.55 0.8 0.4 0.125];
        ax.Legend.FontSize     = 12;

        ax.Legend.String(end)  = {'Mean periodic spectra'};


        nexttile(i_cnt)
        delete(gca)
        i_cnt = i_cnt +1;

    else
        legend off

    end
end

exportgraphics(gcf, [dirs.rcs_pia,'ephy_analysis/spectra_per_sess/group_rcs_s1_spectra_peaks.png'])

%%
i_sense = [1, 3 : 2 : height(hemi_sense_chan)];

%%



%% w/ power spectrum, aperiodic fit, and fooofed spectrum
% explore high/low pain decoding for RCS02R
cfg                 = [];


cfg.raw_dir        = '/Volumes/DBS Pain 3/rcs_device_data/';
cfg.raw_subdir     = '/fft_per_sess/Stage1_only/';


pt_side_id         = pt_sides{i};


%%%
dirs.rcs_preproc= [cfg.raw_dir, pt_side_id(1:end-1),...
            cfg.raw_subdir, pt_side_id, '/'];

load(...
      sprintf('%s/%s_pwrspectra_by_sess.mat', dirs.rcs_preproc, pt_side_id));

load(...
      sprintf('%s/%s_fft_bins_inHz.mat', dirs.rcs_preproc, pt_side_id));

parsed_db_oi   = readtable(...
                   sprintf('%s/%s_parsed_db_oi.xlsx', dirs.rcs_preproc, pt_side_id));

ch_names      = load(...
                    sprintf('%s/%s_ch_names.mat', dirs.rcs_preproc, pt_side_id));
ch_names = ch_names.ch_names;
%%%

for i_ch = 1 : length(ch_names)

    if i_ch == 1
        fooof_freqs ...
            = readmatrix(...
                sprintf('%sexported_fooof_data/fooof_freqs.xlsx', dirs.rcs_preproc));

        fooof_freqs = fooof_freqs(2:end);
    end

    tmp_spectra ...
        = readmatrix(...
        sprintf('%sexported_fooof_data/%s_pwr_spectra_aperiodic_rmv.xlsx', dirs.rcs_preproc, ch_names{i_ch}));

    if i_ch == 1

       pwr_spectra_aperiodic_rmv = ...
           nan([length(ch_names), ...
           size(tmp_spectra(2:end, :))...
           ]);

    end
    pwr_spectra_aperiodic_rmv(i_ch, :, :) = 10.^tmp_spectra(2:end, :);

end

cfg            = [];
cfg.date_range = sess_oi_date_range.(pt_sides{i}(1:end-1));

cfg.z_score    = false;
cfg.freq       = [2, 80];
cfg.pwr_limits = [-.5, 1];

close all
plt_session_spectra(cfg, fooof_freqs, pt_side_id, pwr_spectra_aperiodic_rmv, ...
                                    ch_names, parsed_db_oi, [dirs.rcs_preproc, 'fooofed_spectra/'])


cfg.z_score    = true;
cfg.freq       = [2, 80];
cfg.pwr_limits = [-3, 3];

plt_session_spectra(cfg, fooof_freqs, pt_side_id, pwr_spectra_aperiodic_rmv, ...
                                    ch_names, parsed_db_oi, [dirs.rcs_preproc, 'fooofed_spectra/'])



%%% now see meaningful cut-offs for high versus low pain binarization

redcap  = REDcap.(pt_side_id(1:end-1));

i_dates_oi  = find(...
                    ge(redcap.time, cfg.date_range{1}) &...
                    le(redcap.time, cfg.date_range{2}));

redcap  = redcap(i_dates_oi, :);


figure('Units', 'Inches', 'Position', [0, 0, 14, 10])
histogram(redcap.MPQtotal);      %xlim([0, 45])
ylabel('MPQtotal')


xline(prctile(redcap.MPQtotal, 10:20:90), 'k', compose('%gth percentile',10:20:90))


high_pain_cutoff = prctile(redcap.MPQtotal, 66);

low_pain_cutoff = prctile(redcap.MPQtotal, 33);
% 
% i_low   = redcap.MPQtotal <= low_pain_cutoff;
% i_high  = redcap.MPQtotal >= high_pain_cutoff;
%%%
parsed_db_oi.timeStart.TimeZone= 'America/Los_Angeles';
parsed_db_oi.timeStop.TimeZone = 'America/Los_Angeles';

for i_sess = 1:height(parsed_db_oi)

    i_rcap    = find(ge(redcap.time, parsed_db_oi.timeStart(i_sess) - duration('0:10:00')) &...
                 le(redcap.time, parsed_db_oi.timeStop(i_sess)));

    parsed_db_oi.ind_rcap{i_sess} = i_rcap;

    if ~isempty(i_rcap)
        parsed_db_oi.MPQtotal_mean(i_sess) = mean(redcap.MPQtotal(i_rcap));
        parsed_db_oi.MPQtotal_std(i_sess) = std(redcap.MPQtotal(i_rcap));

    else
        parsed_db_oi.MPQtotal_mean(i_sess) = NaN;
        parsed_db_oi.MPQtotal_std(i_sess) = NaN;

    end
end


i_low  = parsed_db_oi.MPQtotal_mean <= low_pain_cutoff;
i_high = parsed_db_oi.MPQtotal_mean >= high_pain_cutoff;

parsed_db_oi.pain_state = repmat({''}, height(parsed_db_oi), 1);
parsed_db_oi.pain_state(i_high) = {'high_pain'};

parsed_db_oi.pain_state(i_low) = {'low_pain'};

%% work-up to high/low pain spectra across channels

cfg.freq       = [2, 80];
cfg.pwr_limits = [-.5, 1];

pwr_spectra = pwrspectra_by_sess;
%pwr_spectra = pwr_spectra_aperiodic_rmv;

freqs   = fft_bins_inHz;
%freqs   = fooof_freqs;

figg = figure;
figg.Units    = 'Inches';
figg.Position= [0, 3, 18, 12];

i_freq_oi     = find((freqs >= cfg.freq(1) & freqs <= cfg.freq(2)));


tiledlayout(2,2, 'Padding','tight')

sgtitle(sprintf('%s | %s–%s, Stage 1 (prior to analgesic stim testing)\n%g sessions total',...
    pt_side_id, cfg.date_range{1}, cfg.date_range{2}, height(parsed_db_oi)),...
    'FontSize', 24)


i_high = strcmp(parsed_db_oi.pain_state, 'high_pain');
i_low  = strcmp(parsed_db_oi.pain_state, 'low_pain');

for i_ch = 1:4
    nexttile; hold on

    title(sprintf('Ch%g | contacts %s', i_ch - 1, ch_names{i_ch}));

    plt_spectrum =  log10(squeeze(pwr_spectra(i_ch, :, i_freq_oi)));

    [lineout, ~] = stdshade( plt_spectrum(i_low, :), 0.2, 'b', freqs(i_freq_oi));
        lineout.LineStyle = '-';              lineout.LineWidth = 4;

    [lineout, ~] = stdshade( plt_spectrum(i_high, :), 0.2, 'r', freqs(i_freq_oi));
        lineout.LineStyle = '-';              lineout.LineWidth = 4;

%     plot(fooof_freqs(i_freq_oi),...
%         plt_spectrum(i_low, :),  '-b', 'LineWidth', .5); 
% 
%     plot(fooof_freqs(i_freq_oi),...
%         plt_spectrum(i_high, :),  '-r', 'LineWidth', .5, 'Marker','none');

end

