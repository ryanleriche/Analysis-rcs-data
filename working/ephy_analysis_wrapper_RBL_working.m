% see CONFIG_ephy_analysis.m script to specify directories, pull REDcap, 
% and generate subdirectories
CONFIG_ephy_analysis;

%% import REDcap daily, weekly, and monthly surveys from stages 1,2 and 3
% as of Apr. 2023, only daily surveys are analysis-ready/organized
REDcap                = RCS_redcap_painscores(rcs_API_token);

%% import RCS databases, and INS logs per pt side
cfg_rcs                    = [];
cfg_rcs.load_EventLog      = true;

% option to load previous database for efficient processing
cfg_rcs.ignoreold_db                = false;
cfg_rcs.ignoreold_INS_logs          = false;
cfg_rcs.ignoreold_par_db            = false;
cfg_rcs.ignoreold_td_parsing        = false;

cfg_rcs.raw_dir                     = [dirs.rcs_pia, 'raw/'];
cfg_rcs.proc_dir                    = [dirs.rcs_pia, 'processed/'];
cfg_rcs.anal_dir                    = [dirs.rcs_pia, 'ephy_analysis/'];

cfg_rcs.pp_RCS_TD_subset            = 'stage1_only';
cfg_rcs.pp_fft                      = '30s_pre_survey';


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
        cfg_rcs , pt_sides{i});

    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings
    [par_db.(pt_sides{i}), ~]...
        ...
        = makeParsedDatabaseRCS(...
        ...
        cfg_rcs , pt_sides{i}, db);

    %%% find stim settings during each REDcap survey
    [stimLog.(pt_sides{i}), REDcap.(pt_sides{i})]  ...
        ...
        = align_REDcap_to_stimLog(...
        ...
        cfg_rcs , pt_sides{i}, db, REDcap);

    %%% plot API to INS latency
    % (i.e., how long is the INS ahead OR behind internet time [generally behind])
        plt_INS_lat_per_session(cfg_rcs , pt_sides{i}, db);
    
    %%% plot impedance over time (contacts to case during a "Lead Integrity Test")
        plt_impedance_per_session(cfg_rcs , pt_sides{i}, db);


end
%% load RCS TD and format as FFT
for i = 1 : length(pt_sides)
  %%% save RCS time-domain data wrt REDcap survey
   ft_form_TD.(pt_sides{i}) ...
       = rcs_TD_peri_survey(cfg_rcs, dirs.rcs_preproc, pt_sides{i}, par_db, REDcap);

  %%% save FFT per channel across sessions
    % summary spectrograms over sessions (i.e., days to months worth of
    % spectra saved in cfg_rcs.proc_dir subfolders)
        fft_taper_comparison(cfg_rcs, dirs, pt_sides{i}, ft_form_TD);

  %%% pre-process FieldTrip FFT spectra
    % (1) remove noisy sessions by visual inspection
    % (2) interpolate over RC+S spectra artifacts in frequency-domain
        pp_rcs_spectra(cfg_rcs, dirs, pt_sides{i})

end

for i = 1 : length(pt_sides)

    load_dir = fullfile(dirs.rcs_preproc,'spectra_per_sess', pt_sides{i}, [cfg_rcs.pp_RCS_TD_subset, ' (pp work up)'], cfg_rcs.pp_fft);
    
    load(fullfile(load_dir, [pt_sides{i} '_ft_form_pp_FFT_struct.mat'])); %#ok<*LOAD>    
    
    
    writetable(ft_freq_pp.rcs.par_db, fullfile(load_dir, [pt_sides{i}, '_par_db.xlsx']));

end
%% filter out sessions with nonsequitor sense settings
%{
after running 'per_rcs_session_fft" on Chang Lab server, noticed not all
stage 1 streaming sessions had the same sense defintion.

* return N of RCS sessions TD settings per channel
    * have indices from bulk-ran Stage 0

* 

%} 
close all
pt_sides               = {'RCS02R','RCS04R','RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};

set(0,'DefaultFigureVisible','on')
for i = 1 : length(pt_sides)
    plt_rcs_fft_w_same_TD_settings(cfg_rcs , pt_sides{i}, dirs)
end



%% merge RCS and NK spectra 
%{
From FOOOFed RCS stage 1 spectra, identfy band peaks + width while plotting
all sensing hemisphere is single figure.

%}

pt_sides    = {'RCS02R','RCS04R','RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};
pts         = {'RCS02', 'RCS04', 'RCS05', 'RCS06','RCS07'};          

[rcs.sense_chan, rcs.fft_bins_inHz]   = import_fooof(...
                                            fullfile(dirs.rcs_preproc, 'spectra_per_sess/'),...
                                             pt_sides, ...
                                            '/stage1_only (pp work up)/30s_pre_survey/');

[nk.sense_chan, nk.fft_bins_inHz]     = import_fooof(...
                                            fullfile(dirs.nk_preproc, 'fooof_specs (30s_prior_to_survey, updated channels)'),...
                                            pts,...
                                            '');
%% compare peak frequencies between NK and RC+S spectra
%{

% find distance between peaks




%}

plt_nk_rcs_spectra(dirs, pt_sides, nk, rcs, 'ALL_hemispheres_nk_rcs (30s pre survey)')


%%
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

