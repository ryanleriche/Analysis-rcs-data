%% load CONFIG_ephy_wrapper

[dirs,    rcs_API_token,   pcs_API_token, ... -> input and output directories, and API tokens
 PFS,     PFS_sum_stats,...                   -> pain fluctuation study (PFS) data and summary statistics
 pt_META, stage_dates]...                     -> hard-coded pt meta data (RCS Stage dates pulled from patient iternaries, Google Drive folders, etc
...
    = CONFIG_ephy_wrapper;

%% import REDcap daily, weekly, and monthly surveys from stages 1,2 and 3
% as of Apr. 2023, only daily surveys are analysis-ready/organized

REDcap                = RCS_redcap_painscores(rcs_API_token);

%% import RCS databases, and INS logs per pt side
cfg                    = [];
cfg.load_EventLog      = true;

% option to load previous database for efficient processing
cfg.ignoreold_db                = false;
cfg.ignoreold_INS_logs          = false;
cfg.ignoreold_par_db            = true;

cfg.raw_dir                     = [dirs.rcs_pia, 'raw/'];
cfg.proc_dir                    = [dirs.rcs_pia, 'processed/'];

% specify patient hemispheres
%%% pts to update database from scratch locally:
%pt_sides               = {'RCS02R','RCS04R','RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};
pt_sides                = {'RCS06R','RCS06L','RCS07L', 'RCS07R'};

cfg.pp_RCS_TD_subset   = 'stage1_only';

%% main loop (per pt hemisphere)
for i = 1 : length(pt_sides)
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
        = RCS_logs(...
        ...
        cfg, pt_sides{i});

    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings

    [par_db.(pt_sides{i}), ~]...
        ...
        = makeParsedDatabaseRCS(...
        ...
        cfg, pt_sides{i}, db);

    %%% find stim settings during each REDcap survey
    [stimLog.(pt_sides{i}), REDcap.(pt_sides{i})]  ...
        ...
        = align_REDcap_to_stimLog(...
        ...
        cfg, pt_sides{i}, db, REDcap);

    %%% save time-domain LFPs (from RC+S streaming sessions) 
        % across pt hemispheres as clearly labelled .mat
        % first entry allows subsetting of streaming sessions based on
        % criteria within fxn
    par_db_out.(pt_sides{i})  ...
        ...
        = pp_RCS_ss_TD( ...
        ...
        cfg.pp_RCS_TD_subset, pt_sides{i}, dirs.rcs_preproc,...
        par_db, stimLog);

    %%% save FFT per channel across sessions
    % summary spectrograms over sessions (i.e., days to months worth of
    % spectra saved in cfg.proc_dir subfolders)
        per_rcs_session_fft(cfg, dirs.rcs_preproc, pt_sides{i}, par_db_out)

    %%% plot API to INS latency
    % (i.e., how long is the INS ahead OR behind internet time [generally behind])
        plt_INS_lat_per_session(cfg, pt_sides{i}, db);

    %%% plot impedance over time (contacts to case during a "Lead Integrity Test")
        plt_impedance_per_session(cfg, pt_sides{i}, db);
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

pt_sides               = {'RCS06R','RCS06L','RCS07L', 'RCS07R'};


set(0,'DefaultFigureVisible','on')
for i = 1 : length(pt_sides)
    plt_rcs_fft_w_same_TD_settings(cfg, pt_sides{i}, dirs)
end



%%
%{
From FOOOFed RCS stage 1 spectra, identfy band peaks + width while plotting
all sensing hemisphere is single figure.

%}

pt_sides               = {'RCS02R','RCS04R','RCS04L', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};

[hemi_sense_chan, fft_bins_inHz]        = import_fooofed_hemispheres(dirs, pt_sides);

%%
hemi_chan_lbls = hemi_sense_chan.Properties.RowNames;


close all
figure('Units','inches','Position',[0, 3, 8, 10])


tiledlayout(5,2, 'Padding','tight')
sgtitle(sprintf('Power-spectra across RC+S patient sensing hemispheres\nABOVE 1/f (the aperiodic component)'))

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

