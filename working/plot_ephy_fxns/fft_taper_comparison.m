function fft_taper_comparison(cfg, dirs, pt_side_id, ft_form_TD)

ft_TD_struct = ft_form_TD.(pt_side_id);

%%
N_samp = seconds(ft_TD_struct.rcs.td_cfg.epoch_dur) * ft_TD_struct.fsample;

for j = 1 :length(ft_TD_struct.time)

   ft_TD_struct.time{j}  = ft_TD_struct.time{j}(end-(N_samp-1):end);
   ft_TD_struct.trial{j} = ft_TD_struct.trial{j}(:, end-(N_samp-1):end);

end


%%
K_tapers   = 3;
time_width = 2;

tapsmofrq = (K_tapers +1) / (time_width *seconds(ft_TD_struct.rcs.td_cfg.epoch_dur));

cfg_freq            = [];
cfg_freq.output     = 'pow';       % return the power-spectra
cfg_freq.method     = 'mtmfft';    % analyses an entire spectrum for the entire data length, implements multitaper frequency transformation.
cfg_freq.taper      = 'dpss';      % multiple tapers based on discrete prolate spheroidal sequences (DPSS), also known as the Slepian sequence.

cfg_freq.pad        ='nextpow2';
cfg_freq.foi        = 1:.25:100;    % frequencies of interest

cfg_freq.toi              = 'all';
cfg_freq.t_ftimwin        = time_width;
cfg_freq.tapsmofrq   = tapsmofrq;  % number, the amount of spectral smoothing through multi-tapering. 
%                                  % Note that 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
cfg_freq.keeptrials  = 'yes';
cfg_freq.polyremoval = -1;

    [ft_freq_mtm] = ft_freqanalysis(cfg_freq, ft_TD_struct);

ft_freq_mtm.rcs = ft_TD_struct.rcs;

%%

cfg_freq            = [];
cfg_freq.output     = 'pow';       
cfg_freq.method     = 'mtmfft';    
cfg_freq.taper      = 'hanning';

cfg_freq.pad        ='nextpow2';
cfg_freq.foi        = 1:.25:100;    % frequencies of interest

cfg_freq.toi              = 'all';

cfg_freq.keeptrials  = 'yes';
cfg_freq.polyremoval = -1;

    [ft_freq_han] = ft_freqanalysis(cfg_freq, ft_TD_struct);

ft_freq_han.rcs = ft_TD_struct.rcs;

%%
fig_dir = fullfile(cfg.anal_dir,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset,cfg.pp_fft,'/');

cfg.z_score         = false;
cfg.freq            = [0.5, 80];
cfg.pwr_limits      = [-10, -1];
cfg.title           = 'taper_comparison';

cfg.txt_leg         = ["FieldTrip's MTM", "FieldTrip's Hanning"];

    plt_rcs_PSDs(cfg, pt_side_id, fig_dir, ft_freq_mtm, ft_freq_han)

%%
save_dir    = fullfile(dirs.rcs_preproc,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset,cfg.pp_fft,'/');
ft_fft_file = [save_dir, pt_side_id, '_ft_form_FFT_struct.mat'];

    save(ft_fft_file, 'ft_freq_mtm', '-v7.3');

end
