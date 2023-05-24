function fft_taper_comparison(cfg, dirs, pt_side_id, ft_form_TD)

ft_TD_struct= ft_form_TD.(pt_side_id);

%%
K_tapers   = 3;
time_width = ft_TD_struct.sampleinfo(1,2) - ft_TD_struct.sampleinfo(1,1);

tapsmofrq = (K_tapers +1) / (2*time_width);

cfg_freq = [];
cfg_freq.output     = 'pow';       % return the power-spectra
cfg_freq.method     = 'mtmfft';    % analyses an entire spectrum for the entire data length, implements multitaper frequency transformation.
cfg_freq.taper      = 'dpss';      % multiple tapers based on discrete prolate spheroidal sequences (DPSS), also known as the Slepian sequence.

cfg_freq.pad        ='nextpow2';
cfg_freq.foi        = 1:.5:100;    % frequencies of interest

%cfg_freq.t_ftimwin  =
cfg_freq.tapsmofrq  = tapsmofrq;  % number, the amount of spectral smoothing through multi-tapering. 
                                   % Note that 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.

cfg_freq.keeptrials  = 'yes';

[ft_freq_mtm] = ft_freqanalysis(cfg_freq, ft_TD_struct);

%% a "per day" spectrum day is found by the mean of all the sessions that 
% took place on a given day (between 0-3 sessions generally speaking)
fig_dir = fullfile(cfg.anal_dir,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset,cfg.pp_fft,'/');


cfg.plotted_time     = 'since_inital_session';
cfg.z_score          = true;
cfg.freq             = [0.5, 80];
cfg.y_lbl_txt        = ['Days since RCS implant', newline '(implant as day 0)'];

    plt_session_spectra(cfg, ft_freq_mtm.freq, pt_side_id, ft_freq_mtm.powspctrm, ...
                                    ch_names, par_db_oi_out, fig_dir);

cfg.plotted_time    = 'since_inital_session';
cfg.z_score         = false;
cfg.freq            = [0.5, 80];
cfg.pwr_limits      = [-10, -1];
cfg.y_lbl_txt       = ['Days since RCS implant', newline '(implant as day 0)'];

    plt_session_spectra(cfg, fft_bins_inHz, pt_side_id, pwrspectra_by_sess, ...
                                    ch_names, par_db_oi_out, fig_dir);

save([save_dir, pt_side_id, '_ch_names'], ...
'ch_names', '-v7.3');




end
