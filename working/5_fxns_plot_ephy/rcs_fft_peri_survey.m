function rcs_fft_peri_survey(cfg,  rcs_pp_dir, pt_side_id, par_db, REDcap)
%{


%}
\

cfg.min_sess_dur    = duration('00:05:00');

%%%
cfg.fftSize         = 1024;             % N samples to use in FFT
cfg.windowLoad      = '100% Hann';      % 100 percent hanning window
cfg.fftinterval     = 500;              % how often in ms a FFT should be calculated


cfg.r_cap_cent_toi = [duration('0:02:00'), duration('0:00:00')];

cfg.epoch_dur      = duration('0:00:30');

cfg.epoch_nearest_to = 'end';
%%

%%
dir_ss_pt_side  = fullfile(rcs_pp_dir,'streaming_sessions', pt_side_id, cfg.pp_RCS_TD_subset, '/');

if ~isfolder(dir_ss_pt_side);      error('%s | path to TD data was not found', pt_side_id);    end


db_RCSXXX        = par_db.(pt_side_id);
save_dir_struct  = dir([dir_ss_pt_side, '*Session*']);

par_db_oi  = db_RCSXXX(...
                    ge(db_RCSXXX.duration, cfg.min_sess_dur) &...
                    ismember(db_RCSXXX.sess_name, {save_dir_struct.name}) &...
                    db_RCSXXX.TDsampleRates == 250 ...
                    ,:);

par_db_oi_out = table;
redcap        = REDcap.(pt_side_id(1:end-1));

%% circumvent FieldTrip reading functions by matching output of their prepreprocessing
% https://www.fieldtriptoolbox.org/faq/how_can_i_import_my_own_dataformat/

data = struct;
data.dimord = 'chan_time';

cnt = 1;

for i_sess = 1 : height(par_db_oi)
    
    load([dir_ss_pt_side, par_db_oi.sess_name{i_sess},'/', 'rcs_streamed.mat']);
    
    try
        [rcs_streamed_oi, par_db_oi_out(i_sess,:), Fsample_TD] ...
            ...
            =  pull_continuous_TD_rcs_chunk(...
            ...
        cfg, redcap, rcs_streamed, par_db_oi(i_sess,:));

        
    catch;  continue;  
    end

    if cnt == 1
        data.fsample = Fsample_TD;              % sampling frequency in Hz, single number
        data.label   = compose('Ch_TD%g',0:3)'; % cell-array containing strings, Nchan*1
    end

    data.trial{1, cnt}  ...                   % cell-array containing a data matrix for each
        = rcs_streamed_oi{:, data.label}';    % trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix

    data.time{1, cnt}  ...                    % cell-array containing a time axis for each
         = seconds(...
             rcs_streamed_oi.TimeInPST -...
             rcs_streamed_oi.TimeInPST(1))';                % trial (1*Ntrial), each time axis is a 1*Nsamples vector

    data.sampleinfo(cnt, 1:2) = [data.time{1, cnt}(1), data.time{1, cnt}(end)];
    cnt = cnt+1;
end

%%
K_tapers = 3;
time_width =2;

freq_width = (K_tapers +1) / (2 * time_width);

cfg_freq = [];
cfg_freq.output     = 'pow';       % return the power-spectra
cfg_freq.method     = 'mtmfft';    % analyses an entire spectrum for the entire data length, implements multitaper frequency transformation.
cfg_freq.taper      = 'dpss';      % multiple tapers based on discrete prolate spheroidal sequences (DPSS), also known as the Slepian sequence.


cfg_freq.pad        ='nextpow2';
cfg_freq.foi        = 1:.5:100;    % frequencies of interest

%cfg_freq.t_ftimwin  =
cfg_freq.tapsmofrq  = freq_width;  % number, the amount of spectral smoothing through multi-tapering. 
                                   % Note that 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.

cfg_freq.keeptrials  = 'yes';

[ft_freq_struct] = ft_freqanalysis(cfg_freq, data);

%%
pwrspectra_by_sess  = ft_freq_struct.powspctrm;
fft_bins_inHz       = ft_freq_struct.freq;

% remove streaming sessions w/o REDcap survey
i_w_sur            = find(~cellfun(@isempty, par_db_oi_out.sess_name));
par_db_oi_out      = par_db_oi_out(i_w_sur,:);

ch_names = parse_chan_names(pt_side_id, par_db_oi_out);

%%% save parsed database of interest, the power spectrum per session, frequencies, and channel names
save_dir = fullfile(rcs_pp_dir,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset,cfg.pp_fft,'/');

if ~isfolder(save_dir);     mkdir(save_dir);    end

save([save_dir, pt_side_id, '_pwrspectra_by_sess'], ...
    'pwrspectra_by_sess', '-v7.3');

fft_bins_inHz = unique(fft_bins_inHz);

save([save_dir, pt_side_id, '_fft_bins_inHz'], ...
    'fft_bins_inHz', '-v7.3');

save([save_dir, pt_side_id, '_ch_names'], ...
    'ch_names', '-v7.3');


writetable(par_db_oi_out, [save_dir,pt_side_id, '_parsed_db_oi.xlsx']);


%% a "per day" spectrum day is found by the mean of all the sessions that 
% took place on a given day (between 0-3 sessions generally speaking)
fig_dir = fullfile(cfg.anal_dir,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset,cfg.pp_fft,'/');


cfg.plotted_time     = 'since_inital_session';
cfg.z_score          = true;
cfg.freq             = [0.5, 80];
cfg.y_lbl_txt        = ['Days since RCS implant', newline '(implant as day 0)'];

    plt_session_spectra(cfg, fft_bins_inHz, pt_side_id, pwrspectra_by_sess, ...
                                    ch_names, par_db_oi_out, fig_dir);

cfg.plotted_time    = 'since_inital_session';
cfg.z_score         = false;
cfg.freq            = [0.5, 80];
cfg.pwr_limits      = [-10, -1];
cfg.y_lbl_txt       = ['Days since RCS implant', newline '(implant as day 0)'];

    plt_session_spectra(cfg, fft_bins_inHz, pt_side_id, pwrspectra_by_sess, ...
                                    ch_names, par_db_oi_out, fig_dir);
end