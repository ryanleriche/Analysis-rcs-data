function rcs_fft_peri_survey(cfg,  rcs_pp_dir, pt_side_id, par_db, REDcap)
%{


%}
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
%%
pwrspectra_by_sess  = nan(4, height(par_db_oi), cfg.fftSize/2+1);
freq_by_sess        = nan(height(par_db_oi),    cfg.fftSize/2+1);


for i_sess = 1 : height(par_db_oi)
    
    load([dir_ss_pt_side, par_db_oi.sess_name{i_sess},'/', 'rcs_streamed.mat']);
    
    try
        [rcs_streamed_oi, par_db_oi_out(i_sess,:), Fsample_TD] ...
            ...
            =  pull_continuous_TD_rcs_chunk(...
            ...
        cfg, redcap, rcs_streamed, par_db_oi(i_sess,:));
    catch
        continue
    end


    % actual fft size of device for 64, 250, 1024 fftpoints
    switch cfg.fftSize
        case 64,   fftSizeActual = 62;   n_pad   = 4;
        case 256,  fftSizeActual = 250;  n_pad   = 6;
        case 1024, fftSizeActual = 1000; n_pad   = 24;
    end
    
   
    % create Hann window points
    % number of samples since last FFT
    % samp_elapsed  = (samp_rate*interval/1e3);
    % % overlap of current window w/ last window
    % nonoverlap    = (samp_elapsed/fftSizeActual);
    % overlap       = 1-nonoverlap;
    hann_win        = hannWindow(cfg.fftSize, cfg.windowLoad);
    sampInfft       = Fsample_TD * (cfg.fftinterval/1000);

    i_ffts          = 1:sampInfft:(height(rcs_streamed_oi)-cfg.fftSize);

    df              = Fsample_TD/cfg.fftSize;

    freq_by_sess(i_sess, :)   = df * (0:(cfg.fftSize/2));

   for i_ch = 0:3
        % Extract neural channel
        lfp_mv      = rcs_streamed_oi.(['Ch_TD', num2str(i_ch)]);

        pwrspectra_by_time     = nan(length(i_ffts), cfg.fftSize/2+1);

        for i = 1 : length(i_ffts) 
            % zero pad samples of interest based on fftSize
            i_td                   = i_ffts(i) : i_ffts(i)+ fftSizeActual-1;
            tmp_td                 = [lfp_mv(i_td); zeros(n_pad,1)]';
        
            % Apply fft on Hann windowed signal 
                % Q: Should TD should be zero-padded or Hann windowed first?
            tmp_fft                      = fft(tmp_td .* hann_win, cfg.fftSize);
        
            % From double to single-sided FFT
            % divide bin centered on zero by two--ignore negative freq
            tmp_fft                      = [tmp_fft(1)/2, tmp_fft(2:cfg.fftSize/2+1)];

            % get power and normalize by fftSize
            pwrspectra_by_time(i, :)   = (abs(tmp_fft ./cfg.fftSize)).^2;
        end
        pwrspectra_by_sess(i_ch+1, i_sess, :) = mean(pwrspectra_by_time, 1, 'omitnan');
    end
end


%%% save parsed database of interest, the power spectrum per session, frequencies, and channel names
save_dir = fullfile(rcs_pp_dir,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset,cfg.pp_fft,'/');

% remove streaming sessions w/o REDcap survey
i_w_sur  = find(~cellfun(@isempty, par_db_oi_out.sess_name));

par_db_oi_out      = par_db_oi_out(i_w_sur,:);

pwrspectra_by_sess = pwrspectra_by_sess(:, i_w_sur, :);

fft_bins_inHz      = freq_by_sess(~isnan(freq_by_sess))';

ch_names            = parse_chan_names(pt_side_id, par_db_oi_out);



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
                                    ch_names, par_db_oi_out , fig_dir);


end