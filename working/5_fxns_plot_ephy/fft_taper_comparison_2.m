function fft_taper_comparison_2(cfg, dirs, pt_side_id, ft_form_TD)
%{
    
N minutes before REDcap survey, for those epochs w/ at least 3 minutes of
data
->
in 1 s, non-overlapping windows, run MTM per window, and average per
session
->


%}

% based off of Jeremy's MTM analysis of stage 0 spectra
K_tapers   = 3;
time_width = 2;

tapsmofrq = (K_tapers +1) / (time_width *seconds(cfg.cont_chunk_dur));


% pull out FieldTrip formatted time-domain data, and specifically its
% corresponding parsed database
ft_TD_struct = ft_form_TD.(pt_side_id);
par_db        = ft_TD_struct.rcs.par_db;


% use parsed data base to return surveys w/:
    % N minutes of data before REDcap (cfg.min_pre_rcap)
    % w/ 60% of data

N_samp_whole_win   = seconds(cfg.min_pre_rcap) * ft_TD_struct.fsample;
i_rcap_oi          = (par_db.epoch_samp_rcap_wn_TD-1 >= N_samp_whole_win);

ft_TD_struct.time  = ft_TD_struct.time(i_rcap_oi);
ft_TD_struct.trial = ft_TD_struct.trial(i_rcap_oi);
par_db             = par_db(i_rcap_oi, :);

cfg.foi            = 0 ...
                     :0.5:...
                     ft_TD_struct.fsample/2; % frequency band of interest

%%
Nsamp_slide_win    = seconds(cfg.cont_chunk_dur) * ft_TD_struct.fsample;
t_vec              = 0:(1/ft_TD_struct.fsample):seconds(cfg.cont_chunk_dur);

freq_mat_initalized = false;

for j = 1 : height(par_db)

    td_chunk     = par_db.td_chunk{j};

    % only use samples from window of interest before REDcap survey
    i_samp_start = td_chunk.ends(end) - N_samp_whole_win;
    chunk_lat    = i_samp_start - td_chunk.starts;

    % finds window start w/n chunk
    samp_after_start  = min(chunk_lat(chunk_lat > 0));

    % replaces that chunk start WITH the window of interest before REDcap
    % survey
    td_chunk.starts(chunk_lat == samp_after_start) = i_samp_start;

    % recalculate td chunk length since start of one has been changed
    td_chunk.length = td_chunk.ends- td_chunk.starts +1;

    % remove any TD occuring before window of interest starts
    td_chunk(chunk_lat > samp_after_start, :) = [];

    whole_dur = sum(td_chunk.length);

    % only pass TD without packet loss AND large enough chunk length
    td_chunk     = td_chunk(...
                            strcmp(td_chunk.label, 'TD_chunk') &...
                            td_chunk.length >= Nsamp_slide_win...
                            , :);
 
    % only analyze is 60% of data w/n window of interest is available
    prop_avail_TD = sum(...
                        td_chunk.length(strcmp(td_chunk.label, 'TD_chunk'))...
                        ) ...
                        / whole_dur;

    if prop_avail_TD <= 0.6;        continue;       end
       
    ss_td        = ft_TD_struct.trial{j};
    
    % break TD chunks into sub-chunks of identical length
    cnt = 1;
    
    temp_td         = struct;
    temp_td.dimord  = ft_TD_struct.dimord;
    temp_td.label   = ft_TD_struct.label;
    temp_td.fsample = ft_TD_struct.fsample;

    % format into cell-array of sub-chunks of chan by time as FieldTrip
    % struct
   for h = 1:height(td_chunk)
       
       win_ind = td_chunk.starts(h)...
                        :Nsamp_slide_win:...
                 td_chunk.ends(h);

       for k = 1 : length(win_ind) - 1
            temp_td.sampleinfo(cnt, 1) = win_ind(k);
            temp_td.sampleinfo(cnt, 2) = win_ind(k+1);


            temp_td.trial{cnt} = ss_td(:, win_ind(k) : win_ind(k+1));
            temp_td.time{cnt}  = t_vec;
            cnt = cnt + 1;
       end
   end

    cfg_freq            = [];
    cfg_freq.output     = 'pow';       % return the power-spectra
    cfg_freq.method     = 'mtmfft';    % analyses an entire spectrum for the entire data length, implements multitaper frequency transformation.
    cfg_freq.taper      = 'dpss';      % multiple tapers based on discrete prolate spheroidal sequences (DPSS), also known as the Slepian sequence.
    
    cfg_freq.pad        ='nextpow2';
    cfg_freq.foi        = cfg.foi;    % frequencies of interest
    
    cfg_freq.toi         = 'all';
    cfg_freq.tapsmofrq   = tapsmofrq;
    cfg_freq.polyremoval = -1;

        temp_freq_mtm = ft_freqanalysis(cfg_freq, temp_td);


    if freq_mat_initalized == false
        powspctrm_mtm = nan*zeros(height(par_db), length(ft_TD_struct.label),length(temp_freq_mtm.freq));
        powspctrm_han = powspctrm_mtm;

        freq_mat_initalized = true;
    end

    powspctrm_mtm(j, :, :) = temp_freq_mtm.powspctrm;


    cfg_freq.taper      = 'hanning';
        temp_freq_han = ft_freqanalysis(cfg_freq, temp_td);

    powspctrm_han(j, :, :) = temp_freq_han.powspctrm;
end

i_insuff_td   = any(isnan(powspctrm_han), [2,3]);
%i_no_td       = any(powspctrm_han ==0, [2,3]);

ft_freq_mtm_by_win              = temp_freq_mtm;
ft_freq_mtm_by_win.powspctrm    = powspctrm_mtm(~i_insuff_td, : ,:);
ft_freq_mtm_by_win.dimord       = 'rpt_chan_freq';

ft_freq_mtm_by_win.rcs.td_cfg   = ft_TD_struct.rcs.td_cfg;
ft_freq_mtm_by_win.rcs.par_db   = par_db(~i_insuff_td, :);
%%%
ft_freq_han_by_win              = ft_freq_mtm_by_win;
ft_freq_han_by_win.powspctrm    = powspctrm_han(~i_insuff_td, : ,:);


%%
fig_dir = fullfile(cfg.anal_dir,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset,cfg.pp_fft);

cfg.z_score         = false;
cfg.freq            = [1, 100];
cfg.pwr_limits      = [-11, -1];
cfg.title           = 'taper_comparison';

cfg.txt_leg         = ["FieldTrip's MTM", "FieldTrip's Hanning"];

    plt_rcs_PSDs(cfg, pt_side_id, fig_dir, ft_freq_mtm_by_win, ft_freq_han_by_win)


ft_fft_dir  = fullfile(dirs.rcs_preproc,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset, cfg.pp_fft);

if ~isfolder(ft_fft_dir);       mkdir(ft_fft_dir);      end

save(fullfile(ft_fft_dir, [pt_side_id,'_ft_form_FFT_struct.mat']), 'ft_freq_mtm_by_win', '-v7.3');

end
