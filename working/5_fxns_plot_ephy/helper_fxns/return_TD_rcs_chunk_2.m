function [rcs_streamed_oi, par_db_row, fsample_TD] ...
            =  return_TD_rcs_chunk_2(...
            cfg, i_rcap_wn_ss, i_rcap, redcap, rcs_streamed, par_db_row)

%{
(1) find nearest sample point to REDcap survey time

(2) find continuous TD data of specified duration (cfg.epoch_dur) w/n specified 
time window (cfg.r_cap_cent_toi)

(3) return continuous TD data; and REDcap survey index w/n TD data, explicit
epoch start and end (REDcap at 0 s) w/n parsed database
%}
rcap_time = redcap.time(i_rcap_wn_ss(i_rcap));

% sample with survey
[~, i_samp_r_cap] ...
    = min(abs(rcap_time - rcs_streamed.TimeInPST));

% search for continuous epoch closest to epoch of interest
fsample_TD    = par_db_row.TDsampleRates;
samp_back     = seconds(cfg.r_cap_cent_toi(1)) * fsample_TD;
samp_ahead    = seconds(cfg.r_cap_cent_toi(2)) * fsample_TD;

% samp_N        = seconds(cfg.epoch_dur) * fsample_TD;

% search for time before
i_samp_epoch  = i_samp_r_cap - (samp_back):i_samp_r_cap-1+samp_ahead;

% if no time before, return samples from session start
i_samp_oi     = i_samp_epoch(i_samp_epoch >0);
td_vec        = rcs_streamed.Ch_TD0(i_samp_oi);

% explicitly save out epochs of TD chunks, and packet loss (i.e., mV sample
% points vs NaNs)
td_chunk      = table;

[td_chunk.length, td_chunk.vec_ident,...
 td_chunk.starts, td_chunk.ends] ...
    = SplitVec(...
        ~isnan(td_vec), [],'length','firstval','first','last');


td_chunk.label(td_chunk.vec_ident)     = {'TD_chunk'};
td_chunk.label(~td_chunk.vec_ident)    = {'packet_loss'};

i_chunk = strcmp(td_chunk.label, 'TD_chunk');

%%%
par_db_row.td_chunk                 = {td_chunk};
par_db_row.percent_TD_streamed      = 100*sum(td_chunk.length(i_chunk)) / sum(td_chunk.length);


if ~isempty(i_samp_oi)
    rcs_streamed_oi                    = rcs_streamed(i_samp_oi, :);
    
    par_db_row.epoch_samp_rcap_wn_TD  = i_samp_r_cap;
    par_db_row.epoch_start_rcap_at_0s = (i_samp_oi(1) - i_samp_r_cap) /fsample_TD;
    par_db_row.epoch_end_rcap_at_0s   = (i_samp_oi(end) - i_samp_r_cap) /fsample_TD;
    
    par_db_row.epoch_durationInSecs   = length(i_samp_oi) / fsample_TD;
    par_db_row.timeSurvey             = rcap_time;
    
    
    par_db_row.rcap_indices           = i_rcap_wn_ss(i_rcap);
    par_db_row.rcap_latency_timeStart = par_db_row.rcap_latency_timeStart{1}(i_rcap);


    return
end
 


rcs_streamed_oi = {};       fsample_TD = {};

par_db_row.index_samp_rcap_wn_TD  =  nan;
par_db_row.epoch_start_rcap_at_0s =  nan;
par_db_row.epoch_end_rcap_at_0s   =  nan;
par_db_row.epoch_durationInSecs   =  nan;
par_db_row.timeSurvey             =  NaT;

par_db_row.timeSurvey.TimeZone     = rcap_time.TimeZone;

par_db_row.rcap_indices            = nan;
par_db_row.rcap_latency_timeStart  = duration('');
end

