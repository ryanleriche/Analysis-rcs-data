function [rcs_streamed_oi, par_db_row, Fsample_TD] ...
            =  pull_continuous_TD_rcs_chunk(...
            cfg, redcap, rcs_streamed, par_db_row)

%{
(1) find REDcap survey within streaming session

(2) see which sample point is closest to REDcap survey time

(3) find continuous TD data of specified duration (cfg.epoch_dur) w/n specified 
time window (cfg.r_cap_cent_toi)

(4) return continuous TD data; and REDcap survey index w/n TD data, explicit
epoch start and end (REDcap at 0 s) w/n parsed database
%}

i_r_cap       = find(isbetween(redcap.time,...
                              par_db_row.timeStart + cfg.r_cap_cent_toi(1),...
                              par_db_row.timeStop)...
                        );

if ~isempty(i_r_cap)
    % sample with survey
    [~, i_samp_r_cap] ...
        = min(abs(redcap.time(i_r_cap(end)) - rcs_streamed.TimeInPST));
    
    Fsample_TD     = par_db_row.TDsampleRates;
    
    % search for continuous epoch closest to epoch of interest
    samp_back   = seconds(cfg.r_cap_cent_toi(1)) * Fsample_TD;
    samp_ahead  = seconds(cfg.r_cap_cent_toi(2)) * Fsample_TD;
    
    samp_N = seconds(cfg.epoch_dur) * Fsample_TD;
    
    i_samp_epoch  = i_samp_r_cap - (samp_back)...
                    :...
                    i_samp_r_cap-1+samp_ahead;
    
    
    td_vec = rcs_streamed.Ch_TD0(i_samp_epoch);
    
    [td_chunk.length, td_chunk.vec_ident,...
     td_chunk.starts, td_chunk.ends] ...
        = SplitVec(...
            ~isnan(td_vec), [],'length','firstval','first','last');
    
    
    chunk_candi = find(td_chunk.length >= samp_N & td_chunk.vec_ident);
    
    if ~isempty(chunk_candi)

        switch cfg.epoch_nearest_to
            case 'end'
        
                chunk_start = td_chunk.ends(chunk_candi(end))-samp_N+1;
                chunk_end   = td_chunk.ends(chunk_candi(end));
        
            case 'start'
        
        end
        
        i_samp_oi       = i_samp_epoch(chunk_start: chunk_end);
        
        rcs_streamed_oi                    = rcs_streamed(i_samp_oi, :);
    
        par_db_row.index_samp_rcap_wn_TD  = i_samp_r_cap;
        par_db_row.epoch_start_rcap_at_0s = (i_samp_oi(1) - i_samp_r_cap) /Fsample_TD;
        par_db_row.epoch_end_rcap_at_0s   = (i_samp_oi(end) - i_samp_r_cap) /Fsample_TD;
        par_db_row.survey_time            = redcap.time(i_r_cap(end));

        return
    end
end

rcs_streamed_oi = {};       Fsample_TD = {};

par_db_row.index_samp_rcap_wn_TD  =  nan;
par_db_row.epoch_start_rcap_at_0s =  nan;
par_db_row.epoch_end_rcap_at_0s   =  nan;
par_db_row.survey_time            =  NaT;

par_db_row.survey_time.TimeZone = redcap.time.TimeZone;
