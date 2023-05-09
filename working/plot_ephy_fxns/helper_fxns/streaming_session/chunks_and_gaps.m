function [comb_dt_chunks, per_TD_lost] = chunks_and_gaps(comb_dt) 

    per_TD_lost         = sum(isnan(comb_dt.Ch_TD0)) * 100 / height(comb_dt);
    
    i_nan = isnan(comb_dt.Ch_TD0);
    % report gaps 
    diff_nans    = diff(i_nan);
    i_gap_end    = find(diff_nans == -1) + 1; 
    i_gap_start  = find(diff_nans == 1) + 1; 

    if i_nan(1) == 1 % if data start with gap 
        i_gap_start = [1; i_gap_start ];
    end
    if i_nan(end) == 1 % if data ends with gap
        i_gap_end = [i_gap_end; length(i_nan) ];
    end
    localTime = comb_dt.localTime;
    
    gaps = localTime(i_gap_end) - localTime(i_gap_start);
    gaps.Format = 'hh:mm:ss.SSSS';
    
    comb_dt_chunks.gap_avg_duration = mean(gaps);
    comb_dt_chunks.gap_max_duraiton = max(gaps);

    comb_dt_chunks.i_gap_start      = i_gap_start;
    comb_dt_chunks.i_gap_end        = i_gap_end;
    
    % report chunks
%     if i_nan(1) == 1 % if starts w/ gap, first "from gap" is first chunk start
%         from_gap = find(diff_nans == -1);
% 
%         diff_nans(from_gap(1)) = 1;
% 
%     end
    
    i_chunk_start  = find(diff_nans == -1) + 1;
    i_chunk_end    = find(diff_nans == 1); 
    
    if i_nan(1) == 0 % if data start w/ chunk (generally, should)
        i_chunk_start = [1; i_chunk_start];

    end
    if i_nan(end) == 0 % if data ends w/ chunk (generally, should)
        i_chunk_end   = [i_chunk_end ; length(i_nan) ];
    end
    
    chunks = localTime(i_chunk_end) - localTime(i_chunk_start);
    chunks.Format = 'hh:mm:ss.SSSS';
    
    comb_dt_chunks.chunk_avg_duration = mean(chunks);
    comb_dt_chunks.chunk_max_duration = max(chunks);
    
    comb_dt_chunks.idxchunkStart      = i_chunk_start;
    comb_dt_chunks.idxchunkEnd        = i_chunk_end;
end