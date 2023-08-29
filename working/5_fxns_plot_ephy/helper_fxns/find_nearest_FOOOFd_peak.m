function peak_freq_tbl = find_nearest_FOOOFd_peak(fooof_struct_1, fooof_struct_2)


peak_freq_tbl              = table;
peak_freq_tbl.cent_freq_1  = {};
peak_freq_tbl.cent_freq_2  = {};
peak_freq_tbl.bf10         = {};
peak_freq_tbl.p_value      = {};

peak_freq_tbl.bf01         = {};


for i_chan = 1 : height(fooof_struct_1.sense_chan)

    chan_oi_1           = fooof_struct_1.sense_chan(i_chan, :);
    [~, i_peak_sort]    = sort(chan_oi_1.mean_fooof_params.peak_power{1}, 'descend');


    for j_peak = 1 : length(i_peak_sort)

        i_peak_oi = i_peak_sort(j_peak);

        peak_freq         = chan_oi_1.mean_fooof_params.center_freq{1}(i_peak_oi);
        peak_bw           = chan_oi_1.mean_fooof_params.band_width{1}(i_peak_oi);
    
        all_cent_freqs    = chan_oi_1.fooof_params{1}.center_freq;
        peak_freqs_1      = nan(length(all_cent_freqs), 1);

        for j_sess = 1 :length(all_cent_freqs)
    
            j_freqs = all_cent_freqs{j_sess};
    
            if ~isempty(j_freqs)
                [~, i_min]            = min(abs(peak_freq - j_freqs));
                peak_freqs_1(j_sess) = j_freqs(i_min);
            end
    
    
        end

        peak_freqs_1(...
            peak_freqs_1 < (peak_freq -peak_bw/2) ...
            |...
            peak_freqs_1 > (peak_freq + peak_bw/2)...
            )...
            = nan;
        
         peak_freqs_1(isnan(peak_freqs_1)) = [];
        
        %%%
        chan_oi_2       = fooof_struct_2.sense_chan(i_chan, :);
        all_cent_freqs    = chan_oi_2.fooof_params{1}.center_freq;
        peak_freqs_2    = nan(length(all_cent_freqs), 1);
    
        for j_sess = 1 :length(all_cent_freqs)
    
            j_freqs = all_cent_freqs{j_sess};
    
            if ~isempty(j_freqs)
                [~, i_min]            = min(abs(peak_freq - j_freqs));
                peak_freqs_2(j_sess) = j_freqs(i_min);
            end
        end
    
        peak_freqs_2(...
            peak_freqs_2 < (peak_freq -peak_bw/2)...
            |...
            peak_freqs_2 > (peak_freq + peak_bw/2)...
            )...
            = nan;
    
        peak_freqs_2(isnan(peak_freqs_2)) = [];
    
        %%%
        peak_freq_tbl.cent_freq_1{i_chan}{j_peak}  = peak_freqs_1;
        peak_freq_tbl.cent_freq_2{i_chan}{j_peak}  = peak_freqs_2;
        
    
        [peak_freq_tbl.bf10{i_chan}(j_peak),...
         peak_freq_tbl.p_value{i_chan}(j_peak)] ...
         ...
         = bf.ttest2(...
         ...
         peak_freqs_2, peak_freqs_1, 'scale', 5);
    end

    peak_freq_tbl.bf01{i_chan}      = 1 ./ peak_freq_tbl.bf10{i_chan};
end

peak_freq_tbl.labels_2   = fooof_struct_2.sense_chan.Properties.RowNames;
peak_freq_tbl.labels_1   = fooof_struct_1.sense_chan.Properties.RowNames;



end

%%