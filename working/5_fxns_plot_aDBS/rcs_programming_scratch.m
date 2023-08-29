%% from Jackson's "FeatureOverview.csv' of primary and secondary pain 
% decoding features based on sliding AUC

% setting these top two features as channels' powerbands 0 and 1

% peaks from .csv
PbLeft_peak              = [27, 8.5, 26.5, 8.5, 22, 0, 7.5, 24.5];
PbRight_peak             = [25.5, 6.5, 24, 7.5, 17, 34.5, 25, 12];

% band width of 5 Hz surrounding peak frequency (i.e., 2.5 Hz on each side)
band_width              = 5;

%% find closest band start and stop WITHIN desired band width

% FFT bins are based off of time-domain sampling rate (Fsample_TD) and the
% size of the FFT window (fft.sizeInSamp)

Fsample_TD              = 250;
fft.sizeInSamp          = 256;

% preparing channel power-band names based on channels 0, 1, 2, and 3 each
% having power-band 0 and power-band 2 as inputs

ch_power_bands          = {'Ch0_Pb0', 'Ch0_Pb1', 'Ch1_Pb0', 'Ch1_Pb1',...
                           'Ch2_Pb0', 'Ch2_Pb1', 'Ch3_Pb0', 'Ch3_Pb1'};

PbLeft_tbl                = table;
PbRight_tbl               = table;

% where calculation per se takes place

for i = 1 : length(ch_power_bands)

    PbLeft_tbl.(ch_power_bands {i})  ...
        ...
        = rcs_power_band_bin(...
        ...
   Fsample_TD, fft, PbLeft_peak(i), band_width);


   PbRight_tbl.(ch_power_bands {i})  ...
        ...
        = rcs_power_band_bin(...
        ...
   Fsample_TD, fft, PbRight_peak(i), band_width);

end

% labeling rows for completeness
PbLeft_tbl.Properties.RowNames = {'peak', 'desired_width', 'start', 'end', 'actual_width'};
PbRight_tbl.Properties.RowNames = {'peak', 'desired_width', 'start', 'end', 'actual_width'};

%%
function band = rcs_power_band_bin(Fsample_TD, fft, band_peak, band_width)

% frequency resolution (df)
df                 = Fsample_TD/fft.sizeInSamp;

% SCBS shows half the bandwidth for 0 frequency since 0 Hz includes
% negative frequency (essentially cutting the power in half to be on the
% same scale as the rest of the positive frequencies)
freq_vec           = round(...
                            [0, (df/2:df: Fsample_TD/2-df )  ]...
                            ,2)';

band_start         = band_peak - band_width/2;
freq_diff          = freq_vec - band_start;

% band start is max frequency w/n bandwidth
band_start         = freq_vec(...
                            min(freq_vec(freq_diff > 0) - band_start ) == freq_diff...
                            );

band_end         = band_peak+band_width/2;
freq_diff        = freq_vec - band_end;

% band end is min frequency w/n bandwidth
band_end        = freq_vec(...
                            max(freq_vec(freq_diff < 0) - band_end ) == freq_diff...
                            );
% return peak, desired band width, band start, band end, and the actual
% band width
band = [band_peak; band_width; band_start; band_end; band_end - band_start];
end


%%




