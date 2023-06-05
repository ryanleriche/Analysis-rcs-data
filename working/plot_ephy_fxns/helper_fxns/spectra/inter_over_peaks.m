function   ft_struct  = inter_over_peaks(cfg, save_dir, pt_side_id, ft_struct)


%%

i_freq_oi     = find((ft_struct.freq >= cfg.freq(1) & ft_struct.freq <= cfg.freq(2)));

t_start = dateshift( ft_struct.rcs.par_db.timeStart(1), 'start', 'day');
t_end   = dateshift( ft_struct.rcs.par_db.timeStart(end), 'end', 'day');

%%
close all
set(0,'DefaultFigureVisible','on')

figure('Units','inches','Position',[0, 3, 12, 8]);
tiledlayout(2,2, 'Padding','compact')

sgtitle(sprintf('%s | %s–%s, Stage 1 (prior to analgesic stim testing)\n%g mean of sessions',...
            pt_side_id, datestr(t_start, 'dd-mmm-yyyy'), datestr(t_end, 'dd-mmm-yyyy'), height(ft_struct.rcs.par_db)),...
            'FontSize', 12)

art_fft_peaks_tbl =table;

for i_ch = 1:4

    nexttile
    title(sprintf('Ch%g | %s', i_ch - 1, ft_struct.label{i_ch}));

    plt_spectrum =  mean(...
                        log10(squeeze(ft_struct.powspctrm(:, i_ch, i_freq_oi))),...
                        'omitnan');

    input = {plt_spectrum,ft_struct.freq(i_freq_oi), 'MinPeakProminence',.3,'Annotate','extents'};

    [art_fft_peaks_tbl.pks_pwr{i_ch},...
     art_fft_peaks_tbl.pks_freq{i_ch},...
     art_fft_peaks_tbl.half_prom_freq_width{i_ch},...
     art_fft_peaks_tbl.prom_pwr{i_ch}] ...
     ...
        =  findpeaks( input{:});

    %%% nice built-in plotting call
    findpeaks(input{:});

    ylabel('log_{10}(mV^{2}/Hz)');      xlabel('Frequency (Hz)');

end
art_fft_peaks_tbl.Properties.RowNames =compose('Ch%g', 0:3);

ft_struct.rcs.art_fft_peaks_tbl = art_fft_peaks_tbl;

exportgraphics(gcf, fullfile(save_dir, '3_art_spectra_peaks.png'))

%%
fft_bin_size          = min(diff(ft_struct.freq));
freq_oi               = ft_struct.freq(i_freq_oi);

for i_ch = 1:4

    plt_spectrum =  log10(squeeze(ft_struct.powspctrm(:, i_ch, i_freq_oi)));

    tmp_spectrum = plt_spectrum;
                     
    pks_freq     = art_fft_peaks_tbl.pks_freq{i_ch};

    if ~isempty(pks_freq)

        pks_width = art_fft_peaks_tbl.half_prom_freq_width{i_ch};

        for i_pk = 1 : length(pks_freq)

            n_freq_pad     = ceil(pks_width(i_pk) *3 / fft_bin_size);

            min_freq_art   = pks_freq(i_pk) -  n_freq_pad;
            max_freq_art   = pks_freq(i_pk) + n_freq_pad;

            tmp_blank_freq = freq_oi(freq_oi >= min_freq_art & freq_oi <= max_freq_art);

            i_physio = ~ismember(freq_oi, tmp_blank_freq);
            
            % ensure peak has neighboring frequencies to interpolate over
            if sum(i_physio) > 1
                for i_spec = 1 : size(plt_spectrum,1)
    
                    tmp_spectrum(i_spec, ~i_physio)...
                        ...
                        = interp1(...
                            freq_oi(i_physio), ...
                            plt_spectrum(i_spec,i_physio),...
                            tmp_blank_freq);
                end
            end
        end
    end
ft_struct.powspctrm(:, i_ch, i_freq_oi) = 10.^tmp_spectrum;
end


end