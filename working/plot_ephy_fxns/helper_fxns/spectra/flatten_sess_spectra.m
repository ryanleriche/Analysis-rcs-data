function   [pwrspectra_by_sess_out, art_fft_peaks_tbl]...
            ...
            = flatten_sess_spectra(...
            ...
            cfg, fft_bins_inHz, pt_side_id, pwrspectra_by_sess, ch_names, par_db_oi, save_dir)


%%

i_freq_oi     = find((fft_bins_inHz >= cfg.freq(1) & fft_bins_inHz <= cfg.freq(2)));

t_start = dateshift( par_db_oi.timeStart(1), 'start', 'day');
t_end   = dateshift( par_db_oi.timeStart(end), 'end', 'day');

%%
figure('Units','inches','Position',[0, 3, 12, 8]);
tiledlayout(2,2, 'Padding','compact')

sgtitle(sprintf('%s | %s–%s, Stage 1 (prior to analgesic stim testing)\n%g mean of sessions',...
            pt_side_id, datestr(t_start, 'dd-mmm-yyyy'), datestr(t_end, 'dd-mmm-yyyy'), height(par_db_oi)),...
            'FontSize', 12)

art_fft_peaks_tbl =table;

for i_ch = 1:4

    nexttile
    plt_ch_name = unique(ch_names(:, i_ch));

    title(sprintf('Ch%g | %s', i_ch - 1, plt_ch_name{:}));

    plt_spectrum =  mean(...
                        log10(squeeze(pwrspectra_by_sess(i_ch, :, i_freq_oi))),...
                        'omitnan');

    input = {plt_spectrum,fft_bins_inHz(i_freq_oi), 'MinPeakProminence',.25,'Annotate','extents'};

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


exportgraphics(gcf, [save_dir, '3_art_sess_power_spectra_peaks.png'])

%%
fft_bin_size          = min(diff(fft_bins_inHz));
freq_oi               = fft_bins_inHz(i_freq_oi);

pwrspectra_by_sess_out = pwrspectra_by_sess;
for i_ch = 1:4

    plt_spectrum =  log10(squeeze(pwrspectra_by_sess(i_ch, :, i_freq_oi)));

    tmp_spectrum = plt_spectrum;
                     
    pks_freq     = art_fft_peaks_tbl.pks_freq{i_ch};

    if ~isempty(pks_freq)

        pks_width = art_fft_peaks_tbl.half_prom_freq_width{i_ch};

        for i_pk = 1 :length(pks_freq)

            n_freq_pad = ceil(pks_width(i_pk) *3 / fft_bin_size);

            min_freq_art   = pks_freq(i_pk) -  fft_bin_size*n_freq_pad ;
            max_freq_art   = pks_freq(i_pk) + fft_bin_size*n_freq_pad ;
            tmp_blank_freq = min_freq_art:fft_bin_size:max_freq_art;


            i_physio = ~ismember(freq_oi, tmp_blank_freq);
            
            for i_spec = 1 : size(plt_spectrum,1)

                tmp_spectrum(i_spec, ~i_physio)...
                    ...
                    = interp1(...
                    ...
                        freq_oi(i_physio), ...
                        plt_spectrum(i_spec,i_physio),...
                        tmp_blank_freq);


            end
        end
    end
pwrspectra_by_sess_out(i_ch, :, i_freq_oi) = 10.^tmp_spectrum;
end


end