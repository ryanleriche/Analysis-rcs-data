function  i_ch_sess_rej = plt_rej_sess_spectra(cfg, fft_bins_inHz, pt_side_id, pwrspectra_by_sess, ...
                                    ch_names, par_db_oi, save_dir)


%%

i_freq_oi     = find((fft_bins_inHz >= cfg.freq(1) & fft_bins_inHz <= cfg.freq(2)));

t_start = dateshift( par_db_oi.timeStart(1), 'start', 'day');
t_end   = dateshift( par_db_oi.timeStart(end), 'end', 'day');
t_vec   = par_db_oi.timeStart;

c_str = "log_{10}(mV^{2}/Hz)";
preset_colormap(@brewermap, "YlOrRd");

%%

cfg.y_lbl_txt  = 'individual RCS streaming sessions';
i_ch_sess_rej  = zeros(size(pwrspectra_by_sess, [1,2]))';

%%%
figure('Units','inches','Position',[0, 3, 18, 12]);

tiledlayout(2,2, 'Padding','tight')

sgtitle(sprintf('%s | %s–%s, Stage 1 (prior to analgesic stim testing)\n%g sessions total',...
    pt_side_id, datestr(t_start, 'dd-mmm-yyyy'), datestr(t_end, 'dd-mmm-yyyy'), height(par_db_oi)),...
    'FontSize', 24)

for i_ch = 1:4
    nexttile
    plt_ch_name = unique(ch_names(:, i_ch));

    title(sprintf('Ch%g | %s', i_ch - 1, plt_ch_name{:}));

    plt_spectrum =  log10(squeeze(pwrspectra_by_sess(i_ch, :, i_freq_oi)));

    plt_t_vec =  t_vec;

    h_sur = surface(fft_bins_inHz(i_freq_oi), plt_t_vec, plt_spectrum);
    
    format_spectra(cfg, h_sur);
    
    %% on 2nd y axis plot power spectra on top of each other (rather than in spectrogram)
    yyaxis right;         set(gca, 'YColor', 'k'); yticklabels([]);     hold on
    plot(fft_bins_inHz(i_freq_oi), plt_spectrum,  '-k', 'LineWidth', .5, 'Marker','none');

    ylim(cfg.pwr_limits);

    %% select power and freq critera for artifact rejection
    disp('Select a rectangle (min power and freq value) to reject spectra')
    
    rect =  drawrectangle;
    
    %%% reject any session with given freq w/ power above visually
    %%% inspected threshold
    freq_min    = rect.Position(1);
    freq_max    =  freq_min + rect.Position(3);

    pwr_min     = rect.Position(2);
    pwr_max     = pwr_min + rect.Position(4);

    i_rej_foi   = fft_bins_inHz(i_freq_oi) >  freq_min & ...
                  fft_bins_inHz(i_freq_oi) <  freq_max;  
    i_sess_rej  =  any(...
                            plt_spectrum(:, i_rej_foi) > pwr_min & ...
                            plt_spectrum(:, i_rej_foi) < pwr_max...
                        , 2);

    hold on
    if any(i_sess_rej)
        plot(fft_bins_inHz(i_freq_oi), plt_spectrum(i_sess_rej, :), 'LineWidth', .5, 'Marker','none');
    end

    i_ch_sess_rej(:, i_ch) = i_sess_rej;

    if ~isfolder(save_dir);    mkdir(save_dir);  end

end

exportgraphics(gcf, [save_dir, '2_art_rejected_sess_power_spectra.png'])

%%

figure('Units','inches','Position',[0, 3, 18, 12]);

tiledlayout(2,2, 'Padding','tight')

sgtitle(sprintf('%s | %s–%s, Stage 1 (prior to analgesic stim testing)\n%g sessions total',...
    pt_side_id, datestr(t_start, 'dd-mmm-yyyy'), datestr(t_end, 'dd-mmm-yyyy'), sum(any(~i_ch_sess_rej,2))),...
    'FontSize', 24)

for i_ch = 1:4
    nexttile
    plt_ch_name = unique(ch_names(:, i_ch));

    title(sprintf('Ch%g | %s', i_ch - 1, plt_ch_name{:}));

    spectra_oi = squeeze(pwrspectra_by_sess(i_ch, :, i_freq_oi));
    
    spectra_oi(i_ch_sess_rej(:, i_ch)==1, :) = nan;

    plt_spectrum =  log10(spectra_oi);

    plt_t_vec =  t_vec;

    h_sur = surface(fft_bins_inHz(i_freq_oi), plt_t_vec, plt_spectrum);
    
    format_spectra(cfg, h_sur);
    
    %% on 2nd y axis plot power spectra on top of each other (rather than in spectrogram)
    yyaxis right;         set(gca, 'YColor', 'k'); yticklabels([]);     hold on
    plot(fft_bins_inHz(i_freq_oi), plt_spectrum,  '-k', 'LineWidth', .5, 'Marker','none');

    ylim(cfg.pwr_limits);

end

exportgraphics(gcf, [save_dir, '2_pruned_sess_power_spectra.png'])


%%
function format_spectra(cfg, h_sur)
    hold on 

    ylabel(cfg.y_lbl_txt);          xlabel('Frequency (Hz)');
   
    h_sur.LineStyle = 'none';

    colormap(preset_colormap);      
    cb_hdl = colorbar;              cb_hdl.Limits          = cfg.pwr_limits; 
    
    set(gca, 'FontSize', 14);       cb_hdl.Label.String    = c_str;
    cb_hdl.FontSize     = 10;       cb_hdl.Label.FontSize  = 18;

    ax = gca;               ax.YDir = 'Reverse';
end
end