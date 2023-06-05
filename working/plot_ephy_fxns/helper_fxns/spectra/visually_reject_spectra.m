function  i_ch_sess_rej = visually_reject_spectra(cfg, save_dir, pt_side_id, ft_struct)


%%
i_freq_oi     = find((ft_struct.freq >= cfg.freq(1) & ft_struct.freq <= cfg.freq(2)));

t_vec   = ft_struct.rcs.par_db.timeStart;

t_start = dateshift(t_vec(1), 'start', 'day');
t_end   = dateshift(t_vec(end), 'end', 'day');

preset_colormap(@brewermap, "YlOrRd");

%%
set(0,'DefaultFigureVisible','on')

cfg.y_lbl_txt  = 'individual RCS streaming sessions';
i_ch_sess_rej  = zeros(size(ft_struct.powspctrm, [1,2]));

%%%
figure('Units','inches','Position',[0, 3, 18, 12]);

tiledlayout(2,2, 'Padding','tight')

sgtitle(sprintf('%s | %s–%s, Stage 1 (prior to analgesic stim testing)\n%g sessions total',...
    pt_side_id, datestr(t_start, 'dd-mmm-yyyy'), datestr(t_end, 'dd-mmm-yyyy'), length(t_vec)),...
    'FontSize', 24)

for i_ch = 1:4
    nexttile

    title(sprintf('Ch%g | %s', i_ch - 1, ft_struct.label{i_ch}));

    plt_spectrum =  log10(squeeze(ft_struct.powspctrm(:,i_ch, i_freq_oi)));

    plt_t_vec =  t_vec;
    h_sur     = surface(ft_struct.freq(i_freq_oi), plt_t_vec, plt_spectrum);
    
    format_spectra(cfg, h_sur);
    
    %% on 2nd y axis plot power spectra on top of each other (rather than in spectrogram)
    yyaxis right;         set(gca, 'YColor', 'k'); yticklabels([]);     hold on
    plot(ft_struct.freq(i_freq_oi), plt_spectrum,  '-k', 'LineWidth', .5, 'Marker','none');

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

    i_rej_foi   = ft_struct.freq(i_freq_oi) >  freq_min & ...
                  ft_struct.freq(i_freq_oi) <  freq_max;  
    i_sess_rej  =  any(...
                            plt_spectrum(:, i_rej_foi) > pwr_min & ...
                            plt_spectrum(:, i_rej_foi) < pwr_max...
                        , 2);

    hold on
    if any(i_sess_rej)
        plot(ft_struct.freq(i_freq_oi), plt_spectrum(i_sess_rej, :), 'LineWidth', .5, 'Marker','none');
    end

    i_ch_sess_rej(:, i_ch) = i_sess_rej;

end

exportgraphics(gcf, fullfile(save_dir, '2_art_rejected_sess_power_spectra.png'))

%%
set(0,'DefaultFigureVisible','off')
figure('Units','inches','Position',[0, 3, 18, 12]);

tiledlayout(2,2, 'Padding','tight')

sgtitle(sprintf('%s | %s–%s, Stage 1 (prior to analgesic stim testing)\n%g sessions total',...
    pt_side_id, datestr(t_start, 'dd-mmm-yyyy'), datestr(t_end, 'dd-mmm-yyyy'), sum(any(~i_ch_sess_rej,2))),...
    'FontSize', 24)


for i_ch = 1:4
    nexttile

    title(sprintf('Ch%g | %s', i_ch - 1, ft_struct.label{i_ch}));

    spectra_oi =  squeeze(ft_struct.powspctrm(:,i_ch, i_freq_oi));

    spectra_oi(i_ch_sess_rej(:, i_ch)==1, :) = nan;

    plt_spectrum =  log10(spectra_oi);

    plt_t_vec =  t_vec;
    h_sur     = surface(ft_struct.freq(i_freq_oi), plt_t_vec, plt_spectrum);
    
    format_spectra(cfg, h_sur);
    
    %% on 2nd y axis plot power spectra on top of each other (rather than in spectrogram)
    yyaxis right;         set(gca, 'YColor', 'k'); yticklabels([]);     hold on
    plot(ft_struct.freq(i_freq_oi), plt_spectrum,  '-k', 'LineWidth', .5, 'Marker','none');

    ylim(cfg.pwr_limits);


end
exportgraphics(gcf, fullfile(save_dir, '2_pruned_sess_power_spectra.png'))
set(0,'DefaultFigureVisible','off')
end