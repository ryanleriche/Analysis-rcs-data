function plt_rcs_PSDs(cfg, pt_side_id, save_dir, ft_freq_struct_1, ft_freq_struct_2)



if cfg.z_score;    plt_lbl = {'z_every'};
else;              plt_lbl = {'every', 'mean_std'};
end

i_freq_oi     = find((ft_freq_struct_1.freq >= cfg.freq(1) & ft_freq_struct_1.freq <= cfg.freq(2)));
freq_oi       = ft_freq_struct_1.freq(i_freq_oi);
%%
close all
set(0,'DefaultFigureVisible','off')

for j_plt = 1  :  length(plt_lbl)
    figg = figure(j_plt);
    figg.Units    = 'Inches';
    figg.Position= [0, 3, 10, 8];

    tiledlayout(2,2, 'Padding','tight')
    
    sgtitle(sprintf('%s | %g spectra | %s \n', pt_side_id, height(ft_freq_struct_1.rcs.par_db), cfg.title),...
                'FontSize', 24, 'Interpreter', 'none')
    
    for i_ch = 1:4
        nexttile
        plt_ch_name = unique(ft_freq_struct_1.label(i_ch));

        title(sprintf('Ch%g | %s', i_ch - 1, plt_ch_name{:}));

      
        if cfg.z_score
            zscor_xnan      = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
            plt_spectrum_1  = zscor_xnan(log10(squeeze(ft_freq_struct_1.powspctrm(:, i_ch, i_freq_oi))));
            plt_spectrum_2  = zscor_xnan(log10(squeeze(ft_freq_struct_2.powspctrm(:, i_ch, i_freq_oi))));

             cfg.pwr_limits = [-4, 4];
        else
            plt_spectrum_1 =  log10(squeeze(ft_freq_struct_1.powspctrm(:, i_ch, i_freq_oi)));
            plt_spectrum_2 =  log10(squeeze(ft_freq_struct_2.powspctrm(:, i_ch, i_freq_oi)));

        end

        if j_plt == 1
            hold on

            for h = 1 : size(plt_spectrum_1, 1)
                hand_plt = plot(freq_oi, plt_spectrum_1(h,:),  '-k', 'LineWidth', .5, 'Marker','none');

                hand_plt.Color(4) = .1;
            end

            for h = 1 : size(plt_spectrum_1, 1)
                hand_plt = plot(freq_oi, plt_spectrum_2(h,:),  '-b', 'LineWidth', .5, 'Marker','none');

                hand_plt.Color(4) = .1;
            end
            
        else
            [lineout, ~] = stdshade(plt_spectrum_1, freq_oi, 0.2, 'k'); hold on
            
            lineout.LineStyle = '-';              lineout.LineWidth = 2;

            [lineout, ~] = stdshade(plt_spectrum_2, freq_oi, 0.2, 'b');
            
            lineout.LineStyle = '-';              lineout.LineWidth = 2;

            legend(cfg.txt_leg)
        
        end

        title(ft_freq_struct_1.label(i_ch));

        ylim(cfg.pwr_limits);     ylabel("log_{10}(mV^{2}/Hz)")
        xlabel('Frequency (Hz)')

        if ~isfolder(save_dir);    mkdir(save_dir);  end
 
    end
    exportgraphics(gcf, [save_dir, cfg.title, '_', plt_lbl{j_plt}, '_session_power_spectra.png'])
end
end
