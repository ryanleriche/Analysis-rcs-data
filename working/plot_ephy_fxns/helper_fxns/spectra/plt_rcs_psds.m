function plt_session_spectra(cfg, fft_bins_inHz, pt_side_id, pwrspectra_by_sess, ...
                                    ch_names, par_db_oi, save_dir)



t_start = dateshift( par_db_oi.timeStart(1), 'start', 'day');
t_end   = dateshift( par_db_oi.timeStart(end), 'end', 'day');


t_vec   = t_start:duration('24:00:00'):t_end;


plt_pwrspectra_by_sess = nan(length(t_vec), 4, length(fft_bins_inHz));

for i_day = 1 :length(t_vec)-1

    i_ss = find(ge(par_db_oi.timeStart, t_vec(i_day)) &...
           le(par_db_oi.timeStart, t_vec(i_day+1)));


    plt_pwrspectra_by_sess(i_day, :, :) ...
        =...
        mean(pwrspectra_by_sess(i_ss, :, :), 1);

end

if cfg.z_score
    plt_lbl = {'z_every'};

else
    plt_lbl = {'every', 'mean_std'};
end

i_freq_oi     = find((fft_bins_inHz >= cfg.freq(1) & fft_bins_inHz <= cfg.freq(2)));

%%%
set(0,'DefaultFigureVisible','off')
for j_plt = 1:length(plt_lbl)
    figg = figure(j_plt);
    figg.Units    = 'Inches';
    figg.Position= [0, 3, 18, 12];


    tiledlayout(2,2, 'Padding','tight')
    
    sgtitle(sprintf('%s | %s–%s, Stage 1 (prior to analgesic stim testing)\n%g sessions total',...
        pt_side_id, datestr(t_start, 'dd-mmm-yyyy'), datestr(t_end, 'dd-mmm-yyyy'), height(par_db_oi)),...
        'FontSize', 24)
    
    for i_ch = 1:4
        nexttile
        plt_ch_name = unique(ch_names(:, i_ch));

        title(sprintf('Ch%g | %s', i_ch - 1, plt_ch_name{:}));

      
        if cfg.z_score
            zscor_xnan     = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

            plt_spectrum = zscor_xnan(...
                            log10(squeeze(plt_pwrspectra_by_sess(:, i_ch, i_freq_oi)))...
                        );

             c_str = ['z-scored',newline, 'log_{10}(mV^{2}/Hz)', newline,'[ACROSS sessions of given frequency]'];

             preset_colormap(@brewermap, "-Spectral");
             cfg.pwr_limits = [-3, 3];

        else
            plt_spectrum =  log10(squeeze(plt_pwrspectra_by_sess(:, i_ch, i_freq_oi)));
               
            c_str = "log_{10}(mV^{2}/Hz)";
            preset_colormap(@brewermap, "YlOrRd");
        end

        switch cfg.plotted_time
            case 'dates_per_se'
                 plt_t_vec =  t_vec;
            
            case 'since_inital_session'
                plt_t_vec =  0:length(t_vec)-1;     
        end
        h_sur = surface(fft_bins_inHz(i_freq_oi), plt_t_vec, plt_spectrum);
        hold on 

        ylabel(cfg.y_lbl_txt)
       
        h_sur.LineStyle = 'none';
    
        colormap(preset_colormap);      
        cb_hdl = colorbar;          cb_hdl.Limits = cfg.pwr_limits; 
        
        set(gca, 'FontSize', 14);
        cb_hdl.Label.String         = c_str;
        %cb_hdl.Label.Rotation       = 270;      %cb_hdl.Label.Position(1)    = 3.3;
        cb_hdl.FontSize             = 10;         cb_hdl.Label.FontSize       = 18;

        figg.CurrentAxes.YDir = 'Reverse';
        
        yyaxis right

        set(gca, 'YColor', 'k'); yticklabels([]);

        xlabel('Frequency (Hz)')
    
        if ~cfg.z_score
            if j_plt == 1
                hold on
                plot(fft_bins_inHz(i_freq_oi), plt_spectrum,  '-k', 'LineWidth', .5, 'Marker','none');
                
            else
                [lineout, ~] = stdshade(plt_spectrum, fft_bins_inHz(i_freq_oi), 0.2, 'k');
                
                lineout.LineStyle = '-';              lineout.LineWidth = 4;
            end
        end

        ylim(cfg.pwr_limits);
        if ~isfolder(save_dir);    mkdir(save_dir);  end
 
    end
    exportgraphics(gcf, [save_dir, '1_', plt_lbl{j_plt}, '_session_power_spectra.png'])
end
end
