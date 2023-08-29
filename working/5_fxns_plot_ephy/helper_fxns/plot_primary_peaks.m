function plot_primary_peaks(dirs, peak_freq_tbl, fooof_struct_1, fooof_struct_2, fig_name)


for i_peak = 1:3


    figure('Units','inches','Position',[0, 0, 10, 10])
    tiledlayout(3,3, 'TileSpacing','loose');
    nexttile

    sgtitle(...
        sprintf('%s largest peak of inpatient stage (across sessions)', iptnum2ordinal(i_peak))...
        ); hold on

    
    for i_chan = 1 : height(fooof_struct_1.sense_chan)
        
        if i_chan ~= 1; nexttile; end
        try
            peak_freqs_1  = peak_freq_tbl.cent_freq_1{i_chan}{i_peak};
            peak_freqs_2  = peak_freq_tbl.cent_freq_2{i_chan}{i_peak};
        
            
            % Specify bin edges or centers
            min_freq = min([peak_freqs_1; peak_freqs_2]);
            max_freq = max([peak_freqs_1; peak_freqs_2]);
        
        
            binEdges   = linspace(min_freq, max_freq, 15);
            % Plot histograms with equal bin sizes
          
            histogram(peak_freqs_1, binEdges, 'Normalization', 'probability', 'FaceColor', 'blue', 'EdgeColor', 'none');
            hold on;
            histogram(peak_freqs_2, binEdges, 'Normalization', 'probability', 'FaceColor', 'red', 'EdgeColor', 'none');
        
        
            chan_1  = strsplit(fooof_struct_1.sense_chan.Properties.RowNames{i_chan},'_');
            chan_2  = strsplit(fooof_struct_2.sense_chan.Properties.RowNames{i_chan},'_');
        
        
            legend(sprintf('%s (inpatient)', chan_1{2}),...
                   sprintf('%s (ambulatory)', chan_2{2}),...
                   'Interpreter', 'none');

            legend box off
        
            title(sprintf('%s; BF_{01} = %.1f', chan_2{1}, peak_freq_tbl.bf01{i_chan}(i_peak)));

            mean_freq = mean([min_freq,max_freq]);
            xlim([mean_freq-5, mean_freq+5])
        catch
        end
        ylim([0,0.5]); xlabel('Frequency (Hz)');
    
        ylabel('Probability')
    
       
        set(gca, 'FontSize', 10)
    
    end
    
    save_dir = fullfile(dirs.rcs_pia, 'ephy_analysis/staged_spectra_stability/');
    
    save_name = sprintf('%s (%s largest)', fig_name, iptnum2ordinal(i_peak));
    exportgraphics(gcf, [save_dir, save_name, '.png']);

end
end

    
