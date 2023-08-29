function plot_primary_peaks_2(dirs, peak_freq_tbl, fooof_struct_1, fooof_struct_2, fig_name)

figure('Units','inches','Position',[0, 0, 12, 10]);

tiledlayout(3,3, 'TileSpacing','loose');
nexttile

colors = brewermap([],"Paired");

leg_plted = false;

for i_chan = 1 :  height(fooof_struct_1.sense_chan)
    if i_chan ~= 1; nexttile; end


    min_freq = min(fooof_struct_1.fft_bins_inHz);
    max_freq = max(fooof_struct_1.fft_bins_inHz);

    res_freq = max(diff(fooof_struct_1.fft_bins_inHz));
    
    N_peaks  = length(peak_freq_tbl.cent_freq_1{i_chan});

    cnt = 0;

    for i_peak = 1 : N_peaks

        peak_freqs_1  = peak_freq_tbl.cent_freq_1{i_chan}{i_peak};
        peak_freqs_2  = peak_freq_tbl.cent_freq_2{i_chan}{i_peak};
   
        binEdges   = min_freq: res_freq: max_freq;
        % Plot histograms with equal bin sizes
      
        histogram(peak_freqs_1, binEdges,...
            'Normalization', 'probability',...
            'FaceColor', colors(i_peak+cnt, :), 'EdgeColor', 'none');
        hold on;

        histogram(peak_freqs_2, binEdges,...
            'Normalization', 'probability', ...
            'FaceColor', colors(i_peak+cnt+1,:), 'EdgeColor', 'none');
    
        cnt = cnt+1;

        bf = peak_freq_tbl.bf10{i_chan}(i_peak);

        if bf > 1/3
           bf_form =  'k-.';

        else
            bf_form =  'r-.';

        end
        xline(median(peak_freqs_1,"all"), ...
            bf_form,...
            sprintf('BF_{10} = %.0e ', bf),...
            'HandleVisibility', 'off');
       
    end

     if N_peaks ==3 && ~leg_plted

       leg_txt = {'largest inpat. peak',  'largest amb. peak',...
               '2nd largest inpat. peak', '2nd largest amb. peak',...
               '3rd largest inpat. peak', '3rd largest amb. peak'};

        leg = legend(leg_txt, 'Location','northoutside','NumColumns',3, 'FontSize', 12);
        leg.Position(1) = 0.245;
        leg.Position(2) = 0.945;
        leg_plted = true;


     end


    chan_1  = strsplit(fooof_struct_1.sense_chan.Properties.RowNames{i_chan},'_');
    chan_2  = strsplit(fooof_struct_2.sense_chan.Properties.RowNames{i_chan},'_');

    if i_chan < 4
    title(...
        sprintf('\n\n\n%s', chan_2{1}),...
        sprintf('%s (inpat) -> %s (amb)', chan_1{2}, chan_2{2}),...
        ...
        'FontSize', 10)

    else
     title(...
        sprintf('%s', chan_2{1}),...
        sprintf('%s (inpat) -> %s (amb)', chan_1{2}, chan_2{2}),...
        ...
        'FontSize', 10)
    end

  
    ylim([0,.5]);

    xlim([3, 30])

    ylabel('Probability', 'FontSize', 10)

    xlabel('Frequency (Hz)', 'FontSize', 10);
    
end

save_dir = fullfile(dirs.rcs_pia, 'ephy_analysis/staged_spectra_stability/');

exportgraphics(gcf, [save_dir, fig_name, '.png']);

end

    
