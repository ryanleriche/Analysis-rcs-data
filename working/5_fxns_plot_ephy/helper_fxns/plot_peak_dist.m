function plot_peak_dist(dirs, fooof_struct_1, fooof_struct_2, fig_name)

set(0,'DefaultFigureVisible','on')
close all

figure('Units','inches','Position',[0, 0, 10, 10])
tiledlayout(3,3, 'TileSpacing','loose');
nexttile
sgtitle('FOOOFed frequency peaks per session between stages'); hold on

for i_chan = 1 : height(fooof_struct_1.sense_chan)

    tmp_fooof = fooof_struct_1.sense_chan.fooof_params{i_chan};

    freq_order_1 = nan(height(tmp_fooof), 4);

    for i_sess = 1 : height(tmp_fooof)

       [~, i_amp] = sort(tmp_fooof.peak_power{i_sess});

       freq_order_1(i_sess, 1:length(i_amp)) ...
           = ...
           tmp_fooof.center_freq{i_sess}(i_amp);

    end

    tmp_fooof    = fooof_struct_2.sense_chan.fooof_params{i_chan};
    freq_order_2 = nan(height(tmp_fooof), 4);

    for i_sess = 1 : height(tmp_fooof)

       [~, i_amp] = sort(tmp_fooof.peak_power{i_sess});

       freq_order_2(i_sess, 1:length(i_amp)) ...
           = ...
           tmp_fooof.center_freq{i_sess}(i_amp);

    end

    peak_freqs_1  = freq_order_1(:, 1);
    peak_freqs_2  = freq_order_2(:, 1);

    if i_chan ~= 1; nexttile; end
    % Specify bin edges or centers
    min_freq = min([peak_freqs_1; peak_freqs_2]);
    max_freq = max([peak_freqs_1; peak_freqs_2]);

    min_freq = 2;
    max_freq = 40;

    binEdges   = linspace(min_freq, max_freq, 10);
    % Plot histograms with equal bin sizes
  
    histogram(peak_freqs_1, binEdges, 'Normalization', 'probability', 'FaceColor', 'blue', 'EdgeColor', 'none');
    hold on;
    histogram(peak_freqs_2, binEdges, 'Normalization', 'probability', 'FaceColor', 'red', 'EdgeColor', 'none');


    chan_1  = strsplit(fooof_struct_1.sense_chan.Properties.RowNames{i_chan},'_');
    chan_2  = strsplit(fooof_struct_2.sense_chan.Properties.RowNames{i_chan},'_');


    legend(sprintf('%s (inpatient)', chan_1{2}),...
           sprintf('%s (ambulatory)', chan_2{2}),...
           'Interpreter', 'none');

    %title(sprintf('%s; BF_{01} = %.1f', chan_2{1}, peak_freq_tbl.bf01(i_chan)));

    ylim([0,1]); xlabel('Frequency (Hz)');

    ylabel('Probability')

    %mean_freq = mean([min_freq,max_freq]);

    xlim([min_freq, max_freq])
    
    legend box off
    set(gca, 'FontSize', 10)

end

save_dir = fullfile(dirs.rcs_pia, 'ephy_analysis/staged_spectra_stability/');

exportgraphics(gcf, [save_dir, fig_name, '.png']);


end

    
