function peak_sum_tbl ...
    ...
    = plot_interstage_distance_to_spectra(...
    ...
    dirs, peak_sum_tbl,s0_s1_dist_tbl, fig_name)

% loop through pt hemi channels and assign distance table for given
% contacts inputted to channel (2 contacts per channel)
for i = 1 : length(peak_sum_tbl.Row)

    pt_oi     = peak_sum_tbl.Row{i}(1:5);
    i_pt_oi   = strcmpi(s0_s1_dist_tbl.Subject, pt_oi);

    hemi_dist = s0_s1_dist_tbl(i_pt_oi, :);

    % "labels_2" refers to RC+S formatted labels
    electrode_oi =  peak_sum_tbl{i, {'labels_2_1_TW_format', 'labels_2_2_TW_format'}};
    i_label       = contains(hemi_dist.Label, electrode_oi);

    peak_sum_tbl.avg_SlicerEuclidainMin(i) ...
        ...
        = mean(hemi_dist(i_label, :).SlicerEuclidianMin);

    peak_sum_tbl.avg_SlicerRadialError(i) ...
    ...
    = mean(hemi_dist(i_label, :).SlicerRadialError);

    % return whole TW distance table subsetted for given bipolar channel
    % (i.e., two contacts *each of which* have a distance estimate)
    peak_sum_tbl.hemi_dist_tbl{i} = hemi_dist(i_label, :);

end

%% correlate mean distance between channels to absolute difference between mean spectra
figure('Units','inches','Position',[5, 5, 7, 3.5]);

colors = [brewermap(length(peak_sum_tbl.Row)-1, 'Dark2');...
            [0,0,0]];

% loop per hemisphere, such that, each point can have it's own legend entry
for i = 1: length(peak_sum_tbl.Row)
    scatter(peak_sum_tbl.avg_SlicerEuclidainMin(i), ...
            peak_sum_tbl.abs_spec_diff_in_db(i), ...
            125, colors(i, :), 'filled','MarkerFaceAlpha',.5);
    hold on
end

% loop per hemisphere, such that, each point can have it's own legend entry
for i = 1: length(peak_sum_tbl.Row)
    scatter(peak_sum_tbl.avg_SlicerRadialError(i), ...
            peak_sum_tbl.abs_spec_diff_in_db(i), ...
            125, colors(i, :), 'filled','Marker','diamond','MarkerFaceAlpha',.5);
end

% 29-Aug-2023, RBL
% visually saw RCS05R as outlier 
%  --> remove from correlation and list as limitation
i_keep = ~strcmp(peak_sum_tbl.Row, 'RCS05R');

[rho_euc, p_euc] = corr(peak_sum_tbl.abs_spec_diff_in_db(i_keep),...
                peak_sum_tbl.avg_SlicerEuclidainMin(i_keep),...
                'type','Pearson');

[rho_rad, p_rad] = corr(peak_sum_tbl.abs_spec_diff_in_db(i_keep),...
                peak_sum_tbl.avg_SlicerRadialError(i_keep),...
                'type','Pearson');

% plotting niceties
title("Pearson's corr (w/o RCS05R)");

legend(peak_sum_tbl.Row, 'Location','northeastoutside'); legend box off

ylabel(...
    sprintf('|spectra difference| (db)'))

xlabel(...
    sprintf('Euclidean and Radial Error\nbetween bipolar channels'))

grid on

set(gca, 'FontSize', 14)

text(2, 0.03, sprintf(['euclidean | r = %.2f; p = %.0e\n',...
                       'radial    | r = %.2f; p = %.0e'],...
            rho_euc, p_euc, rho_rad, p_rad), 'FontSize',10)

save_dir = fullfile(dirs.rcs_pia, 'ephy_analysis/staged_spectra_stability/');

exportgraphics(gcf, [save_dir, fig_name, '.png']);

end