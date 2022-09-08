function RCSXX_sum_stats = calc_sum_stats(cfg, RCSXX, db_beh_RCSXX )

    [~, RCSXX, ~] = date_parser(cfg, RCSXX, db_beh_RCSXX );

    % Simple sub-function to rearrange RCSXX table (from REDcap),
    % and calculate summary statistics (option for CI).
    stats_to_calc = {'min','max','range','mean','std'};
    sum_stats      = grpstats(RCSXX(:,2:end), [], stats_to_calc);
    
    temp_table     = rows2vars(splitvars(sum_stats(1, 2:end)));
    
    RCSXX_sum_stats = ...
        array2table(...
        reshape(temp_table.All, [length(stats_to_calc), 21]), ...
        'VariableNames', RCSXX.Properties.VariableNames(2:end),...
        'RowNames',stats_to_calc);

end