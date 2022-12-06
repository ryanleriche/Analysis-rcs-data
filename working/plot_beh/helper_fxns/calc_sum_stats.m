function rcap_stats = calc_sum_stats(cfg, redcap, varargin)

    [redcap, ~] = date_parser(cfg, redcap);

    % Simple sub-function to rearrange RCSXX table (from REDcap),
    % and calculate summary statistics (option for CI).
    stats_to_calc = {'min','max','range','mean','std'};
    sum_stats      = grpstats(redcap(:,2:end), [], stats_to_calc);
    
    temp_table     = rows2vars(splitvars(sum_stats(1, 2:end)));
    
    rcap_stats = ...
        array2table(...
        reshape(temp_table.All, [length(stats_to_calc), (size(redcap,2) -1)]), ...
        'VariableNames', redcap.Properties.VariableNames(2:end),...
        'RowNames',stats_to_calc);

    varnames = rcap_stats.Properties.VariableNames;

    for i = 1 : size(rcap_stats,2)

        if iscell(rcap_stats.(varnames{i}))

            if isnumeric(rcap_stats.(varnames{i}){1})

                rcap_stats.(varnames{i}) = cellfun(@(x) x, rcap_stats.(varnames{i}));
            end
        end
    end
end