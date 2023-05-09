function rcap_stats = calc_sum_stats(cfg, redcap, varargin)

    [redcap, ~] = date_parser(cfg, redcap);

    % Simple sub-function to rearrange RCSXX table (from REDcap),
    % and calculate summary statistics (option for CI).
    stats_to_calc = {'min','max','range','mean','std'};

    i_vars         = find(contains(redcap.Properties.VariableNames,...
                          {'NRS', 'VAS', 'MPQ', 'PBD'}));

    sum_stats      = grpstats(redcap(:,i_vars), [], stats_to_calc);
    
    temp_table     = rows2vars(splitvars(sum_stats(1, 2:end)));
    
    rcap_stats = ...
        array2table(...
        reshape(temp_table.All, [length(stats_to_calc), length(i_vars)]), ...
        'VariableNames', redcap.Properties.VariableNames(i_vars),...
        'RowNames', stats_to_calc);

    % add in number of non-NaN surveys per variable and N surveys total
    for i=1:length(i_vars)
        rcap_stats{'N_variable', i}         = sum(~isnan(redcap{:,i_vars(i)}));

        % improvement would be an increase in releif and mood metrics
        if contains(redcap.Properties.VariableNames(i_vars(i)), {'relief', 'mood'})

            rcap_stats{'half_improve', i}     = rcap_stats{'mean', i} + rcap_stats{'mean', i} *.5;   
            rcap_stats{'third_improve', i}    = rcap_stats{'mean', i} + rcap_stats{'mean', i} *(1/3);  
        else
            rcap_stats{'half_improve', i}     = rcap_stats{'mean', i} - rcap_stats{'mean', i} *.5;   
            rcap_stats{'third_improve', i}    = rcap_stats{'mean', i} - rcap_stats{'mean', i} *(1/3);  

        end
    end

    varnames = rcap_stats.Properties.VariableNames;

    for i = 1 : size(rcap_stats,2)

        if iscell(rcap_stats.(varnames{i}))

            if isnumeric(rcap_stats.(varnames{i}){1})

                rcap_stats.(varnames{i}) = cellfun(@(x) x, rcap_stats.(varnames{i}));
            end
        end
    end
end