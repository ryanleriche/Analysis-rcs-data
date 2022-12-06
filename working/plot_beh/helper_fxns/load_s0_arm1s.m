function  s0_redcaps   = load_s0_arm1s(s0_dir, visits_tbl)


pts         = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};

s0_rcaps    = dir([s0_dir, '*.csv']);

for i = 1 : length(s0_rcaps)

    if any(contains(s0_rcaps(i).name, pts))

        pt_id     = pts{cellfun(@(x) contains(s0_rcaps(i).name, x), pts)};

        arm1_tbl  = readtable([s0_dir, s0_rcaps(i).name]);

        arm1_tbl  = arm1_tbl(contains(arm1_tbl.redcap_event_name, 'arm_1'),:);
        arm1_tbl  = removevars(arm1_tbl, {'redcap_event_name', 'redcap_repeat_instance'});
        


        varnames  = arm1_tbl.Properties.VariableNames;

        for j = 1 : length(varnames)
            if isdatetime(arm1_tbl.(varnames{j}))

                arm1_tbl.(varnames{j}).TimeZone = 'America/Los_Angeles';

                % remove columns with every row as NaT
                if all(isnat(arm1_tbl.(varnames{j})))

                    arm1_tbl = removevars(arm1_tbl, varnames{j});

                % fix datetime formatting
                elseif any(le(arm1_tbl.(varnames{j}), datetime('01-Jan-2019','TimeZone','America/Los_Angeles')))

                    % if every but NOT all datetimes are shifted -> toss error
                    if all(le(arm1_tbl.(varnames{j}), datetime('01-Jan-2019','TimeZone','America/Los_Angeles')))

                        error('RBL: inconsistent datetime formatting w/n %s .csv file', pt_id)
                    end
             
                    arm1_tbl.(varnames{j}) = arm1_tbl.(varnames{j}) + calyears(2000);
                    
                end
            end
        end

        i_s0_start = contains(visits_tbl.(pt_id).desc, 's0_implant');
        i_s0_stop  = contains(visits_tbl.(pt_id).desc, 's0_explant');


        % from the entirity of implant day to the entirity of the explant day
        i_s0       = ge(arm1_tbl.stim_onoff_timestamp, visits_tbl.(pt_id).dates(i_s0_start)) &...
                     le(arm1_tbl.stim_onoff_timestamp, visits_tbl.(pt_id).dates(i_s0_stop) + duration('24:00:00'));

        arm1_tbl   = arm1_tbl(i_s0, :);

       
        % remove more NaN columns (NOT redundant code)
        varnames  = arm1_tbl.Properties.VariableNames;
        arm1_vars = {'stim_onoff_timestamp', ...
                   'nrs_s0', 'intensity_vas_s0', 'mood_vas_s0', ...
                    'unpleasantness_vas_s0', 'relief_s0',...
                    'throbbing_s0', 'shooting_s0', 'stabbing_s0',...
                    'sharp_s0', 'cramping_s0', 'gnawing_s0','hot_burning_s0',...
                    'aching_s0', 'heavy_s0', 'tender_s0', 'splitting_s0',...
                    'tiring_exhausting_s0', 'sickening_s0', 'fearful_s0',...
                    'punishing_cruel_s0'};

        for j = 1 : length(varnames)
            if ~contains(varnames{j}, arm1_vars)
                if (isnumeric(arm1_tbl.(varnames{j})) && all(isnan(arm1_tbl.(varnames{j}))))...
                    ||...
                    (isdatetime(arm1_tbl.(varnames{j})) && all(isnat(arm1_tbl.(varnames{j}))))
    
                        arm1_tbl = removevars(arm1_tbl, varnames{j});
                end
            end
        end

        % rename to match Stages 1, 2, and 3 surveys
        arm1_tbl = renamevars(arm1_tbl, arm1_vars,...
                   ...
                   {'time',...
                    'mayoNRS','painVAS', 'moodVAS',...
                    'unpleasantVAS', 'reliefVAS', ...
                    'MPQthrobbing', 'MPQshooting', 'MPQstabbing',...
                    'MPQsharp', 'MPQcramping', 'MPQgnawing', 'MPQhot_burning',...
                    'MPQaching', 'MPQheavy',   'MPQtender',  'MPQsplitting', ...
                    'MPQtiring', 'MPQsickening', 'MPQfearful',...
                    'MPQcruel'});

        varnames  = arm1_tbl.Properties.VariableNames;

        i_MPQ     = find(cellfun(@(x) contains(x, 'MPQ'), varnames));

        arm1_tbl.MPQtotal = sum(arm1_tbl(:, i_MPQ).Variables, 2, 'omitnan');

        % first 11 are somatic subscore 
        arm1_tbl.MPQsom = sum(arm1_tbl(:, i_MPQ(1:11)).Variables, 2, 'omitnan');

        % last 4 are affective subscore
        arm1_tbl.MPQaff = sum(arm1_tbl(:, i_MPQ(12:15)).Variables, 2, 'omitnan');
                   
        if any(arm1_tbl.MPQtotal ~= arm1_tbl.MPQaff + arm1_tbl.MPQsom)
            error("%s | MPQ affective or somatic subscores are calculated incorrectly (RBL message).", pt_id)
        end

        % add in worst VAS and NRS as NaN to match Stages 1, 2, and 3
        if ~contains(varnames, 'worstVAS')
            arm1_tbl.worstVAS = nan(height(arm1_tbl), 1);
        end


        if ~contains(varnames, 'worstNRS')
            arm1_tbl.worstNRS = nan(height(arm1_tbl), 1);
        end

        % test surveys are sprinkled throughout Stage 0 
        % (denoted as differnt number than pt id)
        if any(contains(varnames, 'pt'))
            arm1_tbl   = arm1_tbl(arm1_tbl.pt == str2double(pt_id(end)), :);
        end


        switch pt_id 
            case 'RCS07'

                times = ...
                    datetime({'2022-09-29 14:55:07','2022-09-25 15:54:37', '2022-09-26 10:35:45','2022-09-28 14:17:25'},...
                    "TimeZone","America/Los_Angeles");

                mismatch_VAS = any(arm1_tbl.time == times, 2);
                 
                arm1_tbl     = arm1_tbl(~mismatch_VAS, :);


            otherwise
        end

        s0_redcaps.(['stage0', pt_id]) = arm1_tbl;
    end
end






end