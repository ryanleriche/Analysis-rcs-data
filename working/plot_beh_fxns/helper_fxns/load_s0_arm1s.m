function  s0_redcaps   = load_s0_arm1s(s0_dir, visits_tbl)

% visits_tbl      = visits;
% s0_dir          = [dropbox_dir, ...
%                     '/DATA ANALYSIS/Stage 0 ALL PATIENTS Redcap Records/'];


pts         = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};

s0_rcaps    = dir([s0_dir, '*.csv']);

for i = 1 : length(s0_rcaps)

    if any(contains(s0_rcaps(i).name, pts))

        pt_id     = pts{cellfun(@(x) contains(s0_rcaps(i).name, x), pts)};

        arm1_tbl  = readtable([s0_dir, s0_rcaps(i).name]);

        arm1_tbl  = arm1_tbl(contains(arm1_tbl.redcap_event_name, 'arm_1'),:);

        % go through and properly format datetimes if need be
        varnames  = arm1_tbl.Properties.VariableNames;

        for j = 1 : length(varnames)
            if isdatetime(arm1_tbl.(varnames{j}))

                arm1_tbl.(varnames{j}).TimeZone = 'America/Los_Angeles';

                % remove columns with every row as NaT
                if all(isnat(arm1_tbl.(varnames{j})))

                    arm1_tbl = removevars(arm1_tbl, varnames{j});

                % fix datetime formatting
                elseif any(le(arm1_tbl.(varnames{j}), datetime('01-Jan-2019','TimeZone','America/Los_Angeles')))

%                     % if some but NOT all datetimes are shifted -> toss error
%                     if ~all(le(arm1_tbl.(varnames{j}), datetime('01-Jan-2019','TimeZone','America/Los_Angeles')))
% 
%                         error('RBL: inconsistent datetime formatting w/n %s .csv file', pt_id)
%                     end
             
                    arm1_tbl.(varnames{j}) = arm1_tbl.(varnames{j}) + calyears(2000);
                    
                end
            end
        end


        % only return surveys that occured DURING pt's Stage 0
        i_s0_start = contains(visits_tbl.(pt_id).desc, 's0_implant');
        i_s0_stop  = contains(visits_tbl.(pt_id).desc, 's0_explant');


        % from the entirety of implant day to the entirety of the explant day
        i_s0       = ge(arm1_tbl.stim_onoff_timestamp, visits_tbl.(pt_id).dates(i_s0_start)) &...
                     le(arm1_tbl.stim_onoff_timestamp, visits_tbl.(pt_id).dates(i_s0_stop) + duration('24:00:00'));

        arm1_tbl   = arm1_tbl(i_s0, :);

        % remove columns that only contain NaNs or NaTs
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
            if all(~strcmp(varnames(j), arm1_vars))
                if (isnumeric(arm1_tbl.(varnames{j})) && all(isnan(arm1_tbl.(varnames{j}))))...
                    ||...
                    (isdatetime(arm1_tbl.(varnames{j})) && all(isnat(arm1_tbl.(varnames{j}))))
    
                        arm1_tbl = removevars(arm1_tbl, varnames{j});
                end
            end
        end

        % arm1_tbl  = removevars(arm1_tbl, {'redcap_event_name', 'redcap_repeat_instance'});


        % rename to match Stages 1, 2, and 3 surveys
        arm1_tbl = renamevars(arm1_tbl, arm1_vars,...
                   ...
                   {'time',...
                    'mayoNRS', 'painVAS', 'moodVAS',...
                    'unpleasantVAS', 'reliefVAS', ...
                    'MPQthrobbing', 'MPQshooting', 'MPQstabbing',...
                    'MPQsharp', 'MPQcramping', 'MPQgnawing', 'MPQhot_burning',...
                    'MPQaching', 'MPQheavy',   'MPQtender',  'MPQsplitting', ...
                    'MPQtiring', 'MPQsickening', 'MPQfearful',...
                    'MPQcruel'});

        varnames  = arm1_tbl.Properties.VariableNames;

        i_MPQ     = find(cellfun(@(x) contains(x, 'MPQ'), varnames));

        % first 11 are somatic subscore 
        arm1_tbl.MPQsom = sum(arm1_tbl(:, i_MPQ(1:11)).Variables, 2, 'omitnan');

        % rather than zero, if all entries are NaN, then subscore if NaN
        i_som_nan       = all(isnan(arm1_tbl(:, i_MPQ(1:11)).Variables),2);
        arm1_tbl.MPQsom(i_som_nan) = NaN;

        % last 4 are affective subscore
        arm1_tbl.MPQaff            = sum(arm1_tbl(:, i_MPQ(12:15)).Variables, 2, 'omitnan');
        i_aff_nan                  = all(isnan(arm1_tbl(:, i_MPQ(12:15)).Variables),2);
        arm1_tbl.MPQaff(i_aff_nan) = NaN;


        arm1_tbl.MPQtotal = sum(arm1_tbl(:, i_MPQ).Variables, 2, 'omitnan');
        arm1_tbl.MPQtotal(i_som_nan & i_aff_nan) = NaN;

        
                   
        if any(arm1_tbl.MPQtotal(~(i_som_nan | i_aff_nan)) ~= ...
                arm1_tbl.MPQaff(~(i_som_nan | i_aff_nan)) + ...
                arm1_tbl.MPQsom(~(i_som_nan | i_aff_nan)))
            error("%s | MPQ affective or somatic subscores are calculated incorrectly (RBL message).", pt_id)
        end

        arm1_tbl = movevars(arm1_tbl, {'MPQtotal', 'MPQsom', 'MPQaff'}, ...
                            "Before", varnames(i_MPQ(1)));

        % rename pain metrics NOT consistent across pts
        if any(contains(varnames, 'best_s0'))
            arm1_tbl = renamevars(arm1_tbl, varnames(contains(varnames, 'best_s0')),...
                'bestNRS');
        end

        if any(contains(varnames, 'worst_s0'))
            arm1_tbl = renamevars(arm1_tbl, varnames(contains(varnames, 'worst_s0'))...
                ,'worstNRS');
        end

        if any(contains(varnames, 'nrs_left_arm'))
            arm1_tbl = renamevars(arm1_tbl, varnames(contains(varnames, 'nrs_left_arm'))...
                ,'leftarmNRS');
        end

        if any(contains(varnames, 'nrs_left_leg'))
            arm1_tbl = renamevars(arm1_tbl, varnames(contains(varnames, 'nrs_left_leg'))...
                ,'leftlegNRS');
        end


        if any(contains(varnames, 'nrs_left_face'))
            arm1_tbl = renamevars(arm1_tbl, varnames(contains(varnames, 'nrs_left_face'))...
                ,'leftfaceNRS');
        end


        % test surveys are sprinkled throughout Stage 0 
        % (denoted as differnt number than pt id)
        if any(contains(varnames, 'pt'))
            arm1_tbl   = arm1_tbl(arm1_tbl.pt == str2double(pt_id(end)), :);
        end

% 
%         switch pt_id 
%             case 'RCS07'
% 
% %                 times = ...
% %                     datetime({'2022-09-29 14:55:07','2022-09-25 15:54:37', '2022-09-26 10:35:45','2022-09-28 14:17:25'},...
% %                     "TimeZone","America/Los_Angeles");
% % 
% %                 mismatch_VAS = any(arm1_tbl.time == times, 2);
% %                  
% %                 arm1_tbl     = arm1_tbl(~mismatch_VAS, :);
% % 
% 
%             otherwise
%         end


        var_oi = arm1_tbl.Properties.VariableNames;
        var_oi = var_oi(contains(var_oi, {'NRS', 'VAS', 'MPQtotal'}));

        arm1_tbl = arm1_tbl(~all(isnan(arm1_tbl{:,var_oi}),2), :);
         %arm1_tbl =  arm1_tbl(all(isnan( =arm1_tbl)));
        s0_redcaps.(['stage0', pt_id]) = arm1_tbl;
    end
end

end