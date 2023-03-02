function [stimLog, redcap] ...
    ...
    = align_REDcap_to_stimLog(...
    ...
    cfg, db_RCSXXX, redcap)

% cfg                    = [];
% cfg.load_EventLog      = false;
% cfg.ignoreold          = false;
% cfg.raw_dir            = pia_raw_dir;
% cfg.stage_dates                 = stage_dates{str2double(pt_sides{i}(end-1))};
% cfg.pt_id                       = pt_sides{i}(1:end-1);
% 
% db_RCSXXX = db.(pt_sides{i});
% redcap    = REDcap.(pt_sides{i}(1:end-1));
% 

i_stimlog    = find(cellfun(@(x) ~isempty(x), db_RCSXXX.stimLogSettings));
i_devset     = find(cellfun(@(x) ~isempty(x), db_RCSXXX.stimSettingsOut));

i_metadata   = find(cellfun(@(x) ~isempty(x), db_RCSXXX.metaData));

for i = 1 : height(db_RCSXXX.stimLogSettings)
    if any(i == i_stimlog)
         % easy case, use metaData UTC offset from the SAME streaming session
        if any(i == i_metadata)
            j = i;
        else
            % per non-empty StimLog, find the nearest non-empty metaData log 
            t_diff        = i - i_metadata;
            near_t        = min(t_diff(t_diff >= 0));
    
            if isempty(near_t)
                near_t     = max(t_diff(t_diff < 0));
            end
            j             = i_metadata(t_diff == near_t);   
        end

     UTCoffset     = sprintf('%+03.0f:00', db_RCSXXX.metaData{j}.UTCoffset);
     HostUnixTime  = db_RCSXXX.stimLogSettings{i}.HostUnixTime; 


     db_RCSXXX.stimLogSettings{i}.time_stimLog = ...
            datetime(HostUnixTime /1000,'ConvertFrom','posixTime','TimeZone', ...
            UTCoffset,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
    end


    if any(i == i_devset)
         % easy case, use metaData UTC offset from the SAME streaming session
        if any(i == i_metadata)
            j = i;
        else
            % per non-empty StimLog, find the nearest non-empty metaData log 
            t_diff        = i - i_metadata;
            near_t        = min(t_diff(t_diff >= 0));
    
            if isempty(near_t)
                near_t     = max(t_diff(t_diff < 0));
            end
            j             = i_metadata(t_diff == near_t);   
        end
        
        % also add datetime to 'stimSettingsOut' which are the final
        % settings from a streaming session
        UTCoffset     = sprintf('%+03.0f:00', db_RCSXXX.metaData{j}.UTCoffset);
        HostUnixTime  = db_RCSXXX.stimSettingsOut{i}.HostUnixTime; 
    
        db_RCSXXX.stimSettingsOut{i}.time_devset = ...
            datetime(HostUnixTime /1000,'ConvertFrom','posixTime','TimeZone', ...
            UTCoffset,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
    end
end

% takes stimLogSettings to a single sorted timetable and removes redundant fields
stimLog =   movevars(...
                     vertcat(db_RCSXXX.stimLogSettings{i_stimlog}),...
            'time_stimLog', 'Before', 'HostUnixTime');
   



vars    = {'time_stimLog', 'activeGroup', 'therapyStatus', 'therapyStatusDescription',...
           'stimParams_prog1', 'stimParams_prog2', 'stimParams_prog3', 'stimParams_prog4'};


[~, i_u]  = unique(stimLog(:, vars), 'rows');
stimLog   = stimLog(i_u, :);


% As of Nov. 2022, we have not intentionally used multiple programs w/n a
% group, display all the unique stim settings of Programs 2, 3, and 4
disp(strjoin(...
        [cfg.pt_id;{'| Unique Program 2 Settings ->'}; unique(stimLog.stimParams_prog2)]...
            ));

disp(strjoin(...
        [cfg.pt_id;{'| Unique Program 3 Settings ->'}; unique(stimLog.stimParams_prog3)]...
            ));

disp(strjoin(...
        [cfg.pt_id;{'| Unique Program 4 Settings ->'}; unique(stimLog.stimParams_prog4)]...
            ));


% remove redundant fields
stimLog = removevars(stimLog, ...
                      {'therapyStatus','HostUnixTime','updatedParameters',...
                      'stimParams_prog2', 'stimParams_prog3', 'stimParams_prog4'});
%% assign nearest REDcap to StimLog and vice versa
% use comprehensive stim parameters (from DeviceSettings.json) to add in cycling, ramping, and Active Recharging

stimSettingsOut  = vertcat(db_RCSXXX.stimSettingsOut{i_devset});


i_empty          = cellfun(@isempty, [stimSettingsOut.GroupA]);

stimSettingsOut  = stimSettingsOut(~i_empty, :);
% use of Program 1 is fine assumption given lack of meaningful stim params
% in Programs 2, 3, and 4

for i = 1 : height(stimLog)

    stim_para = strsplit(char(stimLog.stimParams_prog1(i)),',');
    
    if ~strcmp(stim_para, 'Disabled')

        if isempty(stim_para{1})
            stimLog.stimContacts{i}  = '';
        else
            stimLog.stimContacts{i}  = stim_para{1};
        end
    
        stimLog.ampInMilliamps(i)              = str2double(stim_para{2}(2:end -2));
        stimLog.pulseWidthInMicroseconds(i)    = str2double(stim_para{3}(2:end -2));
        stimLog.rateInHz(i)                    = str2double(stim_para{4}(2:end -2));
    
        % from nearest, previous 'stimSettingsOut' 
    
        t_diff        = stimLog.time_stimLog(i) - stimSettingsOut.time_devset;
        near_t        = min(t_diff(t_diff >= 0));
    
        j             = find(t_diff == near_t);
        abcd          = stimLog.activeGroup{i};
    
        ss_group      = stimSettingsOut.(['Group' abcd]){j};
    
    
        stimLog.cycleEnabled(i)    = ss_group.cyclingEnabled;
        stimLog.cycleOnInSecs(i)   = ss_group.cycleOnInSecs;
        stimLog.cycleOffInSecs(i)  = ss_group.cycleOffInSecs;
        
        stimLog.rampInSecs(i)      = ss_group.rampInSecs;
        stimLog.rampRepeat(i)      = ss_group.rampRepeat;
    
        stimLog.achRechRatio(i)    = ss_group.actRechRatio(1);
    
    
        % find REDcap setting that occurs AFTER given stimLog entry, but BEFORE
        % the subsequent entry -> pull that REDcap survey
        if i < height(stimLog)
            btwn_this_and_next = ...
                find(ge(redcap.time, stimLog.time_stimLog(i)) & ...
                le(redcap.time, stimLog.time_stimLog(i + 1)));
        
            if ~isempty(btwn_this_and_next)
            
                stimLog.redcap_btwn_stimLogs{i}     = 'btwn_stimLogs';
                stimLog.redcap_reports{i}           = redcap(btwn_this_and_next,:);
                stimLog.i_redcap(i)                 = {btwn_this_and_next};
    
                % take ALL stim params and add wrt original redcap table
                vars = stimLog.Properties.VariableNames(1:end-3);
                
                for h = 1 : length(vars)
                    for k = 1 : length(btwn_this_and_next)
    
                        redcap.(vars{h})(btwn_this_and_next(k)) = stimLog.(vars{h})(i);
    
                    end
                end
        
            else 
        
                stimLog.redcap_btwn_stimLogs{i}     = 'none_btwn_stimLogs';
                stimLog.i_redcap(i)           = {NaN};
            end
    
        % for last stimLog associate all reports from then to NOW rather than
        % next stimLog
        else
    
            btwn_this_and_next = ...
                find(ge(redcap.time, stimLog.time_stimLog(i)) & ...
                le(redcap.time, datetime('now','TimeZone','America/Los_Angeles')));
        
            if ~isempty(btwn_this_and_next)
            
                stimLog.redcap_btwn_stimLogs{i}     = 'btwn_stimLogs_and_now';
                
                stimLog.redcap_reports{i}     = redcap(btwn_this_and_next,:);
                stimLog.i_redcap(i)           = {btwn_this_and_next};
        
            else 
        
                stimLog.redcap_btwn_stimLogs{i}     = 'none_btwn_stimLogs_and_now';
            end
        end
    end
end

% remove redundant fields
stimLog  = removevars(stimLog, {'stimParams_prog1'});
redcap   = removevars(redcap, {'stimParams_prog1'});


% verify REDcap to StimLog assignment
all_i_redcap      = vertcat(stimLog.i_redcap{:});

if length(all_i_redcap(~isnan(all_i_redcap))) == length(unique(all_i_redcap(~isnan(all_i_redcap)))) 
    
    disp(strjoin([cfg.pt_id, ' | ','All REDcap report(s) assigned to unique stim settings.']))
else
    error('REDcap report(s) are assigned to multiple stim settings (RBL message).')
end

per_assigned = length(unique(all_i_redcap(~isnan(all_i_redcap)))) ./ height(redcap) * 100;


disp(strjoin([cfg.pt_id, ' | ', num2str(per_assigned), '% of REDcap report(s) assigned to stim settings.']));

end