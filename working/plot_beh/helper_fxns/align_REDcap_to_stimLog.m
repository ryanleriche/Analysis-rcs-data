function [stimLog_w_redcap] = align_REDcap_to_stimLog(cfg, db_RCSXXX, redcap)
% cfg                 = [];
% cfg.stage_dates     = stage_dates{2};
% cfg.pt_id           = 'RCS02';
% 
% 
% db_RCSXXX           = db.RCS02R;
% redcap              = REDcap.RCS02;
% 


% only takes REDcap surverys from AFTER Stage I implant
redcap   = redcap(ge(redcap.time, datetime(cfg.stage_dates(1),'TimeZone','America/Los_Angeles')),...
                  :); 

i_stimlog = cellfun(@(x) ~isempty(x), db_RCSXXX.stimLogSettings);

db_RCSXXX  = db_RCSXXX(i_stimlog,:);

if strcmp('RCS07', cfg.pt_id) 
    i_wo_time_stimLog   = cellfun(@(x) ~any(strcmp(x.Properties.VariableNames,'time_stimLog')), db_RCSXXX.stimLogSettings);
    i_w_time_stimLog    = find(~i_wo_time_stimLog);
    i_wo_time_stimLog    = find(i_wo_time_stimLog);
    
    
    [~, i_near] = min(abs(i_wo_time_stimLog' - i_w_time_stimLog));
    
    
    near_UTCoffset = cellfun(@(x) sprintf('%+03.0f:00', x.UTCoffset), db_RCSXXX.metaData(i_w_time_stimLog(i_near)), 'UniformOutput', false);
    
    time_stimLog   = cellfun(@(x) {x.HostUnixTime ./ 1000}, db_RCSXXX.stimLogSettings(i_wo_time_stimLog));
    
    % takes timezone from next StimLog's timezone
    now_w_datetime = cellfun(@(x,y) datetime(x,'ConvertFrom','posixTime','TimeZone', ...
                        y,'Format','dd-MMM-yyyy HH:mm:ss.SSS'),  time_stimLog, near_UTCoffset);
    
    j = 1;
    for i = i_wo_time_stimLog'
    
        db_RCSXXX.stimLogSettings{i,1}.time_stimLog = now_w_datetime(j);
        j = j+1;
    
    end

    for i =  1 : height(db_RCSXXX)
 
        if db_RCSXXX.stimSettingsOut{i,1}.cyclingEnabled{1}

            db_RCSXXX.stimLogSettings{i,1}.cycleOnTime = ...
                repmat(db_RCSXXX.stimSettingsOut{i,1}.cycleOnTime{1}, height(db_RCSXXX.stimLogSettings{i,1}), 1);
    
            db_RCSXXX.stimLogSettings{i,1}.cycleOffTime = ...
                repmat(db_RCSXXX.stimSettingsOut{i,1}.cycleOffTime{1}, height(db_RCSXXX.stimLogSettings{i,1}), 1);

        else
             db_RCSXXX.stimLogSettings{i,1}.cycleOnTime = ...
                NaN(height(db_RCSXXX.stimLogSettings{i,1}), 1);
    
            db_RCSXXX.stimLogSettings{i,1}.cycleOffTime = ...
                NaN(height(db_RCSXXX.stimLogSettings{i,1}), 1);

        end

    end

        
else
    
    for i =  1 : height(db_RCSXXX)
    
        if ~any(strcmp(db_RCSXXX.stimLogSettings{i,1}.Properties.VariableNames, 'time_stimLog'))
    
            % finds proper timezone for StimLog.json
            try
                time_stimLog = db_RCSXXX.stimLogSettings{i,1}.HostUnixTime / 1000;

                if ~isempty(db_RCSXXX.metaData{i+1})
                    % takes timezone from next StimLog's timezone
                    timeFormat = sprintf('%+03.0f:00', ...
                        db_RCSXXX.metaData{i+1}.UTCoffset);

                else
                    timeFormat = sprintf('%+03.0f:00', ...
                        db_RCSXXX.metaData{i-1}.UTCoffset);
                    
                end

            catch

                if i == 1
                    timeFormat = sprintf('%+03.0f:00', ...
                        db_RCSXXX.metaData{i+2}.UTCoffset);
                    
                elseif ~isempty(db_RCSXXX.metaData{i-1}) % for single case in RCS04R w/ empty metaData
                    % takes previous StimLog's timezone
                    timeFormat = sprintf('%+03.0f:00', ...
                        db_RCSXXX.metaData{i-1}.UTCoffset);
                
                else
                    timeFormat = sprintf('%+03.0f:00', ...
                        db_RCSXXX.metaData{i-2}.UTCoffset);
                
                end
            end

            db_RCSXXX.stimLogSettings{i,1}.time_stimLog = datetime(time_stimLog,...
                'ConvertFrom','posixTime','TimeZone',timeFormat,...
                'Format','dd-MMM-yyyy HH:mm:ss.SSS');
              
        end


 
        if db_RCSXXX.stimSettingsOut{i,1}.cyclingEnabled{1}

            db_RCSXXX.stimLogSettings{i,1}.cycleOnTime = ...
                repmat(db_RCSXXX.stimSettingsOut{i,1}.cycleOnTime{1}, height(db_RCSXXX.stimLogSettings{i,1}), 1);
    
            db_RCSXXX.stimLogSettings{i,1}.cycleOffTime = ...
                repmat(db_RCSXXX.stimSettingsOut{i,1}.cycleOffTime{1}, height(db_RCSXXX.stimLogSettings{i,1}), 1);

        else
             db_RCSXXX.stimLogSettings{i,1}.cycleOnTime = ...
                NaN(height(db_RCSXXX.stimLogSettings{i,1}), 1);
    
            db_RCSXXX.stimLogSettings{i,1}.cycleOffTime = ...
                NaN(height(db_RCSXXX.stimLogSettings{i,1}), 1);

        end
    end
end

% takes stimLogSettings to a single sorted timetable and removes redundant fields
stimLog_w_redcap = sortrows(...
                table2timetable(...
                    vertcat(db_RCSXXX.stimLogSettings{:})));

for i = 1 : height(stimLog_w_redcap)

    stim_para = strsplit(char(stimLog_w_redcap.stimParams_prog1(i)),',');
    
    if isempty(stim_para{1})
        stimLog_w_redcap.stimContacts{i} = '';
    else
        stimLog_w_redcap.stimContacts{i}     = stim_para{1};
    end

    stimLog_w_redcap.stimAmp(i)          = str2double(stim_para{2}(2:end -2));
    stimLog_w_redcap.stimPW(i)           = str2double(stim_para{3}(2:end -2));
    stimLog_w_redcap.stimfreq(i)         = str2double(stim_para{4}(2:end -2));


    if i < height(stimLog_w_redcap)
        btwn_this_and_next = ...
            find(ge(redcap.time, stimLog_w_redcap.time_stimLog(i)) & ...
            le(redcap.time, stimLog_w_redcap.time_stimLog(i + 1)));
    
        if ~isempty(btwn_this_and_next)
        
            stimLog_w_redcap.redcap_btwn_stimLogs{i}     = 'btwn_stimLogs';
            
            stimLog_w_redcap.redcap_reports{i}     = redcap(btwn_this_and_next,:);
            stimLog_w_redcap.i_redcap(i)           = {btwn_this_and_next};
    
        else 
    
            stimLog_w_redcap.redcap_btwn_stimLogs{i}     = 'none_btwn_stimLogs';
            stimLog_w_redcap.i_redcap(i)           = {NaN};
        end

    % for last stimLog associate all reports from then to NOW rather than
    % next stimLog
    else

        btwn_this_and_next = ...
            find(ge(redcap.time, stimLog_w_redcap.time_stimLog(i)) & ...
            le(redcap.time, datetime('now','TimeZone','America/Los_Angeles')));
    
        if ~isempty(btwn_this_and_next)
        
            stimLog_w_redcap.redcap_btwn_stimLogs{i}     = 'btwn_stimLogs_and_now';
            
            stimLog_w_redcap.redcap_reports{i}     = redcap(btwn_this_and_next,:);
            stimLog_w_redcap.i_redcap(i)           = {btwn_this_and_next};
    
        else 
    
            stimLog_w_redcap.redcap_btwn_stimLogs{i}     = 'none_btwn_stimLogs_and_now';
    
        end
    end
end

% remove redundant fields
stimLog_w_redcap = removevars(stimLog_w_redcap, {'GroupA','GroupB','GroupC',...
                      'GroupD','therapyStatus','HostUnixTime','stimParams_prog1','stimParams_prog2',...
                      'stimParams_prog3','stimParams_prog4','updatedParameters',...
                       });

all_i_redcap = vertcat(stimLog_w_redcap.i_redcap{:});

if length(all_i_redcap(~isnan(all_i_redcap))) == length(unique(all_i_redcap(~isnan(all_i_redcap)))) 
    
    disp([cfg.pt_id, ': ','All REDcap report(s) assigned to unique stim settings.'])
else
    error('REDcap report(s) are assigned to multiple stim settings.')
end

per_assigned = length(unique(all_i_redcap(~isnan(all_i_redcap)))) ./ height(redcap) * 100;


disp([cfg.pt_id, ': ', num2str(per_assigned), '% of REDcap report(s) assigned to stim settings.']);


end