function [eventLogTable] = createEventLogTable(folderPath)
%%
% Extract information from EventLog.json
%
% Input: Folder path to Device* folder containing json files
%%
% Load in EventLog.json
eventLog = deserializeJSON([folderPath filesep 'EventLog.json']);

if ~isempty(eventLog)
    
    eventLogTable = table();
    numRecords    = length(eventLog);
    
    for iRecord = 1:numRecords
        clear newEntry
        
        newEntry.SessionId          = eventLog(iRecord).RecordInfo.SessionId;
        newEntry.HostUnixTime       = eventLog(iRecord).RecordInfo.HostUnixTime;
        newEntry.EventName          = eventLog(iRecord).Event.EventName;
        newEntry.EventType          = eventLog(iRecord).Event.EventType;
        newEntry.EventSubType       = eventLog(iRecord).Event.EventSubType;
        newEntry.UnixOnsetTime      = eventLog(iRecord).Event.UnixOnsetTime;
        newEntry.UnixOffsetTime     = eventLog(iRecord).Event.UnixOffsetTime;
        
        eventLogTable               = addRowToTable(newEntry,eventLogTable); 

    end

    i_same_as_next = strcmp(eventLogTable.EventType(1:end), [eventLogTable.EventType(2:end);{''}]);

    % no need to return beeping timestamps
    i_beep         = strcmp(eventLogTable.EventType, 'Log Beep');

    beeps          = find(i_beep);
    
    if ~isempty(beeps)
        eventLogTable  = [eventLogTable(~i_beep | ~i_same_as_next,:); eventLogTable([beeps(1), beeps(end)],:)];
    end

    eventLogTable  = sortrows(eventLogTable, 'UnixOnsetTime');
end
end
