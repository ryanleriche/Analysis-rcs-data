function [adaptiveLogTable,...
    rechargeSessions, ...
    group_changes,...
    adaptiveDetectionEvents,...
    mirrorLog] = read_adaptive_txt_log(cfg, fn)



%% RBL trouble-shooting as script rather than function

% fn = fullfile(f.folder, f.name); % specified in 'RCS_logs()'
% 
% cfg;                             % specified in 'RCS_logs()'

%This is used to extract ambulatory program changes, adaptive stim events
%battery recharge etc
%
%
%
% Written by Roee Gilron
% Optimized (removed almost all for loops) by Prasad Shirvalkar
% July 28 2021

% initialize table

group_changes = table();
rechargeSessions=[];

adaptiveDetectionEvents=[];
mirrorLog = [];

adaptiveLogTable = table();

% get time
str = fileread( fn );

if ~isempty(str)  % only continue if the file is not empty
    
    newBlocks = regexp(str, {'\n\r'});
    newBlockLines = newBlocks{1};
    newBlockLines = [1 newBlockLines];
    
    
    % loop on text and get each new block in a cell array
    cntBlock = 1;
    while cntBlock ~= (length(newBlockLines)-1)
        events{cntBlock} = str(newBlockLines(cntBlock) : newBlockLines(cntBlock+1));
        cntBlock = cntBlock + 1;
    end
    eventsRaw = events;
    
    
    %% get all event types
    xpruse1 = '(';
    cac1 = cellfun(@(x) regexp(x, xpruse1),events,'UniformOutput',false);
    xpruse1 = ')';
    cac2 = cellfun(@(x) regexp(x, xpruse1),events,'UniformOutput',false);
    
    strraw = cellfun(@(x,a,b) x(a(2)+1:b(2)-1),events,cac1,cac2,'UniformOutput',false);
    adaptiveLogEvents.EventID = strraw;
    
    allEvents = eventsRaw;



    
    if cfg.pull_adpt_logs
    %% AdaptiveTherapyStateChange
    i_ActiveDeviceChanged = strcmp(adaptiveLogEvents.EventID,'AdaptiveTherapyStateChange');
    
    events = allEvents(i_ActiveDeviceChanged);
    
    adaptiveLogTable = table('Size', [length(events) 10],'VariableTypes',{'datetime','double','double','double','double','double','double','double','double','cell'},'VariableNames',{'time','status','newstate','oldstate','prog0','prog1','prog2','prog3','rateHz','EventID'});
    
    
    startTimeDt  = get_date_from_hexstring(events); % see subfunction below
    adaptiveLogTable.time = startTimeDt;
    
    
    % get status
    xpr = 'AdaptiveTherapyModificationEntry.Status ';
    cac1 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
    xpr = '(EmbeddedActive)';
    cac2 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
    
    group_status = cellfun(@(x,a,b) x(a+68:b-3),events,cac1,cac2,'UniformOutput',false);
    
    statusdec = hex2dec(group_status);
    adaptiveLogTable.status = statusdec;
    
    
    % new state
    xpr = 'AdaptiveTherapyModificationEntry.NewState ';
    cac1 =  cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
    
    newstate = cellfun(@(x,a) x(a+68:a+69),events,cac1,'UniformOutput',false);
    newstate = hex2dec(newstate);
    adaptiveLogTable.newstate = newstate;
    
    % old state
    xpr = 'AdaptiveTherapyModificationEntry.OldState ';
    cac1 =  cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
    
    oldstate = cellfun(@(x,a) x(a+68:a+69),events,cac1,'UniformOutput',false);
    %
    %     for t = 1:length(cac1)
    %         if (cac1(t)+69) > length(str)
    %             oldstate(t,:) = NaN;
    %         else
    %             oldstate(t,:) = str(cac1(t)+68:cac1(t)+69);
    %         end
    %     end
    oldstate = hex2dec(oldstate);
    adaptiveLogTable.oldstate = oldstate;
    
    % loop on programs
    for p = 0:3
        xpruse = sprintf('AdaptiveTherapyModificationEntry.Prog%dAmpInMillamps ',p);
        cac1 =  cellfun(@(x) regexp(x, xpruse),events,'UniformOutput',false);
        
        clear progNum
        
        prog = cellfun(@(x,a) x(a+66:a+71),events,cac1,'UniformOutput',false);
        
        progNum = str2double(prog);
        fnuse = sprintf('prog%d',p);
        adaptiveLogTable.(fnuse) = progNum';
    end
    
    % rate
    xpruse = 'AdaptiveTherapyModificationEntry.RateAtTimeOfModification ';
    cac1 =  cellfun(@(x) regexp(x, xpruse),events,'UniformOutput',false);
    
    
    
    rate = cellfun(@(x,a) x(a+66:a+73),events,cac1,'UniformOutput',false);
    ratenum = str2double(rate);
    adaptiveLogTable.rateHz = ratenum';
    
    % events ID
    xpruse1 = 'CommonLogPayload`1.EventId      = 0x00 (';
    cac1 =  cellfun(@(x) regexp(x, xpruse1),events,'UniformOutput',false);
    xpruse2 = 'CommonLogPayload`1.EntryPayload = ';
    cac2 =  cellfun(@(x) regexp(x, xpruse2),events,'UniformOutput',false);
    
    strraw = cellfun(@(x,a,b) x(a:b-4),events,cac1,cac2,'UniformOutput',false);
    %       strraw = str(cac1:cac2-4);
    strtmp = erase(strraw,{xpruse1,')'});
    adaptiveLogTable.EventID = string(cellfun(@(x) x(1:end-3),strtmp,'UniformOutput',false))';
    end
    
    
    if cfg.pull_recharge_sess
    %% Recharge sessions
    
    i_ActiveDeviceChanged = strcmp(adaptiveLogEvents.EventID,'RechargeSesson');
    allEvents = eventsRaw;
    events = allEvents(i_ActiveDeviceChanged);
    rechargeSessions = table('Size', [length(events) 2], 'VariableTypes', {'datetime','cell'}, 'VariableNames', {'time','status'});
    
    startTimeDt  = get_date_from_hexstring(events); % see subfunction below
    rechargeSessions.time = startTimeDt;
    
    
    % get type
    xpr = 'RechargeSessionEventLogEntry.RechargeSessionStatus = ';
    cac1 =  cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
    xpr = 'RechargeSessionEventLogEntry.Unused';
    cac2 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
    
    group_status =  cellfun(@(x,a,b) x(a+59:b-12),events,cac1,cac2,'UniformOutput',false);
    
    rechargeSessions.status = group_status';
    end
    %%
    
    
    
    if cfg.pull_adpt_logs
    %% Adaptive therapy status
    i_ActiveDeviceChanged = strcmp(adaptiveLogEvents.EventID,'AdaptiveTherapyStatusChanged');
    allEvents = eventsRaw;
    events = allEvents(i_ActiveDeviceChanged);
    adaptiveStatus = table('Size', [length(events) 2], 'VariableTypes', {'datetime','cell'}, 'VariableNames', {'time','status'});
    
    % if ~isempty(events)
    
    startTimeDt  = get_date_from_hexstring(events); % see subfunction below
    adaptiveStatus.time = startTimeDt;
    
    % get type
    xpr = 'AdaptiveTherapyStatusChangedEventLogEntry.Status = ';
    cac1 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
    xpr = 'AdaptiveTherapyStatusChangedEventLogEntry.Unused = ';
    cac2 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
    
    clear group_status
    group_status = cellfun(@(x,a,b) x(a+57:b-12),events,cac1,cac2,'UniformOutput',false);
    adaptiveStatus.status = group_status';
    % end
    end
    
    
    if cfg.pull_event_logs
    %% Group ChangadaptiveLogEvents.EventIDes - i.e. ActiveDeviceChanged
% unique(adaptiveLogEvents.EventID)

    i_ActiveDeviceChanged = strcmp(adaptiveLogEvents.EventID,'ActiveDeviceChanged');
 
    i_TherapyStatus       = strcmp(adaptiveLogEvents.EventID, 'TherapyStatus');

    allEvents = eventsRaw;
    events = allEvents(i_ActiveDeviceChanged);

    ActiveGroup = table('Size', [length(events) 2], 'VariableTypes', {'datetime','string'}, 'VariableNames', {'time','event'});
    % if ~isempty(events)
    
    startTimeDt  = get_date_from_hexstring(events); % see subfunction below
    ActiveGroup.time = startTimeDt;
    %%
    % get type

    xpr = 'TherapyActiveGroupChangedEventLogEntry.NewGroup = ';
    cac1 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
    
    xpr = 'TherapyActiveGroupChangedEventLogEntry.Unused   = ';
    cac2 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);

    clear group_status
    group_status = cellfun(@(x,a,b) x(a+56:b-12), events, cac1, cac2,'UniformOutput',false);
    
    groupUse = replace(group_status,{'Group0','Group1','Group2','Group3'},{'GroupA','GroupB','GroupC','GroupD'});
    ActiveGroup.event = groupUse';


    xpr      = 'TherapyStatusEventLogEntry.TherapyStatus';
    cac_3    = cellfun(@(x) regexp(x, xpr), eventsRaw(i_TherapyStatus), 'UniformOutput',false);

    OnOff = cellfun(@(x,a) x(a(2):a(2)+57),  eventsRaw(i_TherapyStatus), cac_3, 'UniformOutput',false);
 
    TherapyStatus = table;

    % simplifies name and trims whitespace
    TherapyStatus.event = regexprep(replace(OnOff,...
        {'TherapyStatusEventLogEntry.TherapyStatus     = 0x01 (On)',...
        'TherapyStatusEventLogEntry.TherapyStatus     = 0x00 (Off)'},{'on','off'})',...
         '\s', '');

  

    TherapyStatus.time = get_date_from_hexstring(allEvents(i_TherapyStatus));

    group_changes = sortrows([TherapyStatus; ActiveGroup],2);


    % remove events that repeat
    group_changes = [group_changes(1,:);...
                        group_changes(...
                            ~strcmp(group_changes.event(1:end -1),...
                            group_changes.event(2:end)),:)...
                            ];

    % accounts for if first and second row repeat

    if height(group_changes) > 1
        if strcmp(group_changes.event(1), group_changes.event(2))
    
            group_changes(1,:) = [];
    
        end
    end
    % end
    end
    %%
    
    
    if cfg.ld_detection_events

    %% LD detection events
    i_ActiveDeviceChanged = strcmp(adaptiveLogEvents.EventID,'LdDetectionEvent');
    allEvents = eventsRaw;
    events = allEvents(i_ActiveDeviceChanged);
    adaptiveDetectionEvents = table('Size', [length(events) 5], 'VariableTypes', {'datetime','double','cell','double','cell'}, 'VariableNames', {'time','detectionStatus','detectionText','previousDetectionStatus','previousDetectionText'});
    
    % for e = 1:length(events)
    %     str = events{e};
    %     car = regexp(str, '\r');
    
    
    startTimeDt  = get_date_from_hexstring(events); % see subfunction below
    adaptiveDetectionEvents.time = startTimeDt;
    
    %     GET current detection state
    
    % get type
    xprC = 'LdDetectionEntry.CurrentDetectionState  = ';
    cac1 = cellfun(@(x) regexp(x, xprC),events,'UniformOutput',false);
    xprP = 'LdDetectionEntry.PreviousDetectionState = ';
    cac2 = cellfun(@(x) regexp(x, xprP),events,'UniformOutput',false);
    
    % a few states possible
    tempstr = cellfun(@(x,a,b) x(a:b),events,cac1,cac2,'UniformOutput',false);
    detectionNum  = get_detection_Num(tempstr);
    % get string event
    newstr = get_newstr(tempstr);
    
    adaptiveDetectionEvents.detectionStatus = detectionNum;
    adaptiveDetectionEvents.detectionText = newstr';
    
    
    %       GET the previous detection state
    
    % get type
    cac1 = cellfun(@(x) regexp(x, xprP),events,'UniformOutput',false);  %previous Detection state from above
    xpr = 'LdDetectionEntry.Unused';
    cac2 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
    
    % a few states possible
    tempstr = cellfun(@(x,a,b) x(a:b),events,cac1,cac2,'UniformOutput',false);
    detectionNum  = get_detection_Num(tempstr);
    % get string event
    newstr = get_newstr(tempstr);
    
    adaptiveDetectionEvents.previousDetectionStatus = detectionNum;
    adaptiveDetectionEvents.previousDetectionText = newstr';
    % end
 
    end
    %%
    
    % COMMENTED OUT by prasad because this is redundant with above, except for the 'AdaptiveTherapyStateWritten' Index
    % %% Recharge sessions
    % idxuse = strcmp(adaptiveLogEvents.EventID,'AdaptiveTherapyStateWritten');
    % allEvents = eventsRaw;
    % events = allEvents(idxuse);
    % adaptiveStateWritten = table();
    % for e = 1:length(events)
    %     str = events{e};
    %     car = regexp(str, '\r');
    %
    %     xpr = 'Seconds = ';
    %     cac1 = regexp( str, xpr );
    %
    %     xpr = 'DateTime = ';
    %     cac2 = regexp( str, xpr );
    %
    %     clear hexstr
    %     for t = 1:length(cac1)
    %         hexstr(t,:) = str(cac1(t)+12:cac2(t)-3);
    %     end
    %     rawsecs = hex2dec(hexstr);
    %     startTimeDt = datetime(datevec(rawsecs./86400 + datenum(2000,3,1,0,0,0))); % medtronic time - LSB is seconds
    %     rechargeSessions.time(e) = startTimeDt;
    %     %%
    %
    %     %% get type
    %     xpr = 'RechargeSessionEventLogEntry.RechargeSessionStatus = ';
    %     cac1 = regexp( str, xpr );
    %
    %
    %     xpr = 'RechargeSessionEventLogEntry.Unused';
    %     cac2 = regexp( str, xpr );
    %     clear status
    %     status = str(cac1+59:cac2-12);
    %
    %     rechargeSessions.status{e} = status;
    %     %%
    % end
    
    
    
    if size(adaptiveLogTable,1) > 30
        at = adaptiveLogTable(1:20,:);
        idxzero = at.newstate==0;
        unique(at.prog0(idxzero));
    end
else
    fprintf('Empty Log file detected../n')
    adaptiveLogTable=[];
    rechargeSessions=[];
    ActiveGroup=[];
    adaptiveDetectionEvents=[];
    return
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS DECLARED HERE FOR USE IN LD Detection events

    function detectionNum_out  = get_detection_Num(tempstr_input)
        xpr0 = '0x';
        cac11 =  cellfun(@(x) regexp(x, xpr0),tempstr_input,'UniformOutput',false);
        xpr0 = '(';
        cac22 =  cellfun(@(x) regexp(x, xpr0),tempstr_input,'UniformOutput',false);
        
        detHexString = cellfun(@(x,a,b) x(a+2:b-2),tempstr_input,cac11,cac22,'UniformOutput',false);
        detectionNum_out = hex2dec(detHexString);
    end

    function newstr_out = get_newstr(tempstr_input)
        xpr0 = '(';
        cac11 = cellfun(@(x) regexp(x, xpr0),tempstr_input,'UniformOutput',false);
        xpr0 = ')';
        cac22 = cellfun(@(x) regexp(x, xpr0),tempstr_input,'UniformOutput',false);
        
        newstr_out = cellfun(@(x,a,b) x(a:b),tempstr_input,cac11,cac22,'UniformOutput',false);
    end

    function date_out  = get_date_from_hexstring(events_input)
        xpr0 = 'Seconds = ';
        cac11 = cellfun(@(x) regexp(x, xpr0),events_input,'UniformOutput',false);
        xpr0 = 'DateTime = ';
        cac22 = cellfun(@(x) regexp(x, xpr0),events_input,'UniformOutput',false);
        
        hexstr0 = cellfun(@(x,a,b) x(a+12:b-3),events_input,cac11,cac22,'UniformOutput',false);
        
        rawsecs0 = hex2dec(hexstr0);
        date_out = datetime(datevec(rawsecs0./86400 + datenum(2000,3,1,0,0,0))); % medtronic time - least signififant bit (LSB) is seconds
        
    end

%%
end