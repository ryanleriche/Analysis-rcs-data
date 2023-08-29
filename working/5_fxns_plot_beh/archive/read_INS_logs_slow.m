 function output = read_INS_logs_slow(fn)

% Written by Roee Gilron
%{

Configured by Ryan to:
    save INS logs based off of every single line of text

    captures all information saved in AppLogs and EventLogs
    
Nov. 2022
%}


% allows trouble-shooting as script rather than function:
%fn = fullfile(filelist.folder{i}, filelist.name{i});


txtlog = fileread(fn);


%%
if ~isempty(txtlog)  % only continue if the file is not empty
    
    newBlocks = regexp(txtlog, {'\n\r'});
    newBlockLines = newBlocks{1};
    newBlockLines = [1 newBlockLines];
    
    events_raw = {};
    % loop on text and get each new block in a cell array
    cntBlock = 1;
    while cntBlock ~= (length(newBlockLines)-1)
        events_raw{cntBlock,1} = txtlog(newBlockLines(cntBlock) : newBlockLines(cntBlock+1));
        cntBlock = cntBlock + 1;
    end
    
    %% get all event types

    % save in table
    events_tbl            = table;
    events_tbl.time       = get_date_from_hexstring(events_raw);

    cac1 = cellfun(@(x) regexp(x, '('), events_raw,'UniformOutput',false);
    cac2 = cellfun(@(x) regexp(x, ')'), events_raw,'UniformOutput',false);
   
    events_tbl.event_id   = ...
        cellfun(@(x,a,b) x(a(2)+1:b(2)-1), events_raw, cac1, cac2,'UniformOutput',false);


    cac1 = cellfun(@(x) regexp(x, 'CommonLogPayload`1.EntryPayload ='), events_raw,'UniformOutput',false);
    
    % return everything AFTER CommonLogPayload
    entry_info = cellfun(@(x,a) x(a+44:end), events_raw, cac1,'UniformOutput',false);

  %%

  
% RBL, 11/03/22 version commented out
%{
 
    after_date = cellfun(@(x) regexp(x,  'CommonLogPayload`1.EntryPayload ='), events_raw, 'UniformOutput', false);

    temp_event_info = cellfun(@(x, b) x(b+44:end), events_raw, after_date, 'UniformOutput', false);


%     Entry =  cellfun(@(x) regexp(x,  'Entry'), temp_event_info, 'UniformOutput', false);




    events_tbl.event_id   = temp_event_ids';

    events_tbl.event_info = temp_event_info';

%%
    % just from parentheses, mine out all text
    for i = 1 : height(events_tbl)

        if length(cac1{i}) >= 3
            desc_txt = [];
            try
                for j = 3 : length(cac1{i})
        
                    temp = events_raw{i}(cac1{i}(j)+1: cac2{i}(j)-1);
                    if j == 3
                        desc_txt = temp;
                    else
                        desc_txt = [desc_txt,'_' temp];
                    end
                    
                end
            
            catch
                desc_txt = 'mismatched parentheses';
            end
        else
            desc_txt = 'None';
        end

        events_tbl.entry{i} = desc_txt ;
    end

%}  
  


    if endsWith(fn,"EventLog.txt")

        cac1 = cellfun(@(x) regexp(x, 'LogEntry'), entry_info,'UniformOutput',false);

        log_entries = cellfun(@(x,a) (x(1:a+7)), entry_info, cac1,'UniformOutput',false);
        
        % for completeness, include the log entry names
        events_tbl.entry_names = log_entries;
        
        
        % but NOW loop through and get entries themselves--use newline and equal sign characters!
        
        % seperate by newline characters
        i_new_line = cellfun(@(x) regexp(x, newline), entry_info, 'UniformOutput', false);


        for j = 1 : height(events_tbl)
        
            temp_fields = entry_info{j};
            temp_i_lines = [1,i_new_line{j}];
            
            temp_line    = {};
        
            for h = 1 : length(temp_i_lines) - 1
            
                temp_line{h,1} = strtrim(temp_fields(temp_i_lines(h) : temp_i_lines(h+1) -2));
            
            end
    
            % special case for 'ClockChanged', since it does not obey equal sign
            % plus newline motiff
            if any(contains(temp_line, 'ClockChanged'))
    
    
                events_tbl.entries(j) = {struct('time', ...
                                     get_date_from_hexstring(temp_line(2)))};
            else
        
                i_LogEntry = cellfun(@(x) regexp(x, 'LogEntry'),...
                                     temp_line, 'UniformOutput', false);
            
                fields = cellfun(@(x,a) x(a+9:end), ...
                                 temp_line, i_LogEntry, 'UniformOutput', false);
            
            
                i_equals = cellfun(@(x) regexp(x, '='), ...
                                   fields, 'UniformOutput', false);
            
                names = strtrim(...
                            cellfun(@(x,a) x(1:a-1), ...
                                    fields, i_equals, 'UniformOutput', false))';
            
                values = strtrim(...
                             cellfun(@(x,a) x(a+2:end), ...
                                     fields, i_equals, 'UniformOutput', false))';
            
                args    = [names; values];
            
                events_tbl.entries(j) = {struct(args{:})};
            end
        end


    elseif endsWith(fn, "AppLog.txt")
%%
        cac1 = cellfun(@(x) regexp(x, 'Entry'), entry_info,'UniformOutput',false);
        
        log_entries = cellfun(@(x,a) (x(1:a+4)), entry_info, cac1,'UniformOutput',false);
        
        % for completeness, include the log entry names
        events_tbl.entry_names = log_entries;
        
        
        % but NOW loop through and get entries themselves--use newline and equal sign characters!
        
        % seperate by newline characters
        i_new_line = cellfun(@(x) regexp(x, newline), entry_info, 'UniformOutput', false);
    
    
        for j = 1 : height(events_tbl)
        
            temp_fields = entry_info{j};
            temp_i_lines = [1,i_new_line{j}];
            
            temp_line    = {};
        
            for h = 1 : length(temp_i_lines) - 1
            
                temp_line{h,1} = strtrim(temp_fields(temp_i_lines(h) : temp_i_lines(h+1) -2));
            
            end
    
            % special case for 'ClockChanged', since it does not obey equal sign
            % plus newline motiff
            if any(contains(temp_line, 'ClockChanged'))
    
    
                events_tbl.entries(j) = {struct('time', ...
                                     get_date_from_hexstring(temp_line(2)))};
            else
            try
                i_LogEntry = cellfun(@(x) regexp(x, 'Entry'),...
                                     temp_line, 'UniformOutput', false);
            
                fields = cellfun(@(x,a) x(a+6:end), ...
                                 temp_line, i_LogEntry, 'UniformOutput', false);
            
            
                i_equals = cellfun(@(x) regexp(x, '='), ...
                                   fields, 'UniformOutput', false);
            
                names = strtrim(...
                            cellfun(@(x,a) x(1:a-1), ...
                                    fields, i_equals, 'UniformOutput', false))';
            
                values = strtrim(...
                             cellfun(@(x,a) x(a+2:end), ...
                                     fields, i_equals, 'UniformOutput', false))';
            
                args    = [names; values];
            
                events_tbl.entries(j) = {struct(args{:})};
            catch
                events_tbl.entries(j) = {'parsing error'};

            end
            end
        end


    %% 
    %{
    % AdaptiveTherapyStateChange

    i_ActiveDeviceChanged = strcmp(log_events.EventID,'AdaptiveTherapyStateChange');
    
    events = all_events(i_ActiveDeviceChanged);
    
    adaptiveLogTable = table('Size', [length(events) 10],...
        'VariableTypes',{'datetime','double','double','double','double','double','double','double','double','cell'},...
        'VariableNames',{'time','status','newstate','oldstate','prog0mA','prog1mA','prog2mA','prog3mA','rateHz','EventID'});
    
    
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
        fnuse = sprintf('prog%dmA',p);
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
    
    temp_event_ids = cellfun(@(x,a,b) x(a:b-4),events,cac1,cac2,'UniformOutput',false);
    %       strraw = str(cac1:cac2-4);
    strtmp = erase(temp_event_ids,{xpruse1,')'});
    adaptiveLogTable.EventID = string(cellfun(@(x) x(1:end-3),strtmp,'UniformOutput',false))';
    %}
%%
        % Adaptive therapy status
        %{
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
        %}
    
    
    
        % Recharge sessions
        %{
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
        %}
    
    
        % LD detection events
        %{
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
        %}


    elseif endsWith(fn, "MirrorLog.txt")



    end

output = events_tbl;
        
else
    output = {'empty'};
end 

% % SUBFUNCTIONS DECLARED HERE FOR USE IN LD Detection events
    
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

end