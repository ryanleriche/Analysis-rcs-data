 function varargout  = read_INS_logs_fast(fn)
%{
Written by Roee Gilron

Optimized (removed almost all for loops) by Prasad Shirvalkar
July 28 2021
%}


% trouble-shooting as script rather than function:

%fn = fullfile(filelist.folder{i}, filelist.name{i});


str = fileread(fn);


%%
if ~isempty(str)  % only continue if the file is not empty
    
    newBlocks    = regexp(str, {'\n\r'});
    newBlockLines = newBlocks{1};
    newBlockLines = [1 newBlockLines];
    
    events_raw = {};
    % loop on text and get each new block in a cell array
    cntBlock = 1;
    while cntBlock ~= (length(newBlockLines)-1)
        events_raw{cntBlock,1} = str(newBlockLines(cntBlock) : newBlockLines(cntBlock+1));
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

        % groups (A, B, D, or D)
        i_ActiveDeviceChanged = strcmp(events_tbl.event_id, 'ActiveDeviceChanged');
        events                = events_raw(i_ActiveDeviceChanged);

        ActiveGroup = table('Size', [length(events) 2], ...
                             'VariableTypes', {'datetime','string'}, 'VariableNames', {'time','event'});


        ActiveGroup.time  = get_date_from_hexstring(events); 

        xpr = 'TherapyActiveGroupChangedEventLogEntry.NewGroup = ';
        cac1 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
        
        xpr = 'TherapyActiveGroupChangedEventLogEntry.Unused   = ';
        cac2 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
    
        clear group_status
        group_status = cellfun(@(x,a,b) x(a+56:b-12), events, cac1, cac2,'UniformOutput',false);
        
        groupUse = replace(...
             group_status,{'Group0','Group1','Group2','Group3'}, {'GroupA','GroupB','GroupC','GroupD'});
        
        ActiveGroup.event = groupUse;

        % therapy status (On, Off, or Lead Integrity test)
        i_TherapyStatus       = strcmp(events_tbl.event_id, 'TherapyStatus');
        events                = events_raw(i_TherapyStatus);

        xpr      = 'TherapyStatusEventLogEntry.TherapyStatus';
        cac_3    = cellfun(@(x) regexp(x, xpr), events , 'UniformOutput',false);
    
        OnOff    = cellfun(@(x,a) x(a(2):a(2)+57),  events  , cac_3, 'UniformOutput',false);
     
        TherapyStatus = table;


        TherapyStatus.time = get_date_from_hexstring(events);
    
        % simplifies name and trims whitespace
        TherapyStatus.event = regexprep(replace(OnOff,...
            {'TherapyStatusEventLogEntry.TherapyStatus     = 0x01 (On)',...
            'TherapyStatusEventLogEntry.TherapyStatus     = 0x00 (Off)'},{'On','Off'})',...
             '\s', '')';

        TherapyStatus.event(contains(TherapyStatus.event, 'LeadI')) = {'LeadIntegrityTest'};

       % return group changes WITH their therapy status
       group_changes = sortrows([TherapyStatus; ActiveGroup],'time');
    
       % remove events that repeat
       group_changes = [group_changes(1,:);...
                            group_changes(...
                                ~strcmp(group_changes.event(1:end -1),...
                                group_changes.event(2:end)),:)...
                                ];
       % accounts for when first and second row repeat
       if height(group_changes) > 1
          if strcmp(group_changes.event(1), group_changes.event(2))
        
              group_changes(1,:) = [];
        
          end
       end

       if isempty(group_changes)
            group_changes = 'empty';
       end

       %%% Recharge sess[i]ons
       i_rech    = strcmp(events_tbl.event_id, 'RechargeSesson');
       events    = events_raw(i_rech);
    
       if ~isempty(events)
    
            rech_sess = table('Size', [length(events) 2], 'VariableTypes', {'datetime','cell'}, 'VariableNames', {'time','status'});
            
            rech_sess.time  = get_date_from_hexstring(events); % see subfunction below
                         
            
            % get type
            xpr  = 'RechargeSessionEventLogEntry.RechargeSessionStatus = ';
            cac1 =  cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
            xpr  = 'RechargeSessionEventLogEntry.Unused';
            cac2 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
            
            rech_sess.status =...
                cellfun(@(x,a,b) x(a+59:b-12),events,cac1,cac2,'UniformOutput',false);
            
        else
            rech_sess = 'empty';
        end
    

        varargout{1}  = group_changes;
        varargout{2}  = rech_sess;


elseif endsWith(fn, "AppLog.txt")
    %%
    %%% AdaptiveTherapyStateChange
    i_change  = strcmp(events_tbl.event_id,'AdaptiveTherapyStateChange');
    events    = events_raw(i_change);
    
    if ~isempty(events)
        aDBS_log_tbl      = table('Size', [length(events) 10],...
            'VariableTypes',{'datetime','double','double','double','double','double','double','double','double','cell'},...
            'VariableNames',{'time','status','newstate','oldstate','prog0mA','prog1mA','prog2mA','prog3mA','rateHz','event_id'});
        
        
        aDBS_log_tbl.time  = get_date_from_hexstring(events); % see subfunction below
     
        % get status
        xpr  = 'AdaptiveTherapyModificationEntry.Status ';
        cac1 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
        xpr  = '(EmbeddedActive)';
        cac2 = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
        
        aDBS_log_tbl.status =...
            hex2dec(...
                cellfun(@(x,a,b) x(a+68:b-3),events,cac1,cac2,'UniformOutput',false));

        
        
        % new state
        xpr      = 'AdaptiveTherapyModificationEntry.NewState ';
        cac1     =  cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
        
        newstate = cellfun(@(x,a) x(a+68:a+69),events,cac1,'UniformOutput',false);
        newstate = hex2dec(newstate);

        aDBS_log_tbl.newstate = newstate;
        
        % old state
        xpr      = 'AdaptiveTherapyModificationEntry.OldState ';
        cac1     =  cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
        
        oldstate = cellfun(@(x,a) x(a+68:a+69),events,cac1,'UniformOutput',false);
 
        aDBS_log_tbl.oldstate = hex2dec(oldstate);
        
        % loop on programs
        for p = 0:3
            xpruse = sprintf('AdaptiveTherapyModificationEntry.Prog%dAmpInMillamps ',p);
            cac1 =  cellfun(@(x) regexp(x, xpruse),events,'UniformOutput',false);
            
            
            prog = cellfun(@(x,a) x(a+66:a+71),events,cac1,'UniformOutput',false);
            
            progNum = str2double(prog);
            fnuse   = sprintf('prog%dmA',p);
            aDBS_log_tbl.(fnuse) = progNum;

            clear progNum
        end
        
        % rate
        xpruse = 'AdaptiveTherapyModificationEntry.RateAtTimeOfModification ';
        cac1   =  cellfun(@(x) regexp(x, xpruse),events,'UniformOutput',false);
        
        
        
        rate    = cellfun(@(x,a) x(a+66:a+73),events,cac1,'UniformOutput',false);
        aDBS_log_tbl.rateHz  = str2double(rate);

        
        % events ID
        xpruse1 = 'CommonLogPayload`1.EventId      = 0x00 (';
        cac1    =  cellfun(@(x) regexp(x, xpruse1),events,'UniformOutput',false);
        xpruse2 = 'CommonLogPayload`1.EntryPayload = ';
        cac2    =  cellfun(@(x) regexp(x, xpruse2),events,'UniformOutput',false);
        
        temp_event_ids = cellfun(@(x,a,b) x(a:b-4),events,cac1,cac2,'UniformOutput',false);
        %       strraw = str(cac1:cac2-4);
        strtmp = erase(temp_event_ids,{xpruse1,')'});

        aDBS_log_tbl.event_id ...
            = string(cellfun(@(x) x(1:end-3),strtmp,'UniformOutput',false));

    else
        aDBS_log_tbl = {'empty'};
    end

    
    %%% Adaptive therapy status
    %{
    i_stat   = strcmp(events_tbl.event_id,'AdaptiveTherapyStatusChanged');
    events   = events_raw(i_stat);

    if ~isempty(events)

        adapt_stat = table('Size', [length(events) 2], 'VariableTypes', {'datetime','cell'}, 'VariableNames', {'time','status'});
        
        
        adapt_stat.time   = get_date_from_hexstring(events); % see subfunction below
        
        % get type
        xpr     = 'AdaptiveTherapyStatusChangedEventLogEntry.Status = ';
        cac1    = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
        xpr     = 'AdaptiveTherapyStatusChangedEventLogEntry.Unused = ';
        cac2    = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
        
        status  = cellfun(@(x,a,b) x(a+57:b-12),events,cac1,cac2,'UniformOutput',false);
        adapt_stat.status = status';

    else
        adapt_stat = {'empty'};
    end
    %}
    
    %%% LD detection events
    i_ld    = strcmp(events_tbl.event_id,'LdDetectionEvent');
    events  = events_raw(i_ld);

    if ~isempty(events)

        ld_detect       = table('Size', [length(events) 5], ...
            ...
            'VariableTypes',...
                 {'datetime','double','cell','double','cell'}, ...
            'VariableNames', ...
                 {'time','detectionStatus','detectionText','previousDetectionStatus','previousDetectionText'});
               
        ld_detect.time  = get_date_from_hexstring(events); % see subfunction below
        
        % GET current detection state
        % get type
        xprC = 'LdDetectionEntry.CurrentDetectionState  = ';
        cac1 = cellfun(@(x) regexp(x, xprC),events,'UniformOutput',false);
        xprP = 'LdDetectionEntry.PreviousDetectionState = ';
        cac2 = cellfun(@(x) regexp(x, xprP),events,'UniformOutput',false);
        
        % a few states possible
        tempstr         = cellfun(@(x,a,b) x(a:b),events,cac1,cac2,'UniformOutput',false);
        detectionNum    = get_detection_Num(tempstr);
        % get string event
        newstr          = get_newstr(tempstr);
        
        ld_detect.detectionStatus     = detectionNum;
        ld_detect.detectionText       = newstr;
        
        
        % GET the previous detection state
        % get type
        cac1    = cellfun(@(x) regexp(x, xprP),events,'UniformOutput',false);  %previous Detection state from above
        xpr     = 'LdDetectionEntry.Unused';
        cac2    = cellfun(@(x) regexp(x, xpr),events,'UniformOutput',false);
        
        % a few states possible
        tempstr       = cellfun(@(x,a,b) x(a:b),events,cac1,cac2,'UniformOutput',false);
        detectionNum  = get_detection_Num(tempstr);
        % get string event
        newstr        = get_newstr(tempstr);
        
        ld_detect.previousDetectionStatus = detectionNum;
        ld_detect.previousDetectionText   = newstr;
    else

        ld_detect = {'empty'};
    end
    

    varargout{1} = aDBS_log_tbl;
    varargout{2} = ld_detect;

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
    
    else
        
        varargout{1} = 'empty';
        varargout{2} = 'empty';
    end
else

    varargout{1} = 'empty';
    varargout{2} = 'empty';
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