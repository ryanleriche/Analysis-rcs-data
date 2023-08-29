function INS_logs  = RCS_logs(cfg, pt_side_id)

%% RCS INS Log Database
%{
This function loops through all session folders in the path to create a
database of all RC+S metrics of interest including processing of the
adaptive text logs

INPUTS: 
1. cfg.raw_dir is the local root directory pathname for all patient session files
      ! ! Make sure patient name is NOT included (e.g RCS02 would be a subfolder in the input folder)
       This should look something like 'C:/Desktop/'
           ** including:
           ** AppLog.txt: adaptive state changes
           ** EventLog.txt: open loop group changes

2. PATIENTIDside (this should indicate the side (L/R) of which device you are
      analyzing (i.e. RCS02R) - EXCEPT for CPRCS01, there is no letter after
      the name

OUTPUT: 
1. textlog.mat timetable with fields:

   [{'time'}; {'rec'}; {'sessname'  };  {'duration'  }; ...
    {'battery'   };{'TDfs'      };{'fft'};
    {'power'     };{'stim'    };   {'stimName'  }; ...
    {'stimparams'};  {'path'};  {'powerbands'}]; ...
    {'adaptiveLD_mean'}; {'adaptiveLD_std'}; {'untsreamedGroupChanges'};


function dependencies:
   read_adaptive_txt_log (starr lab analysis repo) % by Ro'ee Gilron
   makeDataBaseRCSdata (Chronic Pain RCS repo) % by Prasad Shirvalkar
      Open Mind functions in makeDataBaseRCSData: % by Kristin Sellers
      and Open Mind team
          	createDeviceSettingsTable
              createStimSettingsFromDeviceSettings
              deserializeJSON
              createStimSettingsTable
other dependencies:
https://github.com/JimHokanson/turtle_json
in the a folder called "toolboxes" in the same directory as the processing scripts

Ashlyn Schmitgen May2021

%% Updates %%%
  P Shirvalkar July 28 2021 
  -  Updated redundant calculations and included recharge sessions, group
     changes and detector changes to output
  - created new plotting tools for this analysis (as part of RCS_CL which
  calls this function)
%%

For OpenMind
%}



%% just RCS02R for development
%{
'every_entry' call works for RCS02

* include logic to only include new textlogs
    * load the paths -> load previous textlog_db
    * only load paths that do NOT exist in textlog_db
    

NOT Done:

* same fieldnames for 'every_entry' versus 'blazing_fast' calls
%}
% cfg                    = [];
% cfg.raw_dir            = pia_raw_dir;
% cfg.pt_id              = 'RCS02R';
% cfg.ignoreold_INS_logs          = false;

%% create loop through all text files for adaptive_read_log_txt.m for database

%warning("off", "all"); 
tic
cfg.proc_dir = fullfile(cfg.proc_dir, 'INS_logs/');

if ~exist(cfg.proc_dir, 'dir');     mkdir(cfg.proc_dir);      end


% exception for CPRCS01
if ~ (contains(pt_side_id,'CPRCS01'))
    pt_id = pt_side_id(1:end-1); %remove the L or R letter
else 
    pt_id = pt_side_id;
end

fprintf('%s | compiling INS logs \n', pt_side_id)
    
scbs_dir     = fullfile(cfg.raw_dir, pt_id,'SummitData/SummitContinuousBilateralStreaming', pt_side_id);
adbs_dir     = fullfile(cfg.raw_dir, pt_id,'SummitData/StarrLab', pt_side_id);

filelist     = [dir(fullfile(scbs_dir,'**/*.txt')); dir(fullfile(adbs_dir, '**/*.txt'))]; % all txt files contains within session files
% remove the files that start with ._  (some icloud issue of duplicate files to ignore)
badfiles     = arrayfun(@(x) contains(x.name,'._'),filelist);
filelist(badfiles)=[];
filelist      = struct2table(filelist(~[filelist.isdir]));
filelist.date = datetime(filelist.date);

filelist      = sortrows(filelist,'date');

path          = cellfun(@(x, y) [x,'/', y], filelist.folder, filelist.name,...
                             'UniformOutput', false);

i_event                     = endsWith(path, 'EventLog.txt'); 
EventLog_tbl                = table();
EventLog_tbl.path           = path(i_event);
EventLog_tbl.group_changes  = cell(sum(i_event),1);
EventLog_tbl.rech_sess      = cell(sum(i_event),1);

% initiate cell array of every txt log
AppLog_tbl              = table();
i_app                   = endsWith(path, 'AppLog.txt'); 
AppLog_tbl.path         = path(i_app);

% initiate cell array of every app log
AppLog_tbl.app   = cell(sum(i_app),1);
AppLog_tbl.adapt_stat   = cell(sum(i_app),1);
AppLog_tbl.ld_detect    = cell(sum(i_app),1);

% if needed, make folder to save processed INS logs
if ~exist(cfg.proc_dir, 'dir');     mkdir(cfg.proc_dir);    end

INS_log_dir  = sprintf('%s%s_INS_logs.mat', cfg.proc_dir, pt_side_id);

%%
% save paths in final output as reference so logs are not redundantly ran
% % when new text logs are generated
% INS_logs.AppLog_tbl_path              = AppLog_tbl.path(1 : 10);
% INS_logs.EventLog_tbl_path            = EventLog_tbl.path(1 : 10);


if cfg.ignoreold_INS_logs == false && isfile(INS_log_dir) % handles edge case of no file existing

    fprintf('%s | loading existing logs \n', pt_side_id)
    
    
    % brings in 'INS_logs' struct
    load(INS_log_dir, 'INS_logs');

    raw_applogs  = AppLog_tbl.path;

    i_raw_dirs   = ~contains(raw_applogs, INS_logs.AppLog_tbl_path);
    AppLog_tbl   = AppLog_tbl(i_raw_dirs, :);
 
    % repreat w/ event logs
    raw_eventlogs   = EventLog_tbl.path;
    i_raw_dirs     = ~contains(raw_eventlogs, INS_logs.EventLog_tbl_path);

    EventLog_tbl   = EventLog_tbl(i_raw_dirs, :);

end
% EventLog_tbl   = EventLog_tbl(1:10,:);
% AppLog_tbl     = AppLog_tbl(1:10,:);


%% parse through AppLog.txt files
if ~isempty(AppLog_tbl)
    for j = 1: height(AppLog_tbl)
    
        fn = AppLog_tbl.path{j};
    
        [AppLog_tbl.app{j}, AppLog_tbl.ld_detect{j}, ~ ]...
            = ...
        read_INS_logs_fast(fn);
        
    end
else
    fprintf('%s | reportedly no new AppLog.txt files to parse\n', pt_side_id)

end

%% parse through EventLog.txt files
if ~isempty(EventLog_tbl)
    for j = 1 : height(EventLog_tbl)
    
        fn  = EventLog_tbl.path{j};
    
        [EventLog_tbl.group_changes{j}, EventLog_tbl.rech_sess{j}]...
            = ...
        read_INS_logs_fast(fn);
    end
else

    fprintf('%s | reportedly no new EventLog.txt files to parse\n', pt_side_id)
end

% EventLog_tbl    = EventLog_tbl(cellfun(@(x) ~isempty(x), EventLog_tbl.group_changes),:);
% AppLog_tbl      = AppLog_tbl(cellfun(@(x) ~isempty(x), AppLog_tbl.app),:);
%% organize EventLogs
if ~isempty(EventLog_tbl)
    % group changes
    i_tbl                       = cellfun(@(x) istable(x), EventLog_tbl.group_changes);
    group_changes               = unique(...
                                         vertcat(EventLog_tbl.group_changes{i_tbl}), 'rows');
    
    % timezone to match rest of RCS data
    group_changes.time.TimeZone = 'America/Los_Angeles';
    
    % recharge sessions
    i_tbl                       = cellfun(@(x) istable(x), EventLog_tbl.rech_sess);
    rech_sess                   = unique(...
                                         vertcat(EventLog_tbl.rech_sess{i_tbl}), 'rows');
    
    % timezone to match rest of RCS data
    rech_sess.time.TimeZone = 'America/Los_Angeles';
end

if ~isempty(AppLog_tbl)
    %%% organize AppLogs
    i_tbl                    = cellfun(@(x) istable(x), AppLog_tbl.app);
    aDBS_state_changes               = unique(...
                                     vertcat(AppLog_tbl.app{i_tbl}), 'rows');
    if istable(aDBS_state_changes)
        aDBS_state_changes.time.TimeZone = 'America/Los_Angeles';
    end
    % LD Detections
    i_tbl                    = cellfun(@(x) istable(x), AppLog_tbl.ld_detect);
    ld_detect                = unique(...
                                     vertcat(AppLog_tbl.ld_detect{i_tbl}), 'rows');
    
    if istable(ld_detect)
        ld_detect.time.TimeZone = 'America/Los_Angeles';
    end

end
%%
if cfg.ignoreold_INS_logs == false && isfile(INS_log_dir) % handles edge case of no file existing
    % from already processed INS logs, include only new entries (subsequent
    % INS logs often contain overlapping entries)

    if ~isempty(EventLog_tbl)
        tmp_EventLog_tbl_path  = unique([INS_logs.EventLog_tbl_path; EventLog_tbl.path]);
        tmp_group_changes      = unique([INS_logs.group_changes; group_changes], 'rows');
        tmp_rech_sess          = unique([INS_logs.recharge; rech_sess], 'rows');


    end

    if ~isempty(AppLog_tbl)

        tmp_AppLog_tbl_path    = unique([INS_logs.AppLog_tbl_path; AppLog_tbl.path]);


        if istable(aDBS_state_changes)
            tmp_aDBS_state         = unique([INS_logs.app; aDBS_state_changes], 'rows');
        else
            tmp_aDBS_state         = INS_logs.app;
        end

        tmp_ld_detect          = unique([INS_logs.adaptive; ld_detect], 'rows');


   
    end

else

    tmp_AppLog_tbl_path    = AppLog_tbl.path;
    tmp_EventLog_tbl_path  = EventLog_tbl.path;

    tmp_group_changes      = group_changes;
    tmp_aDBS_state         = aDBS_state_changes;
  
    
    tmp_ld_detect          = ld_detect;
    tmp_rech_sess          = rech_sess; 



end

if ~isempty(EventLog_tbl)
    INS_logs.EventLog_tbl_path  = tmp_EventLog_tbl_path;
    INS_logs.recharge           = tmp_rech_sess; 

    INS_logs.group_changes      = tmp_group_changes;

end

if ~isempty(AppLog_tbl)
    INS_logs.AppLog_tbl_path     = tmp_AppLog_tbl_path;
    INS_logs.app                 = tmp_aDBS_state;
   
    INS_logs.ld_detect           = tmp_ld_detect;
end
%% SAVE The Text Log structure
if ~isempty(AppLog_tbl) && ~isempty(EventLog_tbl)

    save(INS_log_dir,'INS_logs')
    fprintf('%s |.mat of INS Logs (as struct) saved to \n %s \n', pt_side_id, INS_log_dir);
    
end

%% scratch code from slow parsing
%{
proc_dirname      = [cfg.raw_dir(1:end-4), 'processed/INS_logs/'];

outputFileName    = fullfile(proc_dirname,[cfg.pt_id, '_INS_logs.mat']);

if isfile(outputFileName) && ~cfg.ignoreold_INS_logs

    disp('Loading previous INS Logs');
    old = load(outputFileName);

    i_proc_logs   = find(cellfun(@(x) ~isempty(x), old.INSLog_tbl.events));
    proc_INS_logs = old.INSLog_tbl(i_proc_logs,:); 

    filelist(i_proc_logs, :) = [];


end


INSLog_tbl  = table;

INSLog_tbl.folder = filelist.folder;
INSLog_tbl.name   = filelist.name;
INSLog_tbl.date   = filelist.date;
INSLog_tbl.events = cell(height(INSLog_tbl),1);

% loop through all EventLog.txts, and AppLog.txts
for i = 251 : height(filelist)

    fn           = [filelist.folder{i}, '/', filelist.name{i}];

    if endsWith(fn, 'EventLog.txt') 

        temp_events  = read_INS_logs_slow(fn);

    
        if ~iscell(temp_events) && ~isempty(temp_events)

            INSLog_tbl.events{i}  = temp_events;
        else
            INSLog_tbl.events{i}  = 'empty';
        end
    end
end


if isfile(outputFileName) && ~cfg.ignoreold_INS_logs

    i_new_logs      = find(cellfun(@(x) ~isempty(x), INSLog_tbl.events));

    new_INS_logs    = INSLog_tbl(i_new_logs, :);

    temp_INS_logs   = [proc_INS_logs; new_INS_logs];

    temp_INS_logs.folder = cellfun(@(x, y) [x, y],...
                          temp_INS_logs.folder, temp_INS_logs.name,...
                          'UniformOutput', false);

    [~, i_u]        =  unique(temp_INS_logs.folder);

    if length(i_u) ~= height(temp_INS_logs)

        disp('INS logs redundantly loaded in--may be erroneous')
    end

    proc_folders     = temp_INS_logs.folders;


end
%%

save([proc_dirname, cfg.pt_id, '_proc_INS_logs.mat'],...
     'proc_folders', '-v7.3');
%}
% organize EventLogs
%{
EventLogs            = INSLogs(endsWith(INSLogs.name,"EventLog.txt"),:);

i_tbl                = cellfun(@(x) istable(x), EventLogs.events);

events               = vertcat(EventLogs.events{i_tbl});


% TherapyStatus tbl
TherStat      = EventLogs(strcmp(EventLogs.event_id,'TherapyStatus'), : );

exp_entries   = struct2table(cellfun(@(x) x, TherStat.entries));

TherStat      = removevars(TherStat, {'event_id', 'entry_names', 'entries'});
TherStat      = [TherStat, exp_entries];

% ActiveDeviceChanged tbl
ActDev        = EventLogs(strcmp(EventLogs.event_id, 'ActiveDeviceChanged'), : );

exp_entries   = struct2table(cellfun(@(x) x, ActDev.entries));

ActDev        = removevars(ActDev , {'event_id', 'entry_names', 'entries'});
ActDev        = [ActDev, exp_entries];

% rename Groups names to letters
ActDev.NewGroup = replace(ActDev.NewGroup,...
                    {'Group0','Group1','Group2','Group3'}, ...
                    {'GroupA','GroupB','GroupC','GroupD'});


% combine TherapyStatus (On and Off) and ActiveDeviceChanged (Groups) events
temp_tbl = removevars(ActDev, 'Unused');
temp_tbl = renamevars(temp_tbl, 'NewGroup', 'TherapyStatus');

GroupChange_tbl = sortrows( ...
                           [temp_tbl; ...
                           removevars(TherStat, {'TherapyStatusType','Unused'})],...
                    'time');

% timezone to match rest of RCS data
GroupChange_tbl.time.TimeZone = 'America/Los_Angeles';

% remove duplicate rows (time can be the same if there's subsequent events
% w/n a second)
GroupChange_tbl  = unique(GroupChange_tbl, 'rows');



% organize RechargeSess[i]ons --

RechSess       = EventLogs(strcmp(EventLogs.event_id,'RechargeSesson'), : );

exp_entries    = struct2table(cellfun(@(x) x, RechSess.entries));

RechSess       = removevars(RechSess, {'event_id', 'entry_names', 'entries'});
RechSess       = [RechSess, exp_entries];

RechSess       = removevars(unique(RechSess, 'rows'), 'Unused');
%}
% organize AppLogs
%{
AppLogs            = INSLogs(endsWith(INSLogs.name,"AppLog.txt"), 'events');

AppLogs            = vertcat(AppLogs.events{:});

% AdaptiveTherapyStateChange as its own table
AdTherStateChange            = AppLogs(strcmp(AppLogs.event_id, 'AdaptiveTherapyStateChange'), : );


exp_entries         = struct2table(cellfun(@(x) x, AdTherStateChange.entries));

AdTherStateChange            = removevars(AdTherStateChange, {'event_id', 'entry_names', 'entries'});

AdTherStateChange            = [AdTherStateChange, exp_entries];

% remove (redundant) variable names saved in hexidecimal format
var                 = {'Status', 'Prog0Amp', 'Prog1Amp', 'Prog2Amp', 'Prog3Amp',...
                       'RatePeriodAtTimeOfModification'};

AdTherStateChange            = removevars(AdTherStateChange, var);


var                 = ["Prog0AmpInMillamps", "Prog1AmpInMillamps",...
                       "Prog2AmpInMillamps", "Prog3AmpInMillamps",...
                       "RateAtTimeOfModification"];

for i = 1 : length(var)
    AdTherStateChange.(var{i}) = cellfun(@(x) round(str2double(x),2), AdTherStateChange.(var{i}));    
end

% timezone to match rest of RCS data
AdTherStateChange.time.TimeZone = 'America/Los_Angeles';

% remove duplicate rows (time can be the same if there's subsequent events
% w/n a second)
AdTherStateChange              = unique(AdTherStateChange, 'rows');



LdDetEvent      = AppLogs(strcmp(AppLogs.event_id, 'LdDetectionEvent'), : );

exp_entries     = struct2table(cellfun(@(x) x, LdDetEvent.entries));

LdDetEvent      = removevars(LdDetEvent  , {'event_id', 'entry_names', 'entries'});
LdDetEvent      = [LdDetEvent  , exp_entries]; 

LdDetEvent      = unique(...
                    removevars(LdDetEvent, 'Unused'), ...
                  'rows');

% AdTherStateWrite as table -> RBL:  seems redundant w/ 'AdaptiveTherapyChange' event
%{
AdTherStateWrite        = AppLogs(strcmp(AppLogs.event_id, 'AdaptiveTherapyStateWritten'), : );


exp_entries             = struct2table(cellfun(@(x) x, AdTherStateWrite.entries));

AdTherStateWrite        = removevars(AdTherStateWrite, {'event_id', 'entry_names', 'entries'});

AdTherStateWrite        = [AdTherStateWrite, exp_entries];




% timezone to match rest of RCS data
AdTherStateWrite.time.TimeZone = 'America/Los_Angeles';

% remove duplicate rows (time can be the same if there's subsequent events
% w/n a second)
AdTherStateWrite              = unique(AdTherStateWrite, 'rows');
%}



% LfpSenseStateEvent as table -> RBL: just data that is sensed
% 'LDDetectionEvent' tells us what state the INS is in (determining stim)
%{
LfpSSEvent      = AppLogs(strcmp(AppLogs.event_id, 'LfpSenseStateEvent'), : );

j               = cellfun(@(x) ~isempty(x), LfpSSEvent.entries);

LfpSSEvent      = LfpSSEvent(j, :);
exp_entries     = struct2table(cellfun(@(x) x, LfpSSEvent.entries));

LfpSSEvent      = removevars(LfpSSEvent, {'event_id', 'entry_names', 'entries'});
LfpSSEvent      = [LfpSSEvent, exp_entries]; 

LfpSSEvent      = unique(...
                    removevars(LfpSSEvent, 'Unused'), ...
                  'rows');

%}

%%


org_INS_logs.event_logs      = EventLogs;
org_INS_logs.ther_stat       = TherStat;
org_INS_logs.act_dev         = ActDev;
org_INS_logs.group_changes   = GroupChange_tbl;

org_INS_logs.app      = AdTherStateChange;
org_INS_logs.ld_det          = LdDetEvent;

org_INS_logs.rech            = RechSess;
%}
end
