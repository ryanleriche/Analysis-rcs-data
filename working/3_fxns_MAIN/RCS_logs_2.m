function INS_logs  = RCS_logs_2(cfg, pt_side_id)

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
    
        [AppLog_tbl.app{j}, AppLog_tbl.ld_detect{j}, AppLog_tbl.aDBS_state_def{j} ]...
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

    % aDBS state definitions
    i_tbl                    = cellfun(@(x) istable(x), AppLog_tbl.aDBS_state_def);
    aDBS_state_def               = unique(...
                                     vertcat(AppLog_tbl.aDBS_state_def{i_tbl}), 'rows');
    if istable(aDBS_state_def)
         aDBS_state_def.time.TimeZone = 'America/Los_Angeles';
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


        if istable(aDBS_state_def)
            tmp_aDBS_state_def  = unique([INS_logs.state_def; aDBS_state_changes], 'rows');
        else
            tmp_aDBS_state_def  = INS_logs.state_def;
        end
    end

else

    tmp_AppLog_tbl_path    = AppLog_tbl.path;
    tmp_EventLog_tbl_path  = EventLog_tbl.path;

    tmp_group_changes      = group_changes;
    tmp_aDBS_state         = aDBS_state_changes;
    tmp_aDBS_state_def     = aDBS_state_def;
    
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
    INS_logs.app  = tmp_aDBS_state;
    INS_logs.state_def           = tmp_aDBS_state_def;

    INS_logs.ld_detect           = tmp_ld_detect;
end
%% SAVE The Text Log structure
if ~isempty(AppLog_tbl) && ~isempty(EventLog_tbl)
    save(INS_log_dir,'INS_logs')
    fprintf('%s |.mat of INS Logs (as struct) saved to \n %s \n', pt_side_id, INS_log_dir);
    
end
end
