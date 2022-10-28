% THIS WILL RECOMPILE:
% For a specific patient named in cfg.pt_id
% 
% 1. Redcap pain scores 
% 2. Makedatabase for JSON files 
% 3. RCS text logs 
% 
% 
% PS 2022 



%% Compile REDCAP scores 
%  REDCAP API token for Ryan


REDcap                 = RCS_redcap_painscores(API_token);

save(fullfile(processed_dir,'redcap_painscores.mat'), 'REDcap');

%% COMPILE json databases 
% save RCS session summaries as DataBases w/n 'db' structure and
% badsessions w/n 'bs' structure 
cfg.load_EventLog      = false;
cfg.ignoreold          = false;


try 
[db.([cfg.pt_id 'L']), bs.([cfg.pt_id 'L'])] = ...
    makeDatabaseRCS_Ryan(data_raw_dir,[cfg.pt_id 'L'],cfg);
fprintf(' - OK compiled database:  %s \n',[cfg.pt_id 'L'])
catch ME_L
    fprintf(' - NO database compiled for %s  - check variable "ME_L"  \n',[cfg.pt_id 'L']);
end

try
[db.([cfg.pt_id 'R']), bs.([cfg.pt_id 'R'])] = ...
    makeDatabaseRCS_Ryan(data_raw_dir, [cfg.pt_id 'R'],cfg);
fprintf(' - OK Compiled database :  %s \n',[cfg.pt_id 'R'])
catch ME_R
    fprintf(' - NO database compiled for %s check variable "ME_R"  \n',[cfg.pt_id 'R']);
end


%% COMPILE text logs 
cfg.pull_adpt_logs      = false;
cfg.pull_event_logs     = true;
cfg.ld_detection_events = false;
cfg.pull_recharge_sess  = false;
cfg.pull_mirror_logs    = false;
cfg.savedir = processed_dir;

try 
[textlog.([cfg.pt_id 'L'])]         = RCS_logs(data_raw_dir,[cfg.pt_id 'L'],cfg);
fprintf(' - OK Compiled textlogs  :  %s \n',[cfg.pt_id 'L'])
catch TR_L
    fprintf(' - NO textlogs for %s \n',[cfg.pt_id 'L'])
end

try 
[textlog.([cfg.pt_id 'R'])]         = RCS_logs(data_raw_dir,[cfg.pt_id 'R'],cfg);
fprintf(' - OK Compiled textlogs  :  %s \n',[cfg.pt_id 'R'])
catch TR_R
    fprintf(' - NO textlogs for %s \n',[cfg.pt_id 'R'])
end


