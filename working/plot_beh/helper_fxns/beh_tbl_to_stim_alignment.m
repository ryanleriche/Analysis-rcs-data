% function [stim_log, redcap] ...
%     ...
%     = beh_tbl_to_stim_alignment(...
%     ...
%     cfg, pt_side_id, db, INS_logs_API_t_synced)
%%
pt_side_id = pt_sides{i};

%% index db and time aligned INS log for pt hemisphere    
db_RCSXXX     = db.(pt_side_id);
proc_g_chan   = INS_logs_API_t_synced.(pt_side_id).group_changes;
proc_app      = INS_logs_API_t_synced.(pt_side_id).app;

%%


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

stim_log =   movevars(...
                     vertcat(db_RCSXXX.stimLogSettings{i_stimlog}),...
            'time_stimLog', 'Before', 'HostUnixTime');

n_log_per_sess   = cellfun(@(x) height(x), {db_RCSXXX.stimLogSettings{i_stimlog}})';
sess_w_logs      = db_RCSXXX.sess_name(i_stimlog);
   

rep_sess_names = cellfun(@(x,y) ...
                        repmat({x}, y, 1), ...
                        sess_w_logs, num2cell(n_log_per_sess), ...
                        'UniformOutput', false);

rep_sess_names      = vertcat(rep_sess_names{:});
stim_log.sess_name  = rep_sess_names;

% vars    = {'time_stimLog', 'activeGroup', 'therapyStatus', 'therapyStatusDescription',...
%            'stimParams_prog1', 'stimParams_prog2', 'stimParams_prog3', 'stimParams_prog4',...
%            'sess_name'};
% 
% % [~, i_u]  = unique(stim_log(:, vars), 'rows');
% % stim_log   = stim_log(i_u, :);

% As of Nov. 2022, we have not intentionally used multiple programs w/n a
% group, display all the unique stim settings of Programs 1, 2, 3

% **note** reporting using 0->3 (like Medtronic) for programs, but keeping 1->4 variable
% naming for back back compatibility w/ Analysis-rca-data
disp(strjoin(...
        [pt_side_id; '| Unique Program 1 Settings:', newline; char(131); unique(stim_log.stimParams_prog2)]...
            ));

disp(strjoin(...
        [pt_side_id; '| Unique Program 2 Settings:', newline; char(131); unique(stim_log.stimParams_prog3)]...
            ));


disp(strjoin(...
        [pt_side_id; '| Unique Program 3 Settings:', newline; char(131); unique(stim_log.stimParams_prog4)]...
            ));

% remove redundant fields
stim_log = removevars(stim_log, ...
                      {'therapyStatus','HostUnixTime','updatedParameters',...
                      'stimParams_prog2', 'stimParams_prog3', 'stimParams_prog4'});

stim_log = movevars(stim_log, 'sess_name', 'After', 'time_stimLog');
%% include merged INS + StimLog to account for PTM initiated changes
var_names      = stim_log.Properties.VariableNames;

add_var_to_INS = var_names(...
                        ~contains(var_names, ...
                        {'time', 'sess_name', 'activeGroup', 'therapyStatusDescription'}));

% remove entries occuring before first stim_log entry and Lead Integrity 
% (i.e., impedance) tests (impedance test returns to previous Group once
% done)
i_rmv = le(proc_g_chan.time_INS, stim_log.time_stimLog(1)) |...
            strcmp(proc_g_chan.event, 'LeadIntegrityTest');

proc_g_chan(i_rmv, :) = [];

stim_log_b4     = cell(height(stim_log),1);
INS_entr_btwn   = stim_log_b4;

for j   =  1 : height(stim_log) - 1  
    i_btwn = find(...
                    ge(proc_g_chan.time_INS, stim_log.time_stimLog(j)) &...
                    le(proc_g_chan.time_INS, stim_log.time_stimLog(j+1) - duration('0:01:00'))...
                    );
    
    INS_entr_btwn{j} = proc_g_chan(i_btwn, {'time_INS', 'event'});
    stim_log_b4{j}   = repmat(stim_log(j, add_var_to_INS), length(i_btwn), 1);

    %%%infer INS log entries based on inital stim_log entry


end

i_emp  = cellfun(@isempty, INS_entr_btwn);
i_entr = find(~i_emp);

i_PTM = cellfun(@(x) ...
            sum(ge(diff(x.time_INS), duration('0:01:00'))) >1, INS_entr_btwn(i_entr));


INS_oi = INS_entr_btwn(i_entr(i_PTM))

INS_entr_btwn(i_emp(i_PTM ))




i_INS_2_add     = vertcat(stim_log.ind_INS_entries{:});



tmp_INS_tbl = [proc_g_chan(i_INS_2_add, {'time_INS', 'event'}), ...
                  vertcat(stim_log_oi{:})];


%% from DeviceSettings.json add in cycling, ramping, and Active Recharging

stimSettingsOut  = vertcat(db_RCSXXX.stimSettingsOut{i_devset});

i_empty          = cellfun(@isempty, [stimSettingsOut.GroupA]);

stimSettingsOut  = stimSettingsOut(~i_empty, :);
% use of Program 1 is fine assumption given lack of meaningful stim params
% in Programs 2, 3, and 4

for i = 1 : height(stim_log)

    stim_para = strsplit(char(stim_log.stimParams_prog1(i)),',');
    
    if ~strcmp(stim_para, 'Disabled')

        if isempty(stim_para{1})
            stim_log.stimContacts{i}  = '';
        else
            stim_log.stimContacts{i}  = stim_para{1};
        end
    
        stim_log.ampInMilliamps(i)              = str2double(stim_para{2}(2:end -2));
        stim_log.pulseWidthInMicroseconds(i)    = str2double(stim_para{3}(2:end -2));
        stim_log.rateInHz(i)                    = str2double(stim_para{4}(2:end -2));
    
        % from nearest, previous 'stimSettingsOut' 
    
        t_diff        = stim_log.time_stimLog(i) - stimSettingsOut.time_devset;
        near_t        = min(t_diff(t_diff >= 0));
    
        j             = find(t_diff == near_t);

        abcd          = stim_log.activeGroup{i};

        
    
        ss_group      = stimSettingsOut.(['Group' abcd]){j};
    
    
        stim_log.cycleEnabled(i)     = ss_group.cyclingEnabled;
        stim_log.cycleOnInSecs(i)    = ss_group.cycleOnInSecs;
        stim_log.cycleOffInSecs(i)   = ss_group.cycleOffInSecs;

        stim_log.percentDutyCycle(i) =  100* stim_log.cycleOnInSecs(i) ./...
                                          (stim_log.cycleOnInSecs(i) + stim_log.cycleOffInSecs(i));


        stim_log.rampInSecs(i)      = ss_group.rampInSecs;
        stim_log.rampRepeat(i)      = ss_group.rampRepeat;
    
        stim_log.achRechRatio(i)    = ss_group.actRechRatio(1);

        stim_log.cl_stim(i)         = strcmp(abcd, 'D'); 
    
    end  
end
%%
switch pt_side_id
    case 'RCS02R';  stimReg    = [{'RACC ', ["0","1","2","3"]}; {'RThal ', ["8","9","10","11"]}];
             
    case 'RCS04L';  stimReg    = [{'LACC ', ["0","1","2","3"]}; {'LCaud ', ["8","9","10","11"]}];
    case 'RCS04R';  stimReg    = [{'RACC ', ["0","1","2","3"]}; {'RThal ', ["8","9","10","11"]}];
        
        
    case 'RCS05L';  stimReg    = [{'LCaud ', ["0","1","2","3"]}; {'LACC ', ["8","9","10","11"]}];
    case 'RCS05R';  stimReg    = [{'RThal ', ["0","1","2","3"]}; {'RIFG ', ["8","9","10","11"]}];

        
    case 'RCS06L';  stimReg    = [{'LACC ',  ["0","1","2","3"]}; {'LCaud ', ["8","9","10","11"]}];
    case 'RCS06R';  stimReg    = [{'RThal ', ["0","1","2","3"]}; {'RSFG ', ["8","9","10","11"]}];
     
    % SGC (not ACC) was bilaterally implanted (RCS data is wrong)
    case 'RCS07L';  stimReg    = [{'LGPi ',  ["0","1","2","3"]}; {'LSGC ', ["8","9","10","11"]}];
    case 'RCS07R';  stimReg    = [{'RThal ', ["0","1","2","3"]}; {'RSGC ', ["8","9","10","11"]}];    
end

% adding in side + region for unambiguous contacts when comparing both sides
tmp          = stim_log.stimContacts;
i_con        = cellfun(@(x) ~isempty(x) , tmp);

ind_contacts = cellfun(@(x) regexp(x,'\d*','Match'), tmp(i_con), 'UniformOutput', false);
ind_contacts = cellfun(@(x) x(1), ind_contacts);

i_small      = cellfun(@(x) any(strcmp(x, stimReg{1,2})), ind_contacts);
i_large      = cellfun(@(x) any(strcmp(x, stimReg{2,2})), ind_contacts);

i_con        = find(i_con);


tmp(i_con(i_small)) = cellfun(@(x) [stimReg{1,1}, x], ...
                                    tmp(i_con(i_small)), 'UniformOutput', false);

tmp(i_con(i_large)) = cellfun(@(x) [stimReg{2,1}, x], ...
                                    tmp(i_con(i_large)), 'UniformOutput', false);

stim_log.stimContacts = tmp;

%%
for i = 1 : height(stim_log)
        % find REDcap setting that occurs AFTER given stimLog entry, but BEFORE
        % the subsequent entry 
        if i < height(stim_log)
            btwn_this_and_next = ...
                find(ge(redcap.time, stim_log.time_stimLog(i)) & ...
                le(redcap.time, stim_log.time_stimLog(i + 1)));
        
            if ~isempty(btwn_this_and_next)
            
                stim_log.redcap_btwn_stimLogs{i}     = 'btwn_stimLogs';
                stim_log.redcap_reports{i}           = redcap(btwn_this_and_next,:);
                stim_log.i_redcap(i)                 = {btwn_this_and_next};
    
        
            else 
        
                stim_log.redcap_btwn_stimLogs{i}     = 'none_btwn_stimLogs';
                stim_log.i_redcap(i)                 = {NaN};
            end
           
        % for last stimLog associate all reports from then to NOW rather than
        % next stimLog
        else
    
            btwn_this_and_next = ...
                find(ge(redcap.time, stim_log.time_stimLog(i)) & ...
                le(redcap.time, datetime('now','TimeZone','America/Los_Angeles')));
        
            if ~isempty(btwn_this_and_next)
            
                stim_log.redcap_btwn_stimLogs{i}     = 'btwn_stimLogs_and_now';
                stim_log.i_redcap(i)                 = {btwn_this_and_next};
        
            else 
        
                stim_log.redcap_btwn_stimLogs{i}     = 'none_btwn_stimLogs_and_now';
            end
        end
end


% remove redundant fields
stim_log = removevars(stim_log, {'stimParams_prog1'});
redcap   = removevars(redcap, {'stimParams_prog1'});


% verify REDcap to StimLog assignment
all_i_redcap      = vertcat(stim_log.i_redcap{:});

if length(all_i_redcap(~isnan(all_i_redcap))) == length(unique(all_i_redcap(~isnan(all_i_redcap)))) 
    
    fprintf('%s | All REDcap report(s) assigned to unique stim settings\n', pt_side_id)
else
    error('REDcap report(s) are assigned to multiple stim settings (RBL message).')
end

per_assigned = length(unique(all_i_redcap(~isnan(all_i_redcap)))) ./ height(redcap) * 100;


fprintf('%s | %.1f%% of REDcap report(s) assigned to stim settings\n', pt_side_id, per_assigned);


save_dir = [cfg.proc_dir,  '/stimLog/'];
if ~isfolder(save_dir);    mkdir(save_dir);   end

save(...
     [save_dir, pt_side_id, '_stimLog.mat'], ...
     ...
     'stim_log', 'redcap', '-v7.3');
%end
