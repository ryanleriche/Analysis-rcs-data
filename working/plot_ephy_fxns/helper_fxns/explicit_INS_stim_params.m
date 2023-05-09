%function [redcap, group_changes]  = explicit_INS_stim_params(cfg, group_changes, db_RCSXXX, redcap, visits_tbl)


%{

from group_changes (as INS Group changes are written to EventLog.txts), 
find explicit stim paramemters (from StimLog.jsons which are NOT saved on 
the INS itself, but during Streaming Sessions)

%}



%%% uncomment below to troubleshoot as script:

cfg                    = [];
cfg.stage_dates        = stage_dates{2};
cfg.pt_id              = 'RCS02';

group_changes          = INS_logs.RCS02R.group_changes;
db_RCSXXX              = db.RCS02R;

visits_tbl             = visits.RCS02;
redcap                 = REDcap.RCS02;
%
%
%
proc_app = proc_app_log.RCS02R;

%%
% only takes REDcap surveys from AFTER Stage I implant
redcap     = redcap(ge(redcap.time, datetime(cfg.stage_dates(1),'TimeZone','America/Los_Angeles')),...
                  :); 

i_stimSett   = cellfun(@(x) ~isempty(x), db_RCSXXX.stimSettingsOut);

tmp_stim_tbl = db_RCSXXX.stimSettingsOut(i_stimSett);
tmp_stim_tbl = cellfun(@(x) x(end,:), tmp_stim_tbl, 'UniformOutput', false);



stim_tbl     = [db_RCSXXX(i_stimSett, {'sess_name', 'path'}),...
                vertcat(tmp_stim_tbl{:})];

i_entry    = find(~cellfun(@isempty, stim_tbl.GroupA));

Groups              = {'GroupA', 'GroupB', 'GroupC', 'GroupD'};


% replace empty Groups w/ previous non-empty group
for i=1:height(stim_tbl)

    if isempty(stim_tbl{i, 'GroupA'}{1})

        emp_diff = (i - i_entry);

        [~, i_prev] = min(emp_diff(emp_diff>0));

        stim_tbl(i,Groups) = stim_tbl(i_entry(i_prev), Groups);

    end
end



time_API   = datetime(stim_tbl.HostUnixTime /1000,...
                        'ConvertFrom','posixTime',...
                        'TimeZone','America/Los_Angeles',...
                         'Format','dd-MMM-yyyy HH:mm:ss.SSS');


stim_tbl    = addvars(stim_tbl, time_API ,'after',"sess_name");


[~, i_unique]          = unique(stim_tbl.time_API);
stim_tbl    = stim_tbl(i_unique,:);

stim_tbl    = sortrows(stim_tbl, 'time_API');


%
% organize INS logs and then infer contacts, amp, pw, and rate based off of DeviceSettings.txt
%group_changes   = group_changes(~contains(group_changes.event, "Lead"),:);
group_changes   = renamevars(group_changes, 'time', 'time_INS');

% remove instances of identical repeated events
i_rep_events      = diff([inf; findgroups(group_changes.event)])~=0;

group_changes   = group_changes(i_rep_events,:);



i_on_off  = find(contains(group_changes.event, {'Off', 'On'}));
i_group   = find(contains(group_changes.event, 'Group'));

% i_lead_int   = find(contains(group_changes.event, 'Lead'));
% unique(group_changes.event(i_lead_int+1))
tmp_groups       = group_changes.event(i_group, :);
tmp_ther_stat    = cell(length(i_group),1);

% based on nearest previous On/Off explicitly have Group X On or Off
for i = 1 : length(i_group)

    t_diff     = i_group(i) - i_on_off;

    near_t     = min(t_diff(t_diff > 0));
    j          = find(t_diff == near_t);

    tmp_ther_stat{i} ...
        = group_changes.event{i_on_off(j)};
             
end

% w/ previous therapyStatus assigned to Group, remove all therapyStatuses,
% and any repeats
tmp_concat = cellfun(@(x,y) [x,'_',y], tmp_groups, tmp_ther_stat, 'UniformOutput', false);

group_changes.event(i_group) = tmp_concat;

group_changes(i_on_off,:) = [];

i_rep_events    = diff([inf; findgroups(group_changes.event)])~=0;
group_changes   = group_changes(i_rep_events,:);


%% from previous StimLog.json get explicit stim params
for i = 1 : height(group_changes)

    t_diff   =  group_changes.time_INS(i) - stim_tbl.time_API;

    % find the smallest negative time difference to get nearest
    % stimLog.json settings BEFORE EventLog.txt settings

    near_t  = min(t_diff(t_diff > 0));

    % need to have streaming session occuring BEFORE INS log time
    if ~isempty(near_t)
        j       = find(t_diff == near_t);

        group_changes.prev_sess_name(i)         = stim_tbl.sess_name(j);
        group_changes.prev_sess_path(i)         = stim_tbl.path(j);

        % Since INS log can change groups SINCE last stimLog setting, find
        % explicit paramemters regardless of last activeGroup in StimLog
    
        if contains(group_changes.event{i}, 'Group')
            act_group                     = group_changes.event{i}(1:6);
        
            group                         = stim_tbl.(act_group){j(1)};
        
            h                             = find(group.validPrograms);
            
            if ~isempty(h)
                group_changes.stimContacts{i} = cell(length(h), 1);
        
                for k = 1 : length(h)
                    % this is monopolar stimulation
                    if group.contacts(h(k),:).anodes{1} == 16
        
                        group_changes.stimContacts{i}{k} = sprintf('c+%d-%d-%d-%d-%d-%d-',...,...
                                                           group.contacts(h(k),:).cathodes{1}); 
        
                    else 
        
                                               % only one cathode, but MANY
                                               % possible anodes
                    group_changes.stimContacts{i}{k} = sprintf('%d+%d-%d-%d-%d-%d-%d-',...
                                                           group.contacts(h(k),:).anodes{1},...
                                                           group.contacts(h(k),:).cathodes{1}); 
        
                    end
                end
        
        
                group_changes.ampInMilliamps{i}               = group.ampInMilliamps(h);
                group_changes.pulseWidthInMicroseconds{i}     = group.pulseWidthInMicroseconds(h);
                
                group_changes.rateInHz(i)         = group.rateInHz;
                
                group_changes.cycleOnInSecs(i)    = group.cycleOnInSecs;
                group_changes.cycleOffInSecs(i)   = group.cycleOffInSecs;
                group_changes.rampInSecs(i)       = group.rampInSecs;
                group_changes.rampRepeat(i)       = group.rampRepeat;
            
            
            else
        
                group_changes.stimContacts{i} = 'No valid program';
            
            end
        end
    end
end


%% adding in side + region for unambiguous contacts when comparing both sides
tmp_group_changes = group_changes;

switch cfg.pt_id
    
    case 'RCS02'
        stimRegR      = [{'RACC ', ["0","1","2","3"]}; {'RThal ', ["8","9","10","11"]}];
end

i_valid         = cellfun(@(x) ~isempty(x),...
                       tmp_group_changes.ampInMilliamps);

temp_cont       = cellfun(@(x) x, ...
                    tmp_group_changes.stimContacts,'UniformOutput', false);


ind_contacts    = cellfun(@(x) regexp(x,'\d*','Match'), ...
                               temp_cont(i_valid), 'UniformOutput', false);

ind_contacts    = cellfun(@(x) x{1}(1), ind_contacts);


i_small         = cellfun(@(x) any(strcmp(x, stimRegR{1,2})), ind_contacts);
i_large         = cellfun(@(x) any(strcmp(x, stimRegR{2,2})), ind_contacts);

i_con           = find(i_valid);


tmp_group_changes.stimContacts(i_con(i_small)) =...
    ...
    cellfun(@(x) [stimRegR{1,1}, x{1}], ...
    temp_cont(i_con(i_small)), 'UniformOutput', false);

tmp_group_changes.stimContacts(i_con(i_large)) =...
    ...
    cellfun(@(x) [stimRegR{2,1}, x{1}], ...
    temp_cont(i_con(i_large)), 'UniformOutput', false);

% reformat given that all programs are at 0 mA
for i = 1 : height(tmp_group_changes)
    if length(tmp_group_changes.ampInMilliamps{i}) > 1 && all(tmp_group_changes.ampInMilliamps{i} == 0)

        tmp_group_changes.ampInMilliamps{i} = 0;

        tmp_group_changes.pulseWidthInMicroseconds{i}...
            =...
        tmp_group_changes.pulseWidthInMicroseconds{i}(1);

    elseif length(tmp_group_changes.ampInMilliamps{i}) < 1

        tmp_group_changes.ampInMilliamps{i}             = NaN;
        tmp_group_changes.pulseWidthInMicroseconds{i}   = NaN;

        tmp_group_changes.rateInHz(i)                   = NaN;

        tmp_group_changes.cycleOnInSecs(i)              = NaN;
        tmp_group_changes.cycleOffInSecs(i)             = NaN;

        tmp_group_changes.rampInSecs(i)                 = NaN;
        tmp_group_changes.rampRepeat(i)                 = NaN;

    end
end

tmp_group_changes.ampInMilliamps             = vertcat(tmp_group_changes.ampInMilliamps{:});
tmp_group_changes.pulseWidthInMicroseconds   = vertcat(tmp_group_changes.pulseWidthInMicroseconds{:});
 
group_changes  = tmp_group_changes;


%% w/ explicit INS logs, now align to nearest REDcap

for i = 1 : height(redcap)

    t_diff   =  redcap.time(i) - group_changes.time_INS;
  
    near_t  = min(t_diff(t_diff > 0));


    % if no previous INSLog, use StimLog.json
    if isempty(near_t)
    % first INS Log started AFTER first month of REDcap -> use
    % DeviceSettings.json

        t_diff   = redcap.time(i) - stim_tbl.time_API;
        near_t   = min(t_diff(t_diff > 0));
        j        = find(t_diff == near_t);

        act_group       = stim_tbl.activeGroup{j};

        group           = stim_tbl.(['Group', act_group])(j);

        h               = find(group.validPrograms);
    
        redcap.therapyStatusDescription(i) = stim_tbl.therapyStatusDescription(j);
   
        if ~isempty(h)
             redcap.stimContacts{i} = cell(length(h), 1);
    
            for k = 1 : length(h)
                % this is monopolar stimulation
                if group.contacts(h(k),:).anodes{1} == 16
    
                    group_changes.stimContacts{i}{k} = sprintf('c+%d-%d-%d-%d-%d-%d-',...,...
                                                       group.contacts(h(k),:).cathodes{1}); 
    
                else 
    
                                           % only one cathode, but MANY
                                           % possible anodes
                group_changes.stimContacts{i}{k} = sprintf('%d+%d-%d-%d-%d-%d-%d-',...
                                                       group.contacts(h(k),:).anodes{1},...
                                                       group.contacts(h(k),:).cathodes{1}); 
    
                end
            end
    
    
             redcap.ampInMilliamps{i}               = group.ampInMilliamps(h);
             redcap.pulseWidthInMicroseconds{i}     = group.pulseWidthInMicroseconds(h);
            
             redcap.rateInHz(i)         = group.rateInHz;
            
             redcap.cycleOnInSecs(i)    = group.cycleOnInSecs;
             redcap.cycleOffInSecs(i)   = group.cycleOffInSecs;
             redcap.rampInSecs(i)       = group.rampInSecs;
             redcap.rampRepeat(i)       = group.rampRepeat;
            
        else

            redcap.activeGroup{i}     = [];
            redcap.stimContacts(i)    = {'   '};
    
            redcap.ampInMilliamps{i}         = NaN;
            redcap.pulseWidthInMicroseconds{i}          = NaN;
            
            redcap.rateInHz(i)        = NaN;
            
            redcap.cycleOnInSecs(i)    = NaN;
            redcap.cycleOffInSecs(i)   = NaN;
    
            redcap.rampInSecs(i)       = NaN;
            redcap.rampRepeat(i)     = NaN;
        
        
        end
    

    else

       t_diff   =  redcap.time(i) - group_changes.time_INS;
       near_t   = min(t_diff(t_diff > 0));
        j       = find(t_diff == near_t);
    
        vars = group_changes(j, :).Properties.VariableNames;
    
        for h = 1 : length(vars)
            redcap.(vars{h})(i) = group_changes.(vars{h})(j(1));
    
        end
    end

    % see if REDcap report occured during inpatient, inclinic, or home testing
    visit_sess_diff      =  visits_tbl.dates - redcap.time(i);
    
    i_visit_day    = find(le(visit_sess_diff, duration('0:00:00')) &...
                  ge(visit_sess_diff, '-24:00:00'));
        
    if ~isempty(i_visit_day)
    
        redcap.visits(i)  = visits_tbl.desc(i_visit_day);
    else
    
        redcap.visits(i)  = {'   '};
    end   
end

redcap                = movevars(redcap,'time_INS', 'After','MPQcruel');



%%%% uncomment below to troubleshoot as script

REDcap_INSLog.RCS02R        = redcap;
proc_group_changes.RCS02R   = group_changes;
%end