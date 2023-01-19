%function beh_stim  = get_INSLog_stim_params(cfg, group_changes, db_RCSXXX, redcap, visits_tbl)


%{
from group_changes (as written to EventLog.txts), find explicit stim paramemters
(from StimLog.jsons)

%}
cfg                    = [];
cfg.stage_dates        = stage_dates{2};
cfg.pt_id              = 'RCS02';

group_changes          = INS_logs.RCS02R.group_changes;
db_RCSXXX              = db.RCS02R;
redcap                 = REDcap.RCS02;
visits_tbl             = visits.RCS02;

% only takes REDcap surveys from AFTER Stage I implant
redcap   = redcap(ge(redcap.time, datetime(cfg.stage_dates(1),'TimeZone','America/Los_Angeles')),...
                  :); 

i_DevSet   = cellfun(@(x) ~isempty(x), db_RCSXXX.stimSettingsOut);

% takes stimSettingsOut to table
dev_w_rcap = sortrows(...
                    vertcat(db_RCSXXX.stimSettingsOut{i_DevSet}), 'HostUnixTime');

time_API = datetime(dev_w_rcap.HostUnixTime /1000,...
                'ConvertFrom','posixTime',...
                'TimeZone','America/Los_Angeles',...
                 'Format','dd-MMM-yyyy HH:mm:ss.SSS');

dev_w_rcap = addvars(dev_w_rcap, time_API ,'after',"HostUnixTime");

[~, i_unique] = unique(dev_w_rcap.time_API);

dev_w_rcap    = dev_w_rcap(i_unique,:);
%% organize INS logs and then infer contacts, amp, pw, and rate based off of DeviceSettings.txt

group_changes   = group_changes(~contains(group_changes.event, "Lead"),:);
group_changes   = renamevars(group_changes, 'time', 'time_INS');

rep_off         = strcmp(group_changes.event(1:end -1), 'Off') &...
                       strcmp(group_changes.event(2:end), 'Off');

group_changes(rep_off,:) = [];



i_on_off  = find(contains(group_changes.event, {'Off', 'On'}));
i_group   = find(contains(group_changes.event, 'Group'));



ther_stat = strings(1, length(i_group));

% based on nearest previous On/Off explicitly have Group X On or Off
for i = 1 : length(i_group)

    t_diff     = i_group(i) - i_on_off;

    near_t     = min(t_diff(t_diff > 0));
    j          = find(t_diff == near_t);

    ther_stat(i) ...
        = group_changes.event{i_on_off(j)};
             
end

group_changes = group_changes(i_group, :);
group_changes.therapyStatusDescription = ther_stat';


rep_off         = strcmp(group_changes.therapyStatusDescription(1:end -1), 'Off') &...
                       strcmp(group_changes.therapyStatusDescription(2:end), 'Off');

group_changes(rep_off,:) = [];

group_changes   = renamevars(group_changes, 'event', 'activeGroup');

group_changes.therapyStatusDescription   = cellstr(group_changes.therapyStatusDescription);
group_changes.activeGroup                = cellfun(@(x) {x(end)}, group_changes.activeGroup);

%% from previous StimLog.json get explicit stim params
for i = 1 : height(group_changes)

    t_diff   =  group_changes.time_INS(i) - dev_w_rcap.time_API;

    % find the smallest negative time difference to get nearest
    % stimLog.json settings BEFORE EventLog.txt settings

    near_t  = max(t_diff(t_diff < 0));
    j       = find(t_diff == near_t);
    
    % Since INS log can change groups SINCE last stimLog setting, find
    % explicit paramemters regardless of last activeGroup in StimLog

    act_group                     = group_changes.activeGroup{i};

    group                         = dev_w_rcap.(['Group', act_group])(j(1));

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


i_off = find(strcmp(group_changes.therapyStatusDescription, 'Off'));
off_rep = [diff(i_off); 0] == 1;

group_changes(i_off(off_rep), :) = [];




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


%%
tmp_group_changes         = group_changes;

groups                     = tmp_group_changes.activeGroup;
therapyStatusDescription   = tmp_group_changes.therapyStatusDescription;
time_INS                   = tmp_group_changes.time_INS;

group_tbl                  = table(time_INS, groups, therapyStatusDescription);

Didx                       = strcmp(groups, 'D') & strcmp(therapyStatusDescription, 'On');


group_tbl.groupD_ON_diff      = [0; diff(Didx)];

offline_cl_sess            = tmp_group_changes(Didx,:);





%% w/ explicit INS logs, now align to nearest REDcap

for i = 1 : height(redcap)

    t_diff   =  redcap.time(i) - group_changes.time_INS;
  
    near_t  = min(t_diff(t_diff > 0));


    % if no previous INSLog, use StimLog.json
    if isempty(near_t)
    % first INS Log started AFTER first month of REDcap -> use
    % DeviceSettings.json

        t_diff   = redcap.time(i) - dev_w_rcap.time_API;
        near_t   = min(t_diff(t_diff > 0));
        j        = find(t_diff == near_t);

        act_group       = dev_w_rcap.activeGroup{j};

        group           = dev_w_rcap.(['Group', act_group])(j);

        h               = find(group.validPrograms);
    
        redcap.therapyStatusDescription(i) = dev_w_rcap.therapyStatusDescription(j);
   
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


%%



%beh_stim_R.R_stimContacts(i_con) = [beh_stim_R.R_stimContacts{i_con}];


%%
%{

if amplitude is 0, then only return program 1 (0) contacts for simplity
here

%}




group_changes                 = movevars(group_changes,'time_INS', 'After','MPQcruel');
%end