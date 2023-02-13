function [app_ss_tbl, INS_proc, INS_ss_g_changes] ...
    ...
    = align_stream_sess_to_INSLogs(...
    ...
    INS_logs, par_db_RCSXXX, exp_sense_state_vars)



proc_app                    = INS_logs.app;

% td_fft_pm_ld_state_stim_tbl = TD_FFT_PB_LD_State_tbl.RCS02R; 
%%
%%%%%%%%%%%
%{
 with TD, FFT, power-bands, LD, and aDBS States, aligned and parsed from streamming sessions

   --> flesh-out aDBS offline session parameters, from previous streaming session
%}
%%%%%%%%%%%%
%% from nearest previous streaming session to INS log entry

% more explicit naming
proc_app     = renamevars(proc_app,  'time', 'time_INS');
proc_app     = proc_app(...
                    ge(proc_app.time_INS, ...
                    datetime('2001-Mar-01', 'TimeZone', 'America/Los_Angeles')), ...
                :);

tmp_ss_tbl   = par_db_RCSXXX(~isnan(par_db_RCSXXX.INS_lat_mean), :);

%eq(tmp_ss_tbl.INS_lat_mean, duration('0:00:0)'))
% (...
%                    strcmp(td_fft_pm_ld_state_stim_tbl.activeGroup, 'D') &...
%                    ge(td_fft_pm_ld_state_stim_tbl.duration, duration('0:05:00')), :);

%%
for i = 1 : height(proc_app)

    t_diff   =  proc_app.time_INS(i) - tmp_ss_tbl.timeStart + tmp_ss_tbl.INS_lat_mean;

    % find the smallest negative time difference to get nearest
    % stimLog.json settings BEFORE AppLog.txt settings

    near_t  = min(t_diff(t_diff > 0));

    % need to have streaming session occuring BEFORE INS log time
    if ~isempty(near_t)
        j       = find(t_diff == near_t);
        j       = j(1);

        proc_app.prev_sess_name(i)    = tmp_ss_tbl.sess_name(j);

    else
         proc_app.prev_sess_name(i) = {'No Entry'};
    end
end

[~,i_in_ss_tbl]   = ismember(proc_app.prev_sess_name, tmp_ss_tbl.sess_name);

shift_by_time     = tmp_ss_tbl.INS_lat_mean(i_in_ss_tbl);

proc_app.time_INS = proc_app.time_INS - shift_by_time;


%% align group changes to API time
%
%
% organize INS logs and then infer contacts, amp, pw, and rate based off of DeviceSettings.txt
%group_changes   = group_changes(~contains(group_changes.event, "Lead"),:);
proc_g_changes   = renamevars(INS_logs.group_changes, 'time', 'time_INS');
%
%

% remove instances of identical repeated events
i_rep_events      = diff([inf; findgroups(proc_g_changes.event)])~=0;

proc_g_changes   = proc_g_changes(i_rep_events,:);



i_on_off  = find(contains(proc_g_changes.event, {'Off', 'On'}));
i_group   = find(contains(proc_g_changes.event, 'Group'));

% i_lead_int   = find(contains(group_changes.event, 'Lead'));
% unique(group_changes.event(i_lead_int+1))
tmp_groups       = proc_g_changes.event(i_group, :);
tmp_ther_stat    = cell(length(i_group),1);

% based on nearest previous On/Off explicitly have Group X On or Off
for i = 1 : length(i_group)

    t_diff     = i_group(i) - i_on_off;

    near_t     = min(t_diff(t_diff > 0));
    j          = find(t_diff == near_t);

    tmp_ther_stat{i} ...
        = proc_g_changes.event{i_on_off(j)};
             
end

% w/ previous therapyStatus assigned to Group, remove all therapyStatuses,
% and any repeats
tmp_concat = cellfun(@(x,y) [x,'_',y], tmp_groups, tmp_ther_stat, 'UniformOutput', false);

proc_g_changes.event(i_group) = tmp_concat;

proc_g_changes(i_on_off,:) = [];

i_rep_events    = diff([inf; findgroups(proc_g_changes.event)])~=0;
proc_g_changes   = proc_g_changes(i_rep_events,:);

%% w/ parsimonious group changes --> align INS time to API time
proc_g_changes      = proc_g_changes(...
                    ge(proc_g_changes.time_INS, ...
                    datetime('2001-Mar-01', 'TimeZone', 'America/Los_Angeles')), ...
                        :);


for i = 1 : height(proc_g_changes)

    t_diff   =  proc_g_changes.time_INS(i) - tmp_ss_tbl.timeStart + tmp_ss_tbl.INS_lat_mean;

    % find the smallest negative time difference to get nearest
    % stimLog.json settings BEFORE AppLog.txt settings

    near_t  = min(t_diff(t_diff > 0));

    % need to have streaming session occuring BEFORE INS log time
    if ~isempty(near_t)
        j       = find(t_diff == near_t);
        j       = j(1);

        proc_g_changes.prev_sess_name(i)    = tmp_ss_tbl.sess_name(j);

    else
         proc_g_changes.prev_sess_name(i) = {'No Entry'};
    end
end

proc_g_changes    = proc_g_changes(...
                        ~strcmp(proc_g_changes.prev_sess_name, 'No Entry'), :);

[~,i_in_ss_tbl]   = ismember(proc_g_changes.prev_sess_name, tmp_ss_tbl.sess_name);

shift_by_time     = tmp_ss_tbl.INS_lat_mean(i_in_ss_tbl);

proc_g_changes.time_INS = proc_g_changes.time_INS - shift_by_time;
INS_proc.group_changes  = proc_g_changes;

%% merge/validate Group changes from streaming sessions and INS Log

ss_group_on_off_tbl = par_db_RCSXXX(:, {'sess_name','timeStop', 'activeGroup', 'therapyStatusDescription'});


ss_group_on_off_tbl.event = cellfun(@(x,y) ['Group',x,'_',y],...
                             ss_group_on_off_tbl.activeGroup, ...
                             ss_group_on_off_tbl.therapyStatusDescription, ...
                                'UniformOutput', false);

timeStop_vec                = repmat({'timeStop_ss'}, height(ss_group_on_off_tbl),1);

ss_group_on_off_tbl          = renamevars(ss_group_on_off_tbl, 'timeStop', 'time_align');
ss_group_on_off_tbl.t_origin = timeStop_vec ;


proc_g_changes               = renamevars(proc_g_changes, 'time_INS', 'time_align');
proc_g_changes.t_origin      = repmat({'time_INS'}, height(proc_g_changes),1);


INS_ss_g_changes         = sortrows([proc_g_changes(:, {'time_align', 'event', 't_origin'});...
                                    ss_group_on_off_tbl(:, {'time_align', 'event', 't_origin'})],...
                                    'time_align');

% combine ANY Group_Off into single Off 
INS_ss_g_changes.event(contains(INS_ss_g_changes.event, 'Off'))      = {'Off'};


i_rep_events           = diff([inf; findgroups(INS_ss_g_changes.event)])~=0;
INS_ss_g_changes   = INS_ss_g_changes(i_rep_events,:);
%% explitcly save streaming sessions and offline aDBS sessions according to unique
% sensing parameters

var_oi  = {'chanFullStr', ...
         'bandFormationConfig', 'interval', 'size','streamSizeBins','streamOffsetBins', 'windowLoad',...
         'powerBandinHz', 'powerBinInHz','LD0', 'LD1', 'GroupDProg0_contacts',...
         'rise', 'fall','state'};


u_aDBS_Settings = exp_sense_state_vars(...
                      contains(exp_sense_state_vars, var_oi));

par_db_RCSXXX   = par_db_RCSXXX(ge(par_db_RCSXXX.duration, duration('0:05:00')), :);

sess_w_same_settings  = findgroups(par_db_RCSXXX(:,u_aDBS_Settings));


par_db_RCSXXX = addvars(par_db_RCSXXX, sess_w_same_settings,...
                          'After', 'sess_name');

proc_app.sess_w_same_settings =nan(height(proc_app),1);

for i = 1:length(unique(par_db_RCSXXX.sess_w_same_settings))

    i_same_sett = find(contains(proc_app.prev_sess_name, ...
                  par_db_RCSXXX.sess_name(...
                  par_db_RCSXXX.sess_w_same_settings == i)...
                         ));

    proc_app.sess_w_same_settings(i_same_sett) = i;


end

proc_app = proc_app(~isnan(proc_app.sess_w_same_settings),:);

INS_proc.app = proc_app;

% only save streaming sessions (and their expanded parameters) of aDBS
% offline sessions
i_entry      = cellfun(@(x) ~isempty(x), proc_app.prev_sess_name);

i_sess_oi    = find(contains(par_db_RCSXXX.sess_name, ...
                          unique(proc_app.prev_sess_name(i_entry))));

app_ss_tbl   = par_db_RCSXXX(i_sess_oi, :);



% app_SS_tbl.RCS02R                   = app_ss_tbl;
% proc_App_log.RCS02R                 = proc_app;

end