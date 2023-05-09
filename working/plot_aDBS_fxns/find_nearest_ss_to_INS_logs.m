function ins_log = find_nearest_ss_to_INS_logs(par_db_RCSXXX, ins_log)

%% from nearest previous streaming session to INS log entry
% more explicit naming
ins_log    = renamevars(ins_log,  'time', 'time_INS');
ins_log    = ins_log(...
                    ge(ins_log.time_INS, ...
                    datetime('2001-Mar-01', 'TimeZone', 'America/Los_Angeles')), ...
                :);



%% for streaming session durations of sufficent duration
% --> shift INS log entries

tmp1_ss_tbl    = par_db_RCSXXX(~isnan(par_db_RCSXXX.INS_lat_mean), :);

i_ss           = ge(tmp1_ss_tbl.duration, duration('0:01:00'));
tmp1_ss_tbl    = tmp1_ss_tbl(i_ss, :); 

%%
for j = 1 : height(ins_log)

    t_diff   =  ins_log.time_INS(j) - tmp1_ss_tbl.timeStart + tmp1_ss_tbl.INS_lat_mean;

    % find the smallest negative time difference to get nearest
    % streaming session BEFORE AppLog.txt settings

    near_t  = min(t_diff(t_diff > 0));

    % need to have streaming session occuring BEFORE INS log time
    if ~isempty(near_t)
        h       = find(t_diff == near_t);
        h       = h(1);

        ins_log.prev_sess_name(j)    = tmp1_ss_tbl.sess_name(h);

    else
        ins_log.prev_sess_name(j) = {'No Entry'};
    end
end

[~,i_in_ss_tbl]   = ismember(ins_log.prev_sess_name, tmp1_ss_tbl.sess_name);

shift_by_time     = tmp1_ss_tbl.INS_lat_mean(i_in_ss_tbl(i_in_ss_tbl~=0));

ins_log.time_INS(i_in_ss_tbl~=0)...
    = ins_log.time_INS(i_in_ss_tbl~=0) - shift_by_time;

%% w/ shifted INS logs, now find nearest streaming session regardless of duration
tmp2_ss_tbl    = par_db_RCSXXX(~isnan(par_db_RCSXXX.INS_lat_mean), :);

for j = 1 : height(ins_log)

    t_diff   =  ins_log.time_INS(j) - tmp2_ss_tbl.timeStart;

    % find the smallest negative time difference to get nearest
    % stimLog.json settings BEFORE AppLog.txt settings

    near_t  = min(t_diff(t_diff > 0));

    % need to have streaming session occuring BEFORE INS log time
    if ~isempty(near_t)
        h       = find(t_diff == near_t);
        h      = h(1);

        ins_log.prev_sess_name(j)    = tmp2_ss_tbl.sess_name(h);

    else
         ins_log.prev_sess_name(j)   = {'No Entry'};
    end
end

