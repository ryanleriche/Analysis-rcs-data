function [dbs_oi, ss_tbl_oi] =...
    ...
    pull_INSLog_SS_oi(...
    ...
    cfg, pt_hemi_id, INS_logs, app_SS_tbl)



dbs_log       = INS_logs.(pt_hemi_id).ol_cl_changes;
app_ss_tbl    = app_SS_tbl.(pt_hemi_id);


% option to run specific dates
if strcmp(cfg.dates, 'DateRange') == 1

    date_range    = datetime(cfg.date_range, 'TimeZone', 'America/Los_Angeles', 'InputFormat','dd-MMM-uuuu');

    i_entries     = ge(dbs_log.time_INS, date_range(1)- duration('72:00:00')) & ...
                    le(dbs_log.time_INS, date_range(2));
    dbs_oi        = dbs_log (i_entries, :);

    % output percent of old states with MATCHING new states definitions
    ambig_state   = find(dbs_oi.newstate(1:end-1) ~= dbs_oi.oldstate(2:end));
    per_ambig     = 100*length(ambig_state) / height(dbs_oi);

    fprintf(['%s | between %s -> %s | ',...
             '%.3f%% of state changes where newstate(i-1) ~= old state(i)\n'],...
             pt_hemi_id, date_range, per_ambig)


    i_entries   = ge(app_ss_tbl.timeStart, date_range(1) - duration('72:00:00')) & ...
                  le(app_ss_tbl.timeStart, date_range(2));
    ss_tbl_oi   = app_ss_tbl(i_entries,:);

% use all entries if 'AllTime' is specified
elseif strcmp(cfg.dates, 'AllTime') == 1      
    dbs_oi     = dbs_log;
    ss_tbl_oi  = app_ss_tbl;
end

end