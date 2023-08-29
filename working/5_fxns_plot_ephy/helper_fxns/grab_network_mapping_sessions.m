function varargout = grab_network_mapping_sessions(cfg, pt_side_id, stimLog, par_db)

stim_log                 = stimLog.(pt_side_id);
[sess_names, i_u_sess]   = unique(stim_log.sess_name);

stim_log_ss_cell         = cell(length(sess_names), 1);


for i_sess    =  1 : length(sess_names)
    i_stimLog                   = find(strcmp(stim_log.sess_name, sess_names{i_sess}));
    stim_log_ss_cell{i_sess}    = stim_log(i_stimLog, :);
end

fprintf('N=%g sessions w/ stimLogs\n', length(stim_log_ss_cell));


stim_log_tbl                = table;
stim_log_tbl.time_stimLog   = stim_log.time_stimLog(i_u_sess);
stim_log_tbl.sess_name      = sess_names;
stim_log_tbl.logs           = stim_log_ss_cell;
stim_log_tbl.N_entries      = cellfun(@height, stim_log_ss_cell);



stim_log_tbl.N_2Hz_entries = cellfun(@(x) ...
                                        sum(...
                                        x.ampInMilliamps >0 &...
                                        x.rateInHz == 2 &...
                                        strcmp(x.therapyStatusDescription, 'On')...
                                        ), stim_log_ss_cell);

%%

save_dir = [cfg.proc_dir, '/', pt_side_id, '/',cfg.proc_subdir];

if ~isfolder(save_dir);    mkdir(save_dir);     end

save(...
     sprintf('%s/%s_stim_log_tbl.mat', save_dir, pt_side_id), ...
     ...
     'stim_log_tbl', '-v7.3');


network_mapping_tbl  = stim_log_tbl(stim_log_tbl.N_2Hz_entries >1, :);

%%%
par_db            = par_db.(pt_side_id);


i_network_mapping = ismember(par_db.sess_name, network_mapping_tbl.sess_name);

fprintf('N=%g network mapping sessions\n', sum(i_network_mapping));

par_db_oi         = par_db(i_network_mapping, :);


writetable(par_db_oi,  ...
    sprintf('%s/%s_par_db_out.xlsx', save_dir ,pt_side_id))


varargout{1} = par_db_oi;

end