function par_db_oi = putative_network_mapping(pt_side_id, stimLog, par_db)

stim_log                 = stimLog.(pt_side_id);


[sess_names, i_u_sess]   = unique(stim_log.sess_name);

stim_log_ss_cell         = cell(length(sess_names), 1);


for i_sess    =  1 : length(sess_names)
    i_stimLog                   = find(strcmp(stim_log.sess_name, sess_names{i_sess}));
    stim_log_ss_cell{i_sess}    = stim_log(i_stimLog, :);
end


stim_log_tbl                = table;
stim_log_tbl.time_stimLog   = stim_log.time_stimLog(i_u_sess);
stim_log_tbl.sess_name      = sess_names;
stim_log_tbl.logs           = stim_log_ss_cell;
stim_log_tbl.N_entries      = cellfun(@height, stim_log_ss_cell);


stim_log_tbl.N_2Hz_entries = cellfun(@(x) ...
                                        sum(...
                                            x.ampInMilliamps >0 &...
                                             x.rateInHz == 2 &...
                                             strcmp(x.therapyStatusDescription, 'On') &...
                                             [diff(x.time_stimLog) > duration('0:00:05');1]), stim_log_ss_cell);


network_mapping_tbl  = stim_log_tbl(stim_log_tbl.N_2Hz_entries >1, :);

%%
par_db            = par_db.(pt_side_id);


i_network_mapping = ismember(par_db.sess_name, network_mapping_tbl.sess_name);
par_db_oi         = par_db(i_network_mapping, :);

end