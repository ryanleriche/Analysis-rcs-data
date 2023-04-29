function [app_ss_tbl_out, INS_logs_out] ...
    ...
    = align_INSLogs_to_API_time(...
    ...
    pt_side_id, INS_logs, par_db, ss_var_oi)

    
INS_logs_RCSXXX       = INS_logs.(pt_side_id);
par_db_RCSXXX         = par_db.(pt_side_id);

exp_sense_state_vars  = ss_var_oi.(pt_side_id);

% td_fft_pm_ld_state_stim_tbl = TD_FFT_PB_LD_State_tbl.RCS02R; 
%%
%%%%%%%%%%%
%{
 with TD, FFT, power-bands, LD, and aDBS States, aligned and parsed from streamming sessions

   --> flesh-out aDBS offline session parameters, from previous streaming session
%}
%%%%%%%%%%%%
proc_app...
    = find_nearest_ss_to_INS_logs(par_db_RCSXXX, INS_logs_RCSXXX.app);

proc_g_changes...
    = find_nearest_ss_to_INS_logs(par_db_RCSXXX,  INS_logs_RCSXXX.group_changes);

%% explitcly save streaming sessions and offline aDBS sessions according to unique
% sensing parameters
var_oi  = {'chanFullStr', ...
         'bandFormationConfig', 'interval', 'size','streamSizeBins','streamOffsetBins', 'windowLoad',...
         'powerBandinHz', 'powerBinInHz','LD0', 'LD1', 'GroupDProg0_contacts',...
         'rise', 'fall','state'};

u_aDBS_Settings = exp_sense_state_vars(...
                      contains(exp_sense_state_vars, var_oi));

%%% for sake of finding unique parameters, put zeros for NaNs
%--> (NaNs treated as different values in unique fuction)

for i_var = 1 : length(u_aDBS_Settings)

    tmp_var = par_db_RCSXXX{:, u_aDBS_Settings{i_var}};

    if isnumeric(tmp_var)
        tmp_var(isnan(tmp_var)) = 0;
        par_db_RCSXXX{:, u_aDBS_Settings{i_var}} = tmp_var;

    end
end

%%% group sessions with same settings
par_db_RCSXXX.sess_w_same_settings ...
    = findgroups(par_db_RCSXXX(:,u_aDBS_Settings));

par_db_RCSXXX = movevars(par_db_RCSXXX, ...
    'sess_w_same_settings', 'After', 'sess_name');


%%% add sessions with same aDBS settings to app log itself
proc_app.sess_w_same_settings =nan(height(proc_app),1);

for j = 1:height(par_db_RCSXXX)

    i_same_sett = find(contains(proc_app.prev_sess_name, ...
                  par_db_RCSXXX.sess_name(...
                  par_db_RCSXXX.sess_w_same_settings == j)...
                         ));

    proc_app.sess_w_same_settings(i_same_sett) = j;

end

%proc_app = proc_app(~isnan(proc_app.sess_w_same_settings),:);

% only save streaming sessions (and their expanded parameters) of aDBS
% offline sessions
i_entry      = find(~strcmp(proc_app.prev_sess_name, 'No Entry'));
i_sess_oi    = find(contains(par_db_RCSXXX.sess_name, ...
                          unique(proc_app.prev_sess_name(i_entry))));

%%
% save field per pt
INS_logs_out.group_changes      = proc_g_changes;
INS_logs_out.app                = proc_app;

app_ss_tbl_out                  = par_db_RCSXXX(i_sess_oi, :);

end