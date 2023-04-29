function [app_ss_tbl_out, INS_logs_out] ...
    ...
    = align_INSLogs_to_API_time(...
    ...
    pt_side_id, INS_logs, par_db, ss_var_oi)

% i                  = 3;
% pt_sides           = {'RCS02R','RCS05R','RCS05L'};
% pt_side_id         = pt_sides{i};

INS_logs_RCSXXX       = INS_logs.(pt_side_id);
par_db_RCSXXX         = par_db.(pt_side_id);

exp_sense_state_vars  = ss_var_oi.(pt_side_id);

%%  shift INS time based on 'TimeSync.json'--which returns the API-INS latency
%%% per packet--to assume similar latency when NOT streaming

% first for processed hemisphere's AppLog.txt
proc_app...
    ...
    = find_nearest_ss_to_INS_logs(...
    ...
par_db_RCSXXX, INS_logs_RCSXXX.app);

% then for processed hemisphere's EventLog.txt
proc_g_changes...
    ...
    = find_nearest_ss_to_INS_logs(...
    ...
par_db_RCSXXX,  INS_logs_RCSXXX.group_changes);

% given that Group and Status changes occur discretly, return merged Group
% and Status based on previous entry (i.e., the current device state)
proc_g_changes ...
    ...
    = determine_current_group_status(...
    ...
proc_g_changes);
%%  insert EventLog.txt entires btwn AppLog.txt entries
%%% AppLog.txt does NOT reflect Groups or therapy status changes

%{

find EventLog.txt entires occuring between AppLog.txt entries

report time, newstate, oldstate, program ampliudes *just like* proc_app
    AT time of modification
%}

%proc_g_changes(strcmp(proc_g_changes.Group_Status, 'GroupD_On'), :) = [];
%%
to_add_proc_app = cell(height(proc_app)+1,1);

for i_app = 1 : height(proc_app)

    if i_app < height(proc_app)

        i_grp_chan  =   find(...
                            ge(proc_g_changes.time_INS, proc_app.time_INS(i_app)) &...
                            le(proc_g_changes.time_INS, proc_app.time_INS(i_app+1))...
                            );
    else % return all following group changes since AppLog.txt entry

        i_grp_chan  =   find(...
                            ge(proc_g_changes.time_INS, proc_app.time_INS(i_app)) ...
                            );


    end

    if ~isempty(i_grp_chan)  % group changes the occured BETWEEN AppLog.txt entries
        
        % initalize group changes occuring between AppLog.txt entries
        % then flesh out row-by-row
        inter_grp_tbl = table;
        
        for j = 1 :  length(i_grp_chan)
        
            tmp_grp         = proc_g_changes(i_grp_chan(j),:);
            prev_sess_name  = tmp_grp.prev_sess_name;

            tmp_ss          = par_db_RCSXXX(...
                                    strcmp(par_db_RCSXXX.sess_name, prev_sess_name),...
                                    :);

            if isempty(inter_grp_tbl)
                oldstate = proc_app.newstate(i_app);
            else                      
                oldstate = inter_grp_tbl.newstate(j-1);
            end

            if oldstate == 15                    % therapy Off
               prog_amps   = [0 0 0 0 ];
               rateHz      = 0;

            elseif any(oldstate == 20:22) % Group A, B, or C, respectively, therapy On

               activeGroup = inter_grp_tbl.event_id{j-1}(6);
               amp_vars    = compose(...
                                  [sprintf('Group%s', activeGroup), 'Prog%g_ampInMilliamps'], 0:3);

               prog_amps   = tmp_ss{:, amp_vars};
               rateHz      = tmp_ss{:,sprintf('Group%s_rateInHz', activeGroup)};

            elseif oldstate == -1
                %state_str    = sprintf('state%g', oldstate);
                var_oi       = compose('GroupDProg%g_ampInMilliamps',0:3);

                prog_amps    = tmp_ss{:, var_oi};

                rateHz       = tmp_ss{:,'GroupD_rateInHz'};

            else
                state_str    = sprintf('state%g', oldstate);
                var_oi       = compose('_AmpInMilliampsProg%g',0:3);
                
                amp_state_oi = cellfun(@(x) ...
                                        [state_str, x], var_oi,...
                                        'UniformOutput', false);
                prog_amps    = tmp_ss{:, amp_state_oi};

                rateHz       = tmp_ss{:,'GroupD_rateInHz'};

            end
            prog_amps(prog_amps == -1) = 8.5;

            prog0mA = prog_amps(1);        prog1mA = prog_amps(2);
            prog2mA = prog_amps(3);        prog3mA = prog_amps(4);

            if tmp_grp.therapyStatus

                  switch tmp_grp.Group_Status{1}
                    % hard coding as open-loop status
                    case 'GroupA_On';            newstate = 20;  status = 3;  
                    case 'GroupB_On';            newstate = 21;  status = 3;  
                    case 'GroupC_On';            newstate = 22;  status = 3;

                    % unique case where GroupD is On, but BEFORE any state
                    % changes defaults to RLP settings (not necessarily in
                    % line w/ any Adaptive State (pg. 85)

                    case 'GroupD_On';            newstate = -1;  status = 2;  proc_app.oldstate(i_app+1) = -1;

                  end  
            else  % therapy is OFF
                  % see Medtronic.NeuroStim.Olympus.DataTypes.Therapy.Adaptive 
                  % for 'InActive' status
                newstate = 15;   status = 0; 
            end

            time_INS   = tmp_grp.time_INS;
            event_id   = tmp_grp.Group_Status;

            tmp_app    = table(time_INS, status, newstate, oldstate, ...
                            prog0mA, prog1mA, prog2mA, prog3mA, rateHz, ...
                            event_id, prev_sess_name);

            inter_grp_tbl(j,:)       = tmp_app;
        end
       
        % save all inter-group changes within cell-array
        to_add_proc_app{i_app} = inter_grp_tbl;

    end
end



%%% merged EventLog.txt + AppLog.txt output
ol_cl_changes      = sortrows(...
                            [vertcat(to_add_proc_app{:}, proc_app)],...
                            'time_INS'...
                            );
%% explitcly save streaming sessions and offline aDBS sessions according to unique
% sensing parameters
var_oi  = {'chanFullStr', ...
         'bandFormationConfig', 'interval', 'size',...
         'streamSizeBins','streamOffsetBins', 'windowLoad',...
         'powerBandinHz', 'powerBinInHz',...
         ...
         'LD0', 'LD1', 'GroupDProg0_contacts',...
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


%%% add sessions with same aDBS settings to INSLog itself
ol_cl_changes.sess_w_same_settings =nan(height(ol_cl_changes ),1);

ol_cl_changes.prev_sess_name{...
    cellfun(@isempty, ol_cl_changes.prev_sess_name)} = '';


for i_par = 1:height(par_db_RCSXXX)

    i_same_sett = find(contains(ol_cl_changes.prev_sess_name, ...
                  par_db_RCSXXX.sess_name(...
                  par_db_RCSXXX.sess_w_same_settings == i_par)...
                         ));

    ol_cl_changes.sess_w_same_settings(i_same_sett) = i_par;

end

%proc_app = proc_app(~isnan(proc_app.sess_w_same_settings),:);

% only save streaming sessions (and their expanded parameters) of aDBS
% offline sessions
i_entry      = find(~strcmp(ol_cl_changes.prev_sess_name, 'No Entry'));
i_sess_oi    = find(contains(par_db_RCSXXX.sess_name, ...
                          unique(ol_cl_changes.prev_sess_name(i_entry))));

% 
% doi = datetime({'09-Mar-2023', '17-Mar-2023'}, 'TimeZone','America/Los_Angeles');
% 
% i_par_db = ge(par_db_RCSXXX.timeStart, doi(1)) & le(par_db_RCSXXX.timeStart, doi(2));
% 
% 
% parsed_database = par_db_RCSXXX(i_par_db,...
%      {'sess_name','timeStart', 'timeStop', 'GroupDProg0_ampInMilliamps'});

ambig_state   = find(ol_cl_changes.newstate(1:end-1) ~= ol_cl_changes.oldstate(2:end));
per_ambig     = 100*length(ambig_state) / height(ol_cl_changes);

fprintf('%s | %.3f%% of state changes where newstate(i-1) ~= old state(i)\n', pt_side_id, per_ambig)

%%
% save field per pt
INS_logs_out.group_changes      = proc_g_changes;
INS_logs_out.app                = proc_app;

% Apr 2023 RBL attempt to stitch group changes in app changes
INS_logs_out.ol_cl_changes      = ol_cl_changes;

app_ss_tbl_out                  = par_db_RCSXXX(i_sess_oi, :);



end