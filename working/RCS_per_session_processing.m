function par_db_out = RCS_per_session_processing(cfg, pt_side_id, par_db, stimLog)

%%

dir_pt_side  = [cfg.proc_dir,'/', pt_side_id,'/', cfg.proc_subdir];

if ~exist(dir_pt_side, 'dir');      mkdir(dir_pt_side);    end


db_RCSXXX       = par_db.(pt_side_id);
stimLog_RCSXXX  = stimLog.(pt_side_id);


switch cfg.dates

    case 'DateRange'

        i_dates_oi = find(...
                            ge(db_RCSXXX.timeStart, cfg.date_range{1}) &...
                            le(db_RCSXXX.timeStart, cfg.date_range{2}));
        par_db_out  = db_RCSXXX(i_dates_oi, :);

    otherwise
        par_db_out   = db_RCSXXX;
end

i_suff_dur  = find(ge(par_db_out.duration,cfg.dur_range{1}) & ...
              le(par_db_out.duration, cfg.dur_range{2}));

par_db_out  = par_db_out(i_suff_dur, :);

%%

for i_sess = 1 : height(par_db_out)

    % initalize folder per session--> allows flexibility in N outputs
    dir_sess = [dir_pt_side, par_db_out.sess_name{i_sess}, '/'];

    % only pre-process in file does not exist OR explicitly re-running
    if ~isfile([dir_sess, 'rcs_streamed.mat']) || cfg.ignoreold

        sess_dir = par_db_out.path{i_sess};

        try
            % pull in data using Analysis-rcs-data
            [unifiedDerivedTimes,...
                timeDomainData, ~, ~,...
                ~, ~, ~, ...
                PowerData, ~, ~,...
                FFTData, ~, ~,...
                ...
                AdaptiveData, ~, ~, ~, ~,...
                ~, ~, ~, ~, ~,...
                ~, ~, ~, ...
                ~, ~] ...
                ...
                = ProcessRCS(cfg, sess_dir, 2);
            
            dataStreams      = {timeDomainData, PowerData, AdaptiveData, FFTData};
            dataStreams      = dataStreams(~cellfun(@isempty, dataStreams));
            rcs_streamed     = createCombinedTable(dataStreams, unifiedDerivedTimes);
            
       
            if ~isempty(rcs_streamed)
                [comb_dt_chunks, per_TD_lost]    = chunks_and_gaps(rcs_streamed);
                
                rcs_streamed   = renamevars(rcs_streamed, 'localTime', 'TimeInPST');
                
                par_db_out.per_TD_lost(i_sess)    = per_TD_lost;
                par_db_out.comb_dt_chunks(i_sess) = comb_dt_chunks;

                if ~isfolder(dir_sess);      mkdir(dir_sess);    end

                save([dir_sess, 'rcs_streamed.mat'], 'rcs_streamed', '-v7.3');

                par_db_out.loaded_status(i_sess)    = {'saved rcs_streamed'};


                %%% output stim_log.xls per streaming session
                i_stimLog     = find(strcmp(stimLog_RCSXXX.sess_name, par_db_out.sess_name{i_sess}));
                
                if ~isempty(i_stimLog)
                    stim_log  = removevars(stimLog_RCSXXX(i_stimLog, :), ...
                        {'GroupA', 'GroupB','GroupC', 'GroupD', ...
                        'redcap_btwn_stimLogs', 'i_redcap', 'redcap_reports'});
                    stim_log = renamevars(stim_log, 'time_stimLog', 'stimLogTimeInPST');
                    
                    writetable(stim_log,[dir_sess, 'session_stim_log.xls']); 
                end


            else
                par_db_out.proc_status(i_sess)    = {'reportedly_empty_rcs_file'};
            end

        catch
            par_db_out.proc_status(i_sess)    = {'failed_to_load'};

        end
        
    else
        par_db_out.proc_status(i_sess)    = {'saved_rcs_streamed'};
    end
end

save(...
     sprintf('%s/%s_par_db_oi.mat', dir_pt_side, pt_side_id), ...
     ...
     'par_db_out', '-v7.3');

writetable(par_db_out,  ...
    sprintf('%s/%s_par_db_oi.xlsx', dir_pt_side, pt_side_id))
end
