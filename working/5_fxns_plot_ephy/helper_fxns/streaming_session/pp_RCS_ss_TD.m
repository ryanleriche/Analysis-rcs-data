% pre-process (pp) RC+S streaming session time domain data
function par_db_out  ...
        ...
        = pp_RCS_ss_TD( ...
        ...
        cfg_rcs, pt_side_id, proc_base_dir,...
        par_db, stimLog)

% pre-processing all RCS streaming sessions would lead to excessive
% storage (~50 GB per hemisphere per month)
% --> pp only needed subsets of data


ss_subset = cfg_rcs.pp_RCS_TD_subset;

switch ss_subset
    case 'network_mapping'

        cfg                  = [];
        cfg.proc_dir         = fullfile(proc_base_dir, 'streaming_sessions');
        cfg.proc_subdir      = 'network_mapping';
        
        cfg.dates            = 'AllTime';
        
        cfg.ignoreold        = false;
        cfg.textoutputs      = false;
        
        cfg.dur_range        = {duration('0:00:00'), duration('03:00:00')};
    
    
        nw_mapp_ss.(pt_side_id) = grab_network_mapping_sessions(cfg, pt_side_id, stimLog, par_db);
    
        par_db_out = ...
            RCS_per_session_processing(cfg, cfg_rcs, pt_side_id, nw_mapp_ss, stimLog);


    case 'stage1_only'
        %%% day of stage 1 implant -> day before stage 2 (prior to
        %%% analgesic stim testing)
        % Note that RCS06 and RCS07 underwent analgesic stim testing during
        % stage 1 visit (see parsed databases, and CRC notes in Google
        % Drive)
        [s1_start, s1_end] = pull_s1_dates(par_db.(pt_side_id));

        %%% cfg settings across pt hemispheres
        cfg                  = [];
        
        cfg.proc_dir         = fullfile(proc_base_dir, 'streaming_sessions');
        cfg.proc_subdir      = 'stage1_only';
        
        cfg.dates            = 'DateRange';
        
        cfg.ignoreold        = false;
        cfg.textoutputs      = false;

        cfg.date_range       =  {s1_start, s1_end - duration('24:00:00')};

        cfg.dur_range        = {duration('0:02:00'), duration('03:00:00')};
    
        par_db_out ...
            = RCS_per_session_processing(cfg, pt_side_id, par_db, stimLog);

end
end