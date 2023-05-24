%% load CONFIG_server_wrapper

[dirs,    rcs_API_token,   pcs_API_token, ... -> input and output directories, and API tokens
 PFS,     PFS_sum_stats,...                   -> pain fluctuation study (PFS) data and summary statistics
 pt_META, stage_dates]...                     -> hard-coded pt meta data (RCS Stage dates pulled from patient iternaries, Google Drive folders, etc
...
    = CONFIG_server_wrapper;

%% import REDcap daily, weekly, and monthly surveys from stages 1,2 and 3
% as of Apr. 2023, only daily surveys are analysis-ready/organized
REDcap          = RCS_redcap_painscores(rcs_API_token);

%% import RCS databases, and INS logs per pt side
cfg                    = [];
cfg.load_EventLog      = true;

% option to load previous database for efficient processing
cfg.ignoreold_db                = false;
cfg.ignoreold_INS_logs          = false;
cfg.ignoreold_par_db            = true;

cfg.ignoreold_td_parsing        = false;

cfg.raw_dir                     = [dirs.rcs_pia, 'raw/'];
cfg.proc_dir                    = dirs.rcs_preproc;
cfg.anal_dir                    = [dirs.rcs_pia, 'ephy_analysis/'];

% specify patient hemispheres
%%% pts to update database from scratch locally:
%pt_sides        = {'RCS02R', 'RCS05R', 'RCS05L', 'RCS06R','RCS06L','RCS07L', 'RCS07R'};

pt_sides        = {'RCS02R'};

cfg.pp_RCS_TD_subset   = 'stage1_only';
cfg.pp_fft             = '30s_pre_survey';

%%% main loop (per pt hemisphere)
for i = 1 : length(pt_sides)
    %%% process RCS .jsons into searchable database 

    [db.(pt_sides{i}), bs.(pt_sides{i})] ...
        ...
        = makeDatabaseRCS_Ryan(...
        ...
        cfg, pt_sides{i});
    
    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings

    [par_db.(pt_sides{i}), ~]...
        ...
        = makeParsedDatabaseRCS(...
        ...
        cfg, pt_sides{i}, db);

    %%% find stim settings during each REDcap survey
    [stimLog.(pt_sides{i}), REDcap.(pt_sides{i})]  ...
        ...
        = align_REDcap_to_stimLog(...
        ...
        cfg, pt_sides{i}, db, REDcap);

    %%% save time-domain LFPs (from RC+S streaming sessions) 
        % across pt hemispheres as clearly labelled .mat
        % first entry allows subsetting of streaming sessions based on
        % criteria within fxn
    par_db_out.(pt_sides{i})  ...
        ...
        = pp_RCS_ss_TD( ...
        ...
        cfg.pp_RCS_TD_subset, pt_sides{i}, dirs.rcs_preproc,...
        par_db, stimLog);

  %%% save RCS time-domain data wrt REDcap survey
       ft_form_TD.(pt_sides{i}) ...
           = rcs_TD_peri_survey(cfg, dirs.rcs_preproc, pt_sides{i}, par_db_out, REDcap);

  %%% save FFT per channel across sessions
    % summary spectrograms over sessions (i.e., days to months worth of
    % spectra saved in cfg.proc_dir subfolders)
        fft_taper_comparison(cfg, dirs, pt_sides{i}, ft_form_TD)
        

end



