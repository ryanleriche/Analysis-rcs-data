%% see server_CONFIG.m script to specify directories and general set-up
server_CONFIG;

% load REDcap pain surveys and set-up directories
general_setup;

% see "dirs" structure for main directories used, and "rcs_db_cfg" for the
% configuration of the RC+S databasing
close all
set(0,'DefaultFigureVisible','off')
%% import RCS databases, and INS logs per pt side
% option to load previous database for efficient processing
cfg_rcs.ignoreold_db                = false;
cfg_rcs.ignoreold_INS_logs          = false;
cfg_rcs.ignoreold_par_db            = true;
cfg_rcs.ignoreold_td_parsing        = true;

cfg_rcs.pp_RCS_TD_subset            = 'stage1_only';
cfg_rcs.pp_fft                      = '5m_pre_survey';

% specify patient hemispheres
%%% pts to update database from scratch locally:
pt_sides        = {'RCS06R'};
%%
%%% main loop (per pt hemisphere)
for i = 1 : length(pt_sides)
    %%% process RCS .jsons into searchable database 

    [db.(pt_sides{i}), bs.(pt_sides{i})] ...
        ...
        = makeDatabaseRCS_Ryan(...
        ...
        cfg_rcs, pt_sides{i});
    
    %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings

    [par_db.(pt_sides{i}), ~]...
        ...
        = makeParsedDatabaseRCS(...
        ...
        cfg_rcs, pt_sides{i}, db);

    %%% find stim settings during each REDcap survey
    [stimLog.(pt_sides{i}), REDcap.(pt_sides{i})]  ...
        ...
        = align_REDcap_to_stimLog(...
        ...
        cfg_rcs, pt_sides{i}, db, REDcap);

    %%% save time-domain LFPs (from RC+S streaming sessions) 
        % across pt hemispheres as clearly labelled .mat
        % first entry allows subsetting of streaming sessions based on
        % criteria within fxn
    par_db_out.(pt_sides{i})  ...
        ...
        = pp_RCS_ss_TD( ...
        ...
        cfg_rcs.pp_RCS_TD_subset, pt_sides{i}, dirs.rcs_preproc,...
        par_db, stimLog);

  %%% save RCS time-domain data wrt REDcap survey
       ft_form_TD.(pt_sides{i}) ...
           = rcs_TD_peri_survey_2(cfg_rcs, dirs.rcs_preproc, pt_sides{i}, par_db_out, REDcap);
  
end
