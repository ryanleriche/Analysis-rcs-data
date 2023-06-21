% see CONFIG_ephy_analysis.m script to specify directories, pull REDcap, 
% and generate subdirectories
CONFIG_ephy_analysis;

%% import RCS databases, and INS logs per pt side
cfg_rcs                    = [];
cfg_rcs.load_EventLog      = true;

% option to load previous database for efficient processing
cfg_rcs.ignoreold_db                = false;
cfg_rcs.ignoreold_INS_logs          = false;
cfg_rcs.ignoreold_par_db            = false;
cfg_rcs.ignoreold_td_parsing        = false;

cfg_rcs.raw_dir                     = [dirs.rcs_pia, 'raw/'];
cfg_rcs.proc_dir                    = [dirs.rcs_pia, 'processed/'];
cfg_rcs.anal_dir                    = [dirs.rcs_pia, 'ephy_analysis/'];

cfg_rcs.pp_RCS_TD_subset            = 'stage1_only';
cfg_rcs.pp_fft                      = '30s_pre_survey';

%%
% ***name of the benchtop INS--not implanted w/n pt***
pt_side = 'RCS00L';

%%% process RCS .jsons into searchable database 
[db.(pt_side), bs.(pt_side)] ...
    ...
    = makeDatabaseRCS_Ryan(...
    ...
    cfg_rcs , pt_side);

%%% unpack all sense, LD, and stimulation settings as own variable in table
    % allows for programmatic discernment of unique RC+S settings
[par_db.(pt_side), ~]...
    ...
    = makeParsedDatabaseRCS(...
    ...
    cfg_rcs , pt_side, db);

%%
sess_dir = db.RCS00L.path{2};
cfg.textoutputs = true;

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

% combine data streams based off of packet loss and sampling rates
dataStreams      = {timeDomainData, PowerData, AdaptiveData, FFTData};
dataStreams      = dataStreams(~cellfun(@isempty, dataStreams));
rcs_streamed     = createCombinedTable(dataStreams, unifiedDerivedTimes);









