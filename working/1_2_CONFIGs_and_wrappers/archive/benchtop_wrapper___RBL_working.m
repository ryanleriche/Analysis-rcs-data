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



%%
N_hemi    = 2;      % 4/5 pts have bilateral implants
N_group   = 1;      % use ONE stim contact per shank-hemisphere


t_motif   = 2;      % 2s = (1.9s at 0 mA + 0.1 s at 2 mA) 
N_pulse   = 26;     % 25 pulses + 1 pulse padded at end to have 0mA before and after every pulse
N_u_amp   = 3;      % N unique amplitudes tested
N_iter    = 1;      % N iterations

t_trans   = 30;     % seconds of transition between groups


t_half_Hz = (t_motif * N_pulse + t_trans) * N_u_amp *N_iter;

%
t_motif  = 25;     % 25 s = (12.5 s at 0 mA + 12.5 s at 2 mA)
N_pulse  = 1;
N_u_amp  = 3;
N_iter   = 1;


t_2_Hz     = (t_motif * N_pulse + t_trans) * N_u_amp *N_iter;

t_ERP_sess = ((t_2_Hz + t_half_Hz) * N_hemi * N_group) / 60




