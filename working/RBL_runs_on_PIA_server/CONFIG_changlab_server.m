% input and output directories, and API tokens
% call pain fluctuation study (PFS) data and summary statistics
% call hard-coded pt meta data (RCS Stage dates pulled from patient iternaries, Google Drive folders, etc)

%% input directories unique to user w/n 'dirs' structure
dirs = struct;

% where RCS files are saved from PIA server
dirs.rcs_pia           = '/datastore_spirit/human/rcs_chronic_pain/rcs_device_data/';

% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
% this a fork off of 'Analysis-rcs-data" as off April 2023
dirs.rcs_analysis      = '/home/rleriche/Analysis-rcs-data/';

% where processed RCS streaming sessions are saved
dirs.rcs_preproc       = '/home/rleriche/rcs_device_data/processed/';

% application programming interface (API) token which is essentially a
% password to access REDcap remotely, and is unique per researcher per
% study (e.g., Ryan has a unique token for the RCS and PCS studies)

rcs_API_token   = '95FDE91411C10BF91FD77328169F7E1B';
pcs_API_token   = 'DB65F8CB50CFED9CA5A250EFD30F10DB';


% paths where raw, processed, and analysis data/figures are saved
cfg_rcs.raw_dir     = [dirs.rcs_pia, 'raw/'];
cfg_rcs.proc_dir    = dirs.rcs_preproc;
cfg_rcs.anal_dir    = '/home/rleriche/rcs_device_data/ephy_analysis/';

%%
% below loads pain fluctuation study and PT meta data
%-------------------------------------------------------------------------%
%------------------------↓↓↓ no input needed ↓↓↓--------------------------%
%-------------------------------------------------------------------------%

%%% set-up working directories
cd(fullfile(dirs.rcs_analysis, 'working/'));         addpath(genpath(cd));

%%% set-up original 'Analysis-rcs data' repo
addpath(genpath(fullfile(dirs.rcs_analysis, 'code/')))

% where top-level FieldTrip folder is saved (do not add recursively)
dirs.ft     = fullfile(dirs.rcs_analysis, 'fieldtrip-20220912/');
addpath(dirs.ft);

% initalize FieldTrip
ft_defaults;


%%%  load REDcap from pain fluctuation study, and stages 1, 2, and 3
save_dir = [dirs.rcs_preproc, '/REDcap/'];

if isfile([save_dir, 'REDcap_PFS_RCS_pts.mat'])  
     tmp = load(...
         [save_dir, 'REDcap_PFS_RCS_pts.mat']);

     PFS           = tmp.PFS;
     PFS_sum_stats = tmp.PFS_sum_stats;


else  %%% per RCS pt, organize pain fluctuation study (pain reproting prioir to stage 0 [pre-sEEG trial period]) 
   
    PFS             = RCS_redcap_painscores(rcs_API_token, pcs_API_token, {'FLUCT'});
    
    % pain fluctuation study descriptive stats provides basis for clinical
    % outcomes of interest
    cfg               = [];
    cfg.dates         = 'AllTime';
    pts               = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};
    for i = 1:length(pts)
        PFS_sum_stats.(pts{i})  = calc_sum_stats(cfg, PFS.(pts{i}));
    end
    
    if ~isfolder(save_dir);    mkdir(save_dir);   end
    
    save([save_dir, 'REDcap_PFS_RCS_pts.mat'], ...
         ...
         'PFS', 'PFS_sum_stats', '-v7.3');

end

clear tmp save_dir

%%% load pt meta data (dates of stage starts/stops, and home/clinic visits for RCS pts)
CONFIG_pt_meta;

close all
set(0,'DefaultFigureVisible','off')