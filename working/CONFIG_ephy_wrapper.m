function  [dirs, rcs_API_token, pcs_API_token, ...  input and output directories, and API tokens
           ...
           PFS, PFS_sum_stats,...                   pain fluctuation study (PFS) data and summary statistics
           ...
           pt_META, stage_dates]...                 hard-coded pt meta data (RCS Stage dates pulled from patient iternaries, Google Drive folders, etc
           ...
           = CONFIG_ephy_wrapper

%% input directories unique to user w/n 'dirs' structure
dirs = struct;


% where RCS files are saved from PIA server
dirs.rcs_pia    = '/Users/Leriche/pia_server/datastore_spirit/human/rcs_chronic_pain/rcs_device_data/';

% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
% this a fork off of 'Analysis-rcs-data" as off April 2023
dirs.rcs_analysis      = '/Users/Leriche/Github/';


% where processed RCS streaming sessions are saved
dirs.rcs_preproc        = '/Volumes/DBS Pain 3/rcs_device_data/processed/';

% where DropBox desktop is saved locally
dirs.dropbox     = ['/Users/Leriche/Dropbox (UCSF Department of Neurological Surgery)/',...
                   'SUBNETS Dropbox/Chronic Pain - Activa and Summit 2.0'];


% application programming interface (API) token which is essentially a
% password to access REDcap remotely, and is unique per researcher per
% study (e.g., Ryan has a unique token for the RCS and PCS studies)

rcs_API_token   = '95FDE91411C10BF91FD77328169F7E1B';
pcs_API_token   = 'DB65F8CB50CFED9CA5A250EFD30F10DB';

%% below loads pain fluctuation study and PT meta data
%%% set-up working directories
cd([dirs.rcs_analysis, 'Analysis-rcs-data/working']);         

addpath(genpath([dirs.rcs_analysis, 'Analysis-rcs-data/']));

%%%  load REDcap from pain fluctuation study, and stages 1, 2, and 3
save_dir = [dirs.rcs_pia,'processed/REDcap/'];

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

%% %stage dates and, home/clinic visits for RCS pts 1-7 w/ brief descriptions
%%% make_visit_dates was periodically updated throughout RCS trial

[pt_META, stage_dates]   = make_visit_dates;


end