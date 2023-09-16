% the 'raw', 'proc' and 'ephy_analysis' directories are explicitly set
% depending on use case
rcs_db_cfg                             = [];
rcs_db_cfg.raw_dir                     = fullfile(dirs.rcs_pia, 'raw/');
rcs_db_cfg.proc_dir                    = fullfile(dirs.rcs_pia, 'processed/');
rcs_db_cfg.anal_dir                    = fullfile(dirs.rcs_pia, 'ephy_analysis/');

rcs_db_cfg.load_EventLog               = true;


%%% stage dates and, home/clinic visits for RCS pts 1-7 w/ brief descriptions
CONFIG_pt_meta;

pts               = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};
%% below loads pain surveys pre-trial and peri-trial
%%% (i.e., REDcap surveys from the pain fluctuation study (PFS), and stages 1, 2, and 3
save_dir = fullfile(rcs_db_cfg.proc_dir, 'REDcap/');

if isfile([save_dir, 'REDcap_PFS_RCS_pts.mat'])  
     tmp = load([save_dir, 'REDcap_PFS_RCS_pts.mat']);

     PFS           = tmp.PFS;
     PFS_sum_stats = tmp.PFS_sum_stats;

else  %%% per RCS pt, organize pain fluctuation study (pain reproting prioir to stage 0 [pre-sEEG trial period]) 
   
    PFS             = RCS_redcap_painscores(rcs_API_token, pcs_API_token, {'FLUCT'});
    
    % pain fluctuation study descriptive stats provides basis for clinical
    % outcomes of interest
    cfg               = [];
    cfg.dates         = 'AllTime';
    for i = 1:length(pts)
        PFS_sum_stats.(pts{i})  = calc_sum_stats(cfg, PFS.(pts{i}));
    end
    
    if ~isfolder(save_dir);    mkdir(save_dir);   end
    
    save([save_dir, 'REDcap_PFS_RCS_pts.mat'], 'PFS', 'PFS_sum_stats', '-v7.3');

end

%% import REDcap daily, weekly, and monthly surveys from stages 1,2 and 3
%  as of Apr. 2023, only daily surveys are analysis-ready/organized
REDcap   = RCS_redcap_painscores(rcs_API_token, [], [], pts, save_dir);

clear rcs_API_token pcs_API_token tmp save_dir