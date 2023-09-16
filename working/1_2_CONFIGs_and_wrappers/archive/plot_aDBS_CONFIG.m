

%% below loads pain fluctuation study and PT meta data

%%% set-up working directories
cd(fullfile(dirs.rcs_analysis, 'working/'));         addpath(genpath(cd));

addpath(genpath(fullfile(dirs.rcs_analysis, 'code/')));

%%%  load REDcap from pain fluctuation study, and stages 1, 2, and 3
save_dir = fullfile(dirs.rcs_preproc , 'REDcap/');

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

%% stage dates and, home/clinic visits for RCS pts 1-7 w/ brief descriptions
%%% make_visit_dates was periodically updated throughout RCS trial

%%% load pt meta data (dates of stage starts/stops, and home/clinic visits for RCS pts)
CONFIG_pt_meta;
