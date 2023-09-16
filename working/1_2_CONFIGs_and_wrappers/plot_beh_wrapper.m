% see XXX_CONFIG script to specify directories, pull REDcap, 
% and generate subdirectories
plot_beh_CONFIG;

%% Behavioral plotting 

% Create a config file to define patients and dates of interest
cfg_rcap                     = [];

% Define patient of interest
cfg_rcap.pt_id               = 'RCS05';

% Define dates of interest
cfg_rcap.dates               = 'PreviousDays';
cfg_rcap.ndays               = 10;

% Define stim parameter
cfg_rcap.stim_parameter      = '';

% Define stage dates
cfg_rcap.stage_dates         = stage_dates;

% Use the following line to return every pain score ever reported:
%cfg.dates        = 'AllTime';

cfg_rcap.subplot             = true;
cfg_rcap.stim_parameter      = '';

% Generate a graph of behavioral reports (with pain fluctuation study aka
% PFS)
%     plot_timeline(cfg_rcap, REDcap, PFS_sum_stats);

% Generate a graph of behavioral reports (without PFS information)
    plot_timeline_no_PFS(cfg_rcap, REDcap);    

%% Visualize additional patients of interest

% Create a config file to define patients and dates of interest
cfg_rcap                     = [];

% Define patient of interest
cfg_rcap.pt_id               = 'RCS02';

% Define dates of interest
cfg_rcap.dates               = 'PreviousDays';
cfg_rcap.ndays               = 10;

% Define stim parameter
cfg_rcap.stim_parameter      = '';

% Define stage dates
cfg_rcap.stage_dates         = stage_dates;

% Use the following line to return every pain score ever reported:
%cfg.dates        = 'AllTime';

cfg_rcap.subplot             = true;
cfg_rcap.stim_parameter      = '';

% Generate a graph of behavioral reports (with pain fluctuation study aka
% PFS)
%     plot_timeline(cfg_rcap, REDcap, PFS_sum_stats);

% Generate a graph of behavioral reports (without PFS information)
    plot_timeline_no_PFS(cfg_rcap, REDcap);    

    % Create a config file to define patients and dates of interest
cfg_rcap                     = [];

% Define patient of interest
cfg_rcap.pt_id               = 'RCS02';

% Define dates of interest
cfg_rcap.dates               = 'PreviousDays';
cfg_rcap.ndays               = 10;

% Define stim parameter
cfg_rcap.stim_parameter      = '';

% Define stage dates
cfg_rcap.stage_dates         = stage_dates;

% Use the following line to return every pain score ever reported:
%cfg.dates        = 'AllTime';

cfg_rcap.subplot             = true;
cfg_rcap.stim_parameter      = '';

% Generate a graph of behavioral reports (with pain fluctuation study aka
% PFS)
%     plot_timeline(cfg_rcap, REDcap, PFS_sum_stats);

% Generate a graph of behavioral reports (without PFS information)
    plot_timeline_no_PFS(cfg_rcap, REDcap);    