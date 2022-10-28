% RCSconfig defines all general paths, variables, etc...
% 
%  RUN this general pain RC+S config file at the beginning of every script or wrapper for plotting
%  or behavioral analysis 
% 
%  
%  Define stage 0,1,2,3 start dates and visit dates for each subject 

%% CHANGE  THIS CELL FOR YOUR LOCAL/ REMOTE PATHS
% where RCS files are saved locally 

data_raw_dir = '/Volumes/PrasadX5/spiritdata/raw/';  
    processed_dir = [data_raw_dir(1:end-4) '/processed/'];

% where 'Analysis-rcs-data' Github repo is saved locally
github_dir            = '/Users/pshirvalkar/Documents/GitHub/RyanAnalysis-rcs-data/';

API_token             = '95FDE91411C10BF91FD77328169F7E1B';



% goto working path
cd([github_dir, 'working']);        
addpath(genpath(github_dir));  


%% Define stage start dates, home, and clinic visits for RCS pts 1-7
stage_dates.RCS02 = {'08-Sep-2020'; '31-Jan-2021'; '31-May-2022'}; % RCS02
stage_dates.RCS04 =  {'13-May-2021'; '12-Jul-2021'}; % RCS04
stage_dates.RCS05 =  {'21-Jul-2021';'07-Sep-2021'}; % RCS05
stage_dates.RCS06 =   {'07-Jun-2022','18-Aug-2022'}; % RCS06
stage_dates.RCS07 =   {'19-Oct-2022'}; % RCS07

cfg.stage_dates = stage_dates;

visits              = struct;
visits.RCS02        = table;
visits.RCS02.dates  = datetime(...
   {'08-Sep-2020'; '09-Sep-2020'; '10-Sep-2020'; '11-Sep-2020';...
    '13-Oct-2020';...
    '18-Oct-2020'; ...
    '02-Nov-2021'; '09-Nov-2021';...
    '01-Feb-2021'; '02-Feb-2021';  '03-Feb-2021';...
    '13-Apr-2021'; '14-Apr-2021';...
    '27-Sep-2021'; '28-Sep-2021';...
    '31-May-2022';...
    '28-Jun-2022'; '29-Jun-2022'}, ...
    ...
    'TimeZone', 'America/Los_Angeles');

visits.RCS02.desc   = ...
    {'s1_implant';     's1_inpatient';    's1_inclinic_d1';  's1_inclinic_d2';...
     'remove_hardware_inpatient';...
     's2_start';
     'washout_testing'; 'washout_testing';
     's2_inclinic_d1'; 's2_inclinic_d2' ; 's2_inclinic_d3';...
     's2_home_d1'    ; 's2_home_d2';...
     's2_inclinic_d1'; 's2_inclinic_d2' ;...
     's3_start';...
     's3_inclinic_d1'; 's3_inclinic_d2'};


visits.RCS04        = table;
visits.RCS04.dates  = datetime(...
   {'13-May-2021'; '14-May-2021'; '15-May-2021'; '16-May-2021';...   
    '12-Jul-2021'; '13-Jul-2021'; '14-Jul-2021';...
    '17-Aug-2021'; '18-Aug-2021'; '19-Aug-2021';...
    '13-Jul-2022'}, ...
    ...
    'TimeZone', 'America/Los_Angeles');

visits.RCS04.desc   = ...
    {'s1_implant';      's1_inpatient_d1';  's1_inpatient_d2'; 's1_inpatient_d3'; ...
     's2_inclinic_d1';  's2_inclinic_d2' ;  's2_inclinic_d3';...
     's2_home_d1'    ;  's2_home_d2';       's2_home_d2';...
     's2_home_d1'    ...
     };


visits.RCS05        = table;
visits.RCS05.dates  = datetime(...
   {'21-Jul-2021';   '22-Jul-2021';    '23-Jul-2021';...
    '07-Sep-2021';   '08-Sep-2021';    '09-Sep-2021';...
    '01-Dec-2021';   '02-Dec-2021';...
    '10-Aug-2022';   '11-Aug-2022';...
    '05-Sep-2022'},...
    ...
    'TimeZone', 'America/Los_Angeles');

visits.RCS05.desc   = ...
    {'s1_implant';      's1_inpatient_d1';   's1_inpatient_d2';...
    's2_inclinic_d1';   's2_inclinic_d2' ;   's2_inclinic_d3';...
    's2_home_d1'    ;   's2_home_d2';...
    's2_inclinic_d1';   's2_inclinic_d2';...
    's3_start'};


visits.RCS06        = table;
visits.RCS06.dates  = datetime(...
   {'01-Mar-2022';...   
    '07-Jun-2021'}, ...
    ...
    'TimeZone', 'America/Los_Angeles');

visits.RCS06.desc   = ...
    {'s1_implant';...
     's2'};


cfg.visits = visits;
clear stage_dates visits
fprintf('Stage_dates and visits added to cfg \n')
disp(cfg)

