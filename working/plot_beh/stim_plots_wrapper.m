% USE THIS TO PLOT PAIN METRICS by stim group
%  Prasad Shirvalkar Oct 25, 2022

% 
%% Set config variables and load databases, RCS logs, redcap etc 
clc
clear
cfg.pt_id                   = 'RCS02';

% Run the config script based on patient 
RCSconfig 

%% Load data 
%  This will load the Json databases, REDcap metrics and textlogs
RCSloaddata

%% RUN THIS IF you want to make the databases/ text logs again 
RCScompiledata 

% need to fix RCS logs to get all data including adaptive 
%  and load old log file and add to it

%% concatenate StimLog.json outputs from RCS database and then MAKE the stim groups
fprintf('Stim log concatenation:\n')
matchStimGroupToRedcap

%% PLOT pain metrics by stim group
   
plot_stim_groups