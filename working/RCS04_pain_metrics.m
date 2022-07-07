%% user-inputs
% the patient's ID and specific side as outputted by RCS
PATIENTIDside         = 'RCS04Lt';

% where the RCS files are outputted (as saved on a Portable SSD)
rootdir               = '/Volumes/RBL_X5/';

% where the 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
github_dir            = '/Users/Leriche/Github/';


% application programming interface (API) token which is essentially a
% password to access REDcap remotely, and is unique per researcher per
% study (e.g., Ryan has a unique token for the RCS and PCS studies)

API_token             = '95FDE91411C10BF91FD77328169F7E1B';


%% pulls/organizes arms from REDcap (go into fxn to add new arms)
cd(github_dir);         addpath(genpath(github_dir));

pt_pain                 = RCS_redcap_painscores(API_token);

%% (RBL 06/28/22--in progress)
% import RCS files
patientroot_dir = fullfile(rootdir,char(regexp(PATIENTIDside,...
                        '\w*\d\d','match'))); %match the PATIENTID up to 2 digits: ie RCS02

scbs_dir        = fullfile(patientroot_dir,...
                    '/SummitData/SummitContinuousBilateralStreaming/', PATIENTIDside);



[RCSdatabase_out, badsessions] = ...
    makeDataBaseRCSdata(scbs_dir,'RCS04','ignoreold');




%% RCS04 plot daily metrics

%{
Specify, 'AllTime' for a 5 nearest report sliding average of all daily
reported pain metrics. Alternatively, to see the last n days of pain metrics
specify 'PreviousDays' followed by n days (without the sliding average).
%}

%{ 
Add overlay of stim, days on patch, contacts, etc., as optional arguements

(see introduction comment)
% 07/02/22 RBL: consider adding 'cfg' input w/ stim specifications, and pt
% ID as title

%}

plot_timeline(pt_pain.RCS04, 'AllTime');

% visually inspect to see distributions of the various pain metrics

plot_hist(pt_pain.RCS04,'AllTime');


% pain metric A "versus" pain metric B --> see how metrics visually covary

plot_versus(pt_pain.RCS04, 'AllTime');


