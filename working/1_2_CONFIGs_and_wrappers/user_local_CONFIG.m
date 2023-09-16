%% input directories unique to user w/n 'dirs' structure
%%% further directories will be programatically generated
dirs = struct;

% where RCS files are saved from PIA server
dirs.rcs_pia        = '/Users/Leriche/pia_server/datastore_spirit/human/rcs_chronic_pain/rcs_device_data/';

% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
% this a fork off of 'Analysis-rcs-data" as off April 2023
dirs.rcs_analysis   = '/Users/Leriche/Github/Analysis-rcs-data/';

%%% where DropBox desktop is saved locally
dirs.dropbox     = ['/Users/Leriche/Dropbox (UCSF Department of Neurological Surgery)/',...
                   'SUBNETS Dropbox/Chronic Pain - Activa and Summit 2.0'];

%%% If working w/ RC+S or NK neural spectra
% where processed RCS streaming sessions are saved
dirs.rcs_preproc_ss = '/Volumes/DBS Pain 3/rcs_device_data/processed/';

% where processed NK data are saved locally
dirs.nk_preproc     = '/Volumes/DBS Pain 3/nk_device_data/processed/';

%%% application programming interface (API) tokens
% essentially a password to access REDcap remotely, and is unique per 
% researcher per study (e.g., Ryan has a unique token for the RCS and PCS studies)

rcs_API_token   = '95FDE91411C10BF91FD77328169F7E1B';
pcs_API_token   = 'DB65F8CB50CFED9CA5A250EFD30F10DB';

%-----------↓↓↓ set-up working directories, no input needed ↓↓↓------------%
% done seperately to NOT add the FieldTrip toolbox recursively
cd(fullfile(dirs.rcs_analysis, 'working/'));         addpath(genpath(cd));

addpath(genpath(fullfile(dirs.rcs_analysis, 'code/')));