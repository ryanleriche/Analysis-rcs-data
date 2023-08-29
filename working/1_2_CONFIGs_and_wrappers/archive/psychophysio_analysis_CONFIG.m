%% input directories unique to user w/n 'dirs' structure
dirs = struct;

% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
% this a fork off of 'Analysis-rcs-data" as off April 2023
dirs.rcs_analysis      = '/Users/Leriche/Github/Analysis-rcs-data/';


% where DropBox desktop is saved locally
dirs.dropbox     = ['/Users/Leriche/Dropbox (UCSF Department of Neurological Surgery)/',...
                   'SUBNETS Dropbox/Chronic Pain - Activa and Summit 2.0'];



%-------------------------------------------------------------------------%
%------------------------↓↓↓ no input needed ↓↓↓--------------------------%
%-------------------------------------------------------------------------%
%% loads pain fluctuation study and PT meta data
% application programming interface (API) token which is essentially a
% password to access PRISM_r_cap remotely, and is unique per researcher per
% study (e.g., Ryan has a unique token for the RCS and PCS studies)

rcs_API_token   = '95FDE91411C10BF91FD77328169F7E1B';
pts             = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};

%%% set-up working directories
cd(fullfile(dirs.rcs_analysis, 'working/'));         addpath(genpath(cd));


dirs.working       = fullfile(dirs.dropbox,'DATA ANALYSIS',...
                             '(Ryan) Stage1_2_3_group_level_behavioral_analysis',...
                             'REDcap_FitBit_only/');
%%% load pt meta data (dates of stage starts/stops, and home/clinic visits for RCS pts)
CONFIG_pt_meta;

%%% import REDcap of daily, weekly, and monthly surveys from stages 1,2 and 3
% as of Apr. 2023, only daily surveys are analysis-ready
PRISM_r_cap                  = RCS_redcap_painscores(rcs_API_token);

save_dir = fullfile(dirs.working, 'prism_data/REDcap/');

% save each pt's pain surveys as easily shared .xlsx
imp_field = fieldnames(PRISM_r_cap);

% save as source data as .csv
for i = 1 : length(imp_field)

    switch imp_field{i}
        case pts

            source_xlsx = [save_dir, imp_field{i}, '_stages123.xlsx'];
            
            if exist(source_xlsx, 'file') == 2;        delete(source_xlsx);   end
            
            writetable(PRISM_r_cap.(imp_field{i}), source_xlsx);

    end
end

clear save_dir imp_field cfg
