% RCS LOAD DATA 
%  and loads necessary datasets.
% 
% Datasets to be loaded:
        % REDCAP pain metrics
        % RC+S json databases 
        % RCS textlogs 


%% Load data

%  Load REDCap painscores 
load(fullfile(processed_dir,'redcap_painscores.mat')); 
fprintf('\nREDcap painscores loaded from file \n')





 
% Load Neurophys RCS data
fprintf('Databases:\n')
database_dir = fullfile(processed_dir,'databases'); 
try
   temp1 = load(fullfile(database_dir,[cfg.pt_id 'L_database.mat']));
   db.([cfg.pt_id 'L']) = temp1.RCSdatabase_out;
   clear temp1
   fprintf(' - OK loaded:  %s \n',[cfg.pt_id 'L'])
catch 
    fprintf(' - NO database for %s \n',[cfg.pt_id 'L'])
end

try 
    temp2 = load(fullfile(database_dir,[cfg.pt_id 'R_database.mat']));
    db.([cfg.pt_id 'R']) = temp2.RCSdatabase_out;
    clear temp2
    fprintf(' - OK loaded:  %s \n',[cfg.pt_id 'R'])
catch 
    fprintf(' - NO database for %s  \n',[cfg.pt_id 'R'])
end







% Load the RCS textlogs (state changes from INS)
load(fullfile(processed_dir,'textlogs.mat'));
textlog = log; 
clear log
fprintf('Textlogs loaded OK\n')