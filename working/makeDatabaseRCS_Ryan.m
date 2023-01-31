function [RCSdatabase_out, varargout] = makeDatabaseRCS_Ryan(cfg)
%{
function database_out = makeDatabaseRCS_Ryan(cfg.raw_dir)


This function creates a database of rcs data

INPUT:   cfg.raw_dir should be the root folder of the Session Files (e.g. SCBS folder)
              e.g. cfg.raw_dir = '/Desktop/[PATIENTID]/'
              e.g. cfg.raw_dir = '/Volumes/Prasad_X5/RCS02/;


          PATIENTIDside = should specify what the Patient ID name is, and
          should be same as subfolder (e.g. 'RCS02R')

          (OPTIONAL INPUT) 'ignoreold' will ignore old databases and
          start fresh. Third input should be lowercase as written.



OUTPUT:   RCSdatabase_out is a table with fields ordered in importance, which will
           be saved as a mat file and csv to cfg.raw_dir

          (OPTIONAL)
          SECOND OUTPUT could be provided to collect list of bad sessions
          (e.g. those with no data in Jsons)

             Included fields are:
                'rec',[],...
                'timeStart',[],...    % <-- 07/20/22 RBL: replacing time with             
                'timeStop',[],...     % <-- timeStart and timeStop 
                'sessname',[],...
                'duration',[],...
                'battery',[],...
                'TDfs',[],...
                'fft',[],...
                'power',[],...
                'stim',[],...
                'stimName',[],...
                'stimparams',[],...
                'matExist',[],...
                'path',[],...
                'powerbands',[],...
                'aDBS',[],...
                'adaptive_threshold',[],...
                'adaptive_onset_dur',[],...
                'adaptive_termination_dur',[],...
                'adaptive_states',[],...
                'adaptive_weights',[],...
                'adaptive_updaterate',[],...
                'adaptive_pwrinputchan',[]);

******* NOTE THAT ONLY LD0 data is populated in the adaptive fields *******


 USING CELL DATA IN THE TABLE:
      to concatenate all cell variables in the table (such as duration)
      use:
          alldurations =cat(1, database_out.duration{:})


 **** This will check to see if there is an existing database, if so, it will
 update that database table.
 *****


Depedencies:
https://github.com/JimHokanson/turtle_json
in the a folder called "toolboxes" in the same directory as the processing scripts


Prasad Shirvalkar Sep 13,2021
%}

%% 
%{
Note throughout fxn a "RCS database", "database", etc., is technically a 
MATLAB table containing description of a given recording session WITHOUT 
time series data (voltage, FFT, power).

___________________________________________________________________________

Optimize logic for checking for old databse

Add to DataBase:
    X session starts and stops
    X parsed stim parameters

    X duty cycle
        X figure out unit conversion into seconds 
        - see Metronic.SummitAPI
    X find impedance wrt to each stimparams output
        - reported in EventLog.json each Lead Integrity check
    X see below, but figure out sessions w/ multiple stim parameters and
    effective time stamping
        - return as table w/ each change per session
    X pain report timestamps
        - REDcap timestamps are darn close to API logs (see
        'db_sort_beh()')

Ryan Leriche Aug 20, 2022

%}


tic

% match the PATIENTID up to 2 digits: ie RCS02
pt_rootdir  = fullfile(cfg.raw_dir,char(regexp(cfg.pt_id_side,...
                        '\w*\d\d','match'))); 


%  Define the directories to search in (SCBS and aDBS)
scbs_dir_root       = fullfile(pt_rootdir,'/SummitData/SummitContinuousBilateralStreaming/', cfg.pt_id_side);
adbs_dir_root       = fullfile(pt_rootdir, '/SummitData/StarrLab/', cfg.pt_id_side);

scbs_sess_dirs      = findFilesBVQX(scbs_dir_root,'Sess*',struct('dirs',1,'depth',1));
adbs_sess_dirs      = findFilesBVQX(adbs_dir_root,'Sess*',struct('dirs',1,'depth',1));
sess_dirs           = [scbs_sess_dirs; adbs_sess_dirs];


db_out      = struct('rec',[],...
                'timeStart',[],...
                'timeStop',[],...
                'sess_name','',...
                'duration',[],...
                'timeDomainSettings',[],...
                'fftSettings', [],...
                'powerSettings',[],...
                'metaData',[],...
                'stimSettingsOut',[],...
                'stimMetaData', [],...
                'stimLogSettings', [],...
                'eventLogTable', [],...
                'path', [],...
                'stimReg', []...
                );

%%
% insert section here to load old database, and just add rows to it if
% needed, so as not to replicate whole thing.
% Can be turned off with third input 'ignoreold'

proc_cfg.raw_dir      = [cfg.raw_dir(1:end-4), 'processed/databases/'];

outputFileName    = fullfile(proc_cfg.raw_dir,[cfg.pt_id_side '_database.mat']);


if isfile(outputFileName) && ~cfg.ignoreold

    disp('Loading previously saved database');
    D = load(outputFileName,'RCSdatabase_out','badsessions');

    old_database         = D.RCSdatabase_out;
    old_badsessions      = D.badsessions;
    old_sess             = D.RCSdatabase_out.sess_name;
    old_badsess          = D.badsessions.sess_name;
    
    olddirs = contains(sess_dirs, old_sess) | contains(sess_dirs, old_badsess) ;
    sess_dirs(olddirs)= [];

    if isempty(sess_dirs)

        fprintf("No new data to add!  Existing database returned \n")
        RCSdatabase_out = old_database;
        varargout{1}= old_badsessions;

        return
    end

elseif isfile(outputFileName) && cfg.ignoreold

    disp('Ignoring previous database(s) and compiling from scratch...')
    old_database= [];

else

    disp('Compiling database from scratch...')
    old_database= [];

end

%d = find(contains(sess_dirs, 'Session1673830359519'))

%%
for d = 1: length(sess_dirs)

    % run through all of sessions (aDBS then SCBS sessions)
    given_sess = findFilesBVQX(sess_dirs{d},'Device*',struct('dirs',1,'depth',1));


    fprintf('Reading folder %d of %d  \n', d, length(sess_dirs))

    json_files = dir([given_sess{:}, '/*.json']);

    if isempty(json_files) % no data exists inside

        db_out(d).timeStart = [];
        [~,fn] = fileparts(sess_dirs{d});
        db_out(d).sess_name = fn;
        disp('no data.. moving on');

    else % data may exist, check for time domain data

        %clear devicepath settingsfile 
        td_file      = findFilesBVQX(sess_dirs{d},'EventLog.json');
        dev_file     = findFilesBVQX(sess_dirs{d},'DeviceSettings.json');

        if ~isempty(td_file) && ~isempty(dev_file)  % time data file doesn't exist

            %
            [~,fn]                      = fileparts(sess_dirs{d});
            db_out(d).sess_name         = fn;
            [path,~,~]                  = fileparts(td_file{1});
            db_out(d).path              = path;


            % extract times and .mat status
            % load device settings file
            settingsfile                 = findFilesBVQX(sess_dirs{d},'DeviceSettings.json');
            [devicepath,~,~]            = fileparts(settingsfile{1});
            try

                if isfile(settingsfile)
    
                    [TD_settings_raw, timeDomainSettings, powerSettings, fftSettings, metaData] = ...
                        createDeviceSettingsTable(devicepath);
    
                    % save all outputs from 'DeviceSettings.json'

                    db_out(d).TD_settings_raw       = TD_settings_raw;
                    db_out(d).timeDomainSettings    = timeDomainSettings;
                    db_out(d).fftSettings           = fftSettings;
                    db_out(d).powerSettings         = powerSettings;
                    db_out(d).metaData              = metaData;
    
                % IF there is time/ power domain data
                    if ~isempty(timeDomainSettings) && ~isempty(powerSettings)
        
                        % Get recording start time/ duration
                        timeStart = timeDomainSettings.timeStart / 1000;
                        timeStop = timeDomainSettings.timeStop / 1000;
                        

                    % edge case where SCBS was logged into, but streaming
                    % session never starts
                    % --> pull TD settings anyways (01/24/23; RBL)
                    elseif ~isempty(TD_settings_raw)

                        timeStart = TD_settings_raw(1,:).time /1000;
                        timeStop  = TD_settings_raw(end,:).time /1000;


                    end
                    % accounts for daylight savings and timezones based
                    % off of UTC offset
                    timeFormat = sprintf('%+03.0f:00', metaData.UTCoffset);
                    
                    timeStart = datetime(timeStart,...
                    'ConvertFrom','posixTime','TimeZone',timeFormat,...
                    'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                    
                    timeStop = datetime(timeStop,...
                    'ConvertFrom','posixTime','TimeZone',timeFormat,...
                    'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                    
                    
                    
                    db_out(d).timeStart = timeStart;
                    db_out(d).timeStop  = timeStop;
                    db_out(d).duration  = duration(timeStop - timeStart,'Format','hh:mm:ss.SSS');


    
                else
                warning('No DeviceSettings.json file')
                end
            
            catch
            end
        end

        
        try

            % Get stim settings
            [stimSettingsOut, stimMetaData] ...
                ...
                = ...
                ...
            createStimSettingsFromDeviceSettings(devicepath);

            db_out(d).stimSettingsOut = stimSettingsOut;
            db_out(d).stimMetaData    = stimMetaData;


            % load StimLog.json (stim parameters, active group and therapy status)
            stimfile          = findFilesBVQX(sess_dirs{d},'StimLog.json');
            [stimpath,~,~]    = fileparts(stimfile{1});

        catch
        end

        try
            [stimLogSettings] = createStimSettingsTable(stimpath, stimMetaData);

             if ~isempty(timeDomainSettings) && ~isempty(powerSettings)

                 % Get recording start time/ duration
                time_stimLog = stimLogSettings.HostUnixTime / 1000;

                % accounts for daylight savings and timezones based
                % off of UTC offset
                timeFormat = sprintf('%+03.0f:00', metaData.UTCoffset);
                
                stimLogSettings.time_stimLog = datetime(time_stimLog,...
                'ConvertFrom','posixTime','TimeZone',timeFormat,...
                'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        
             end

            db_out(d).stimLogSettings = stimLogSettings;

            db_out(d).stimLogSettings.stimParams_prog1 = ...
                stimLogSettings.stimParams_prog1;
                
            stimnamegroup        = {'A','B','C','D'; '1' , '5', '9','13'};

            [~,j]                = find(contains(stimnamegroup,stimLogSettings.activeGroup));
            stimname             =  metaData.stimProgramNames(str2double(stimnamegroup{2,j(1)}));

            db_out(d).stimReg    =  stimname{1};    

          
    
        catch
        end

        if cfg.load_EventLog  == true
            try
                prt_folder                 = findFilesBVQX(sess_dirs{d},'EventLog.json');
                [eventlog_path,~,~]        = fileparts(prt_folder{1});
    
                if isfile(prt_folder)
                    [db_out(d).eventLogTable] = createEventLogTable(eventlog_path);
                else
                    warning('No EventLog.json file')
                end
    
            catch
            end         
        end

        try

            % Get Adaptive settings info
            [db_out(d).DetectorSettings, ...
                db_out(d).AdaptiveStimSettings, ...
                db_out(d).AdaptiveEmbeddedRuns_StimSettings] ...
                ...
                = ...
                ...
             createAdaptiveSettingsfromDeviceSettings(devicepath);

        catch
        end
        
    end
end


database_out                = struct2table(db_out,'AsArray',true);

if ~cfg.ignoreold
    database_out.rec        = (1:height(database_out))' + height(old_database);
else
    database_out.rec        = (1:height(database_out))';
end

% delete all rows with empty session names (WHY DOES THIS OCCUR?)
database_out            = database_out(cellfun(@(x) ~isempty(x),database_out.sess_name),:);

RCSdatabase_out         = sortrows(database_out, 'sess_name'); %sorting by session name

% when database does not have any missing fields it is returned as a struct
% array rather than a cell array of structures

if isstruct(RCSdatabase_out.metaData)

% -> meaning we need to convert to a cell array of structures before we
% concatenate with the previous database

    for i = 1: height(RCSdatabase_out)

        RCSdatabase_out.metaData_new(i) = {RCSdatabase_out.metaData(i)};
    end

    RCSdatabase_out = removevars(RCSdatabase_out , 'metaData');

    RCSdatabase_out = renamevars(RCSdatabase_out ,'metaData_new', 'metaData');
end

%% clear empty session rows and assign to new variable 'badsessions'
if iscell(RCSdatabase_out.timeStart)
    loc = cellfun('isempty', RCSdatabase_out{:,'timeStart'});
else
    loc = isempty(RCSdatabase_out.timeStart);
end

badsessions = RCSdatabase_out(loc,:);

%% COMBINE WITH OLD DATABASE
% IF the old database existed, recombine with new database and sort it
% but first fix cell/ mat class issues

if ~isempty(old_database)
    disp('combining with old database...');

    if exist('RCSdatabase_out', 'var')

        if ~isempty(RCSdatabase_out)

            varnames = old_database.Properties.VariableNames;
            
            % if the new sessions don't have variable names as old sessions
            % initalize them as blank--can happen if no new sessions were aDBS
            for i = 1 : length(varnames)
                if ~any(strcmp(RCSdatabase_out.Properties.VariableNames, varnames{i}))

            
                   RCSdatabase_out.(varnames{i}) = cell(height(RCSdatabase_out), 1);
                end
            end

            RCSdatabase_out = [old_database; RCSdatabase_out];
        else
            RCSdatabase_out  = old_database;

            disp(['Database as of ', datestr(old_database.timeStart{end}(1)), '; ',...
            old_database.sess_name{end},'---no new files since ALL were empty'])
        end
    
        if ~isempty(badsessions)


           varnames = old_badsessions.Properties.VariableNames;
            
            % if the new sessions don't have variable names as old sessions
            % initalize them as blank--can happen if no new sessions were aDBS
            for i = 1 : length(varnames)
                if ~any(strcmp(badsessions.Properties.VariableNames, varnames{i}))

                   badsessions.(varnames{i}) = cell(height(badsessions), 1);
                end
            end


            badsessions = [old_badsessions; badsessions];
        else
            badsessions = old_badsessions;

        end
        
    else
        RCSdatabase_out  = old_database;
        varargout{1}     = old_badsessions;

        disp(['Old database as of ', datestr(old_database.timeStart{end}), '; ',...
            old_database.sess_name{end}])

    end
end

% OUTPUTS =========================================================
if nargout == 2
    varargout{1} = badsessions;
end

% Rename file to include patient ID
% writetable(RCSdatabase_out,fullfile(proc_cfg.raw_dir,[cfg.pt_id_side...
%     '_database.csv']));

save(fullfile(proc_cfg.raw_dir,[cfg.pt_id_side '_database.mat']),...
    'RCSdatabase_out','badsessions');

fprintf('mat of database saved as %s to %s \n',...
    [cfg.pt_id_side '_database.mat'],proc_cfg.raw_dir);

%==================================================================

end