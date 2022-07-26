function [RCSdatabase_out, varargout] = makeDataBaseRCSdata(dirname, ptIDside, varargin)
%{
function database_out = makeDataBaseRCSdata(dirname)


This function creates a database of rcs data

INPUT:   DIRNAME should be the root folder of the Session Files (e.g. SCBS folder)
              e.g. DIRNAME = '/Desktop/[PATIENTID]/'
              e.g. DIRNAME = '/Volumes/Prasad_X5/RCS02/;


          PATIENTIDside = should specify what the Patient ID name is, and
          should be same as subfolder (e.g. 'RCS02R')

          (OPTIONAL INPUT) 'ignoreold' will ignore old databases and
          start fresh. Third input should be lowercase as written.



OUTPUT:   RCSdatabase_out is a table with fields ordered in importance, which will
           be saved as a mat file and csv to DIRNAME

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

    * duty cycle
        * (close, but) figure out unit conversion into seconds 
        - "12, 12" could be unix time

    * find impedance wrt to each stimparams output

    * see below, but figure out sessions w/ multiple stim parameters and
    effective time stamping

    
    * pain report timestamps
        - search within the 'log.json' for when the "Complete Pain Report"
        button w/n SCBS was pressed

        - search for further time stamps

        -> eventually, compare to REDcap timestamps (start, vas start, mpq
        start (?))

read through all information contained w/n the adaptive.json, log.json, and
mirror.json


Ryan Leriche Jul 20, 2022

%}


tic

% match the PATIENTID up to 2 digits: ie RCS02
pt_rootdir  = fullfile(dirname,char(regexp(ptIDside,...
                        '\w*\d\d','match'))); 


%  Define the directories to search in (SCBS and aDBS)
scbs_dir_root       = fullfile(pt_rootdir,'/SummitData/SummitContinuousBilateralStreaming/', ptIDside);
adbs_dir_root       = fullfile(pt_rootdir, '/SummitData/StarrLab/', ptIDside);

scbs_sess_dirs      = findFilesBVQX(scbs_dir_root,'Sess*',struct('dirs',1,'depth',1));
adbs_sess_dirs      = findFilesBVQX(adbs_dir_root,'Sess*',struct('dirs',1,'depth',1));
sess_dirs           = [scbs_sess_dirs; adbs_sess_dirs];


db_out      = struct('rec',[],...
                'timeStart',[],...
                'timeStop',[],...
                'sess_name','',...
                'duration',[],...
                'battery',[],...
                'TDfs',[],...
                'fft',[],...
                'power',[],...
                'stim',[],...
                'stimName','',...
                'stimparams','',...
                'matExist',[],...
                'path','',...
                'powerbands',[],...
                'aDBS',[],...
                'adaptive_threshold',[],...
                'adaptive_onset_dur',[],...
                'adaptive_termination_dur',[],...
                'adaptive_states',[],...
                'adaptive_weights',[],...
                'adaptive_pwrinputchan',[],...
                'adaptive_updaterate',[]);

%%
% insert section here to load old database, and just add rows to it if
% needed, so as not to replicate whole thing.
% Can be turned off with third input 'ignoreold'

proc_dirname      = [dirname(1:end-3), 'processed/databases/'];

outputFileName    = fullfile(proc_dirname,[ptIDside '_database.mat']);


if isfile(outputFileName) && nargin < 3

    disp('Loading previously saved database');
    D = load(outputFileName,'RCSdatabase_out','badsessions');

    old_database = D.RCSdatabase_out;
    old_badsessions = D.badsessions;
    oldsess = D.RCSdatabase_out.sessname;
    oldbadsess = D.badsessions.sessname;
    
    olddirs = contains(sess_dirs,oldsess) | contains(sess_dirs,oldbadsess) ;
    sess_dirs(olddirs)= [];

    if isempty(sess_dirs)

        fprintf("No new data to add!  Existing database returned \n")
        RCSdatabase_out = old_database;
        varargout{1}= old_badsessions;

        return
    end

elseif isfile(outputFileName) && strcmp(varargin{1},'ignoreold')

    disp('Ignoring previous database(s) and compiling from scratch...')
    old_database= [];

else

    disp('Compiling database from scratch...')
    old_database= [];


end





%%
for d = 1 : length(sess_dirs)

    % run through all of sessions (aDBS then SCBS sessions)
    given_sess = findFilesBVQX(sess_dirs{d},'Device*',struct('dirs',1,'depth',1));

     
    if nargin == 2 &&  d > numel(scbs_sess_dirs)
        db_out(d).aDBS = 1;
    else
        db_out(d).aDBS= 0;
    end

    fprintf('Reading folder %d of %d  \n', d, length(sess_dirs))

    if isempty(given_sess) % no data exists inside

        db_out(d).time = [];
        db_out(d).matExist  = 0;
        [~,fn] = fileparts(sess_dirs{d});
        db_out(d).sessname = fn;
        disp('no data.. moving on');

    else % data may exist, check for time domain data

        clear devicepath settingsfile 
        td_file      = findFilesBVQX(sess_dirs{d},'EventLog.json');
        dev_file     = findFilesBVQX(sess_dirs{d},'DeviceSettings.json');

        if ~isempty(td_file) && ~isempty(dev_file)  % time data file doesn't exist

            %
            [~,fn] = fileparts(sess_dirs{d});
            db_out(d).sessname = fn;
            [path,~,~] = fileparts(td_file{1});
            db_out(d).path = path;


            % extract times and .mat status
            % load device settings file
            try
                settingsfile        = findFilesBVQX(sess_dirs{d},'DeviceSettings.json');
                [devicepath,~,~]    = fileparts(settingsfile{1});

                [timeDomainSettings, powerSettings, ...
                    fftSettings, metaData] ...
                    ...
                    = ...
                    ...
                 createDeviceSettingsTable(devicepath);


                %IF there is time/ power domain data
                if ~isempty(timeDomainSettings) && ~isempty(powerSettings)


                    %   Get recording start time/ duration
                    sessTime = [timeDomainSettings.timeStart, ...
                                timeDomainSettings.timeStop] /1000;

                    timeFormat = sprintf('%+03.0f:00',metaData.UTCoffset);

                    sessDt = datetime(sessTime,...
                        'ConvertFrom','posixTime','TimeZone',timeFormat,...
                        'Format','dd-MMM-yyyy HH:mm:ss.SSS');


                    db_out(d).timeStart = sessDt(1);
                    db_out(d).timeStop  = sessDt(2);
                    db_out(d).duration  = duration(sessDt(2)-sessDt(1),'Format','hh:mm:ss.SSS');




                    % Get time domain sensing info
                    db_out(d).TDfs = timeDomainSettings.samplingRate;
                    db_out(d).TDchan0= timeDomainSettings.chan1{1};
                    db_out(d).TDchan1= timeDomainSettings.chan2{1};
                    db_out(d).TDchan2= timeDomainSettings.chan3{1};
                    db_out(d).TDchan3= timeDomainSettings.chan4{1};
                    db_out(d).battery = metaData.batteryLevelPercent;



                    %  Get FFT length info
                    if ~isnan(fftSettings.recNum)
                        db_out(d).fft = fftSettings.fftConfig.size;
                    end


                    % Get powerbands and whether recorded info

                    if ~isnan(powerSettings.recNum)
                        db_out(d).power = 1 ;
                        db_out(d).powerbands = powerSettings.powerBands.powerBandsInHz;
                    end



               end
            catch
            end



            try

                % Get Adaptive settings info
                [DetectorSettings, ~, ...
                    AdaptiveEmbeddedRuns_StimSettings] ...
                    ...
                    = ...
                    ...
                 createAdaptiveSettingsfromDeviceSettings(devicepath);

                % Look for adaptive Embedded Run
                if ~isempty(AdaptiveEmbeddedRuns_StimSettings)
                    db_out(d).adaptive_states=AdaptiveEmbeddedRuns_StimSettings.states(end);
                end

                db_out(d).adaptive_onset_dur = DetectorSettings.Ld0.onsetDuration;
                db_out(d).adaptive_termination_dur = DetectorSettings.Ld0.terminationDuration;
                db_out(d).adaptive_weights{1}(1:4) = cat(1,DetectorSettings.Ld0.features.weightVector);
                db_out(d).adaptive_pwrinputchan = DetectorSettings.Ld0.detectionInputs_BinaryCode  ;
                db_out(d).adaptive_threshold = DetectorSettings.Ld0.biasTerm;
                db_out(d).adaptive_updaterate = DetectorSettings.Ld0.updateRate;
            catch
            end



            try

                % Get stim settings
                [stimSettingsOut, stimMetaData] ...
                    ...
                    = ...
                    ...
                createStimSettingsFromDeviceSettings(devicepath);

                db_out(d).stim = stimSettingsOut.therapyStatus;

            catch
            end



            %Get stim information if STIM is on
            if sum(db_out(d).stim) > 0

                stimfile          = findFilesBVQX(sess_dirs{d},'StimLog.json');
                [stimpath,~,~]    = fileparts(stimfile{1});
                [stimLogSettings] = createStimSettingsTable(stimpath,stimMetaData);

                db_out(d).stimparams = stimLogSettings.stimParams_prog1;
                    
                stimnamegroup = {'A','B','C','D'; '1' , '5', '9','13'};
                [~,j]= find(contains(stimnamegroup,stimLogSettings.activeGroup));
                stimname =  metaData.stimProgramNames(str2double(stimnamegroup{2,j(1)}));
                db_out(d).stimName =  stimname{1};

                db_out(d).cyclingEnabled = stimSettingsOut.cyclingEnabled{1};
                % db_out(d).cycleOnOff     = stimSettingsOut.cycleOnOff{1} ;


                if height(stimLogSettings) > 1
                    db_out(d).stimparams = stimLogSettings;

                else
                 
                    stim_para = strsplit(char(db_out(d).stimparams),',');
                    
                    db_out(d).contacts = stim_para{1};
                    db_out(d).amp      = str2double(stim_para{2}(2:end -2));
                    db_out(d).PW       = str2double(stim_para{3}(2:end -2));
                    db_out(d).freq     = str2double(stim_para{4}(2:end -2));


                end
            end

            % load event file - not in use for now (PS)
            %             eventData = createEventLogTable(tdfile{1});
            %             dbout(d).eventData = eventData;


            % does mat file exist?
            matfile = findFilesBVQX(sess_dirs{d},'combinedDataTable.mat');

            if isempty(matfile) % no matlab data loaded
                db_out(d).matExist = false;
                %                 dbout(d).fnm = [];
            else
                db_out(d).matExist = true;
                %                 dbout(d).fnm = matfile{1};
            end
        end
    end
end

database_out = struct2table(db_out,'AsArray',true);

% delete all rows with empty session names ( WHY DOES THIS OCCUR?)
database_out = database_out(cellfun(@(x) ~isempty(x),database_out.sessname),:);

sorted_database = sortrows(database_out, 'sessname'); %sorting by session name
sorted_database.rec = (1:size(sorted_database,1))';
% 
% timeStart = [database_out.timeStart{:,1}];
% timeStop  = [database_out.timeStop{:,1}];

%% clear empty session rows and assign to new variable 'badsessions'
if iscell(sorted_database.timeStart)
    loc = cellfun('isempty', sorted_database{:,'timeStart'});
else
    loc= isempty(sorted_database.time);
end

badsessions = sorted_database(loc,:);
varargout{1} = badsessions;

sorted_database(loc,:) = [];

% format datetimes

% sorted_database = removevars(sorted_database, {'timeStart','timeStop'});
% 
% sorted_database.timeStart = timeStart';
% sorted_database.timeStop = timeStop';


%% expanding all fields within each struct

expanded_database = [];

for rowidx = 1:size(sorted_database, 1)
    tmp_row = sorted_database(rowidx,:);  %tmp_row is the row with multiple entries

    if size(tmp_row.timeStart{1}, 1) > 1  % duplicating entire row if there are multiple entries per session

        for new_row = 1:size(tmp_row.timeStart{1}, 1)
            expanded_database = [expanded_database; tmp_row];

            for col_name = ["time", "duration", "TDfs"]
                expanded_database{end, col_name}{1} = expanded_database{end, col_name}{1}(new_row);
            end

            %make the first subsession an integer (like 2), and  all subsessions
            %decimals like  2.01, 2.02, etc.
            if new_row ==1
                expanded_database.rec(end) = tmp_row.rec;
            else
                expanded_database.rec(end) = tmp_row.rec + ((new_row-1)/100);
            end

        end
    else  % print the single value  if only one entry per session

        expanded_database = [expanded_database; tmp_row];
        for col_name = ["timeStart", "timeStop", "duration", "TDfs"]
            expanded_database{end, col_name}(1) = expanded_database{end, col_name}(1);
        end
    end
end

% expand all variables for each row and make 'Disabled' values in TDfs to NaN
if ~isempty(expanded_database)

    idx_disabled = strcmp(expanded_database.TDfs,'Disabled');
    expanded_database.TDfs(idx_disabled)={nan};
    
    idx_emptyfft = cellfun(@isempty, expanded_database.fft);
    expanded_database.fft(idx_emptyfft)={nan};
    
    %  convert cells to string or double to remove cell structure
    cellvars = {'timeStart', 'timeStop', 'duration','battery','fft'};
    for n = 1:numel(cellvars)
    
        if n >= 4
            expanded_database.(cellvars{n}) = cell2mat(expanded_database.(cellvars{n}));
        else
            expanded_database.(cellvars{n}) = [expanded_database.(cellvars{n}){:}]';
        end
    
    end
    
    expanded_database = movevars(expanded_database, {'TDchan0', 'TDchan1', 'TDchan2', 'TDchan3'}, 'After', 'TDfs');
    RCSdatabase_out = table2timetable(expanded_database); % rename output for clarity

end




%% COMBINE WITH OLD DATABASE
% IF the old database existed, recombine with new database and sort it
% but first fix cell/ mat class issues

if ~isempty(old_database)
    disp('combining with old database...');

    if exist('RCSdatabase_out', 'var')

        %make cells to mat for some fields
        if iscell(RCSdatabase_out.matExist)
            % format some columns so they are not cells
            RCSdatabase_out.matExist = cell2mat(RCSdatabase_out.matExist);
            badsessions.matExist = cell2mat(badsessions.matExist);
        end
    
        if iscell(old_database.matExist)
            old_database.matExist = cell2mat(old_database.matExist);
            old_badsessions.matExist = cell2mat(old_badsessions.matExist);
        end
    
    
        if iscell(old_database.TDfs)
            idx_disabled=strcmp(old_database.TDfs,'Disabled');
            old_database.TDfs(idx_disabled)={nan};
    
    
            old_database.TDfs = old_database.TDfs;
    
        end
    
    
        if isa(RCSdatabase_out.adaptive_onset_dur,'double')
            RCSdatabase_out.adaptive_onset_dur =  num2cell(RCSdatabase_out.adaptive_onset_dur);
            RCSdatabase_out.adaptive_termination_dur =  num2cell(RCSdatabase_out.adaptive_termination_dur);
            RCSdatabase_out.adaptive_updaterate =  num2cell(RCSdatabase_out.adaptive_updaterate);
    
    
            badsessions.adaptive_onset_dur =  num2cell(badsessions.adaptive_onset_dur);
            badsessions.adaptive_termination_dur =  num2cell(badsessions.adaptive_termination_dur);
            badsessions.adaptive_updaterate =  num2cell(badsessions.adaptive_updaterate);
    
    
        end
    
    
        %     COMBINE HERE
        RCSdatabase_out.rec = RCSdatabase_out.rec + old_database.rec(end);
    
        new_database_out = [old_database;RCSdatabase_out];
    
        if ~isempty(badsessions)
            badsessions = [old_badsessions;badsessions];
        else
            badsessions = old_badsessions;
        end
    
        clear RCSdatabase_out
        RCSdatabase_out = new_database_out;  %already a timetable


       
    else
        RCSdatabase_out = old_database;
        varargout{1}     = old_badsessions;

        disp(['Old database as of ', datestr(old_database.time(end)), '; ',...
            old_database.sessname{end}])

    end

else

    % OUTPUTS =========================================================
    if nargout == 2
        varargout{1} = badsessions;
    end
    
    % Rename file to include patient ID
    writetimetable(RCSdatabase_out,fullfile(proc_dirname,[ptIDside...
        '_database.csv']));
    
    save(fullfile(proc_dirname,[ptIDside '_database.mat']),...
        'RCSdatabase_out','badsessions');
    
    fprintf('csv and mat of database saved as %s to %s \n',...
        [ptIDside '_database.mat'],proc_dirname);
    
    %==================================================================

end
end