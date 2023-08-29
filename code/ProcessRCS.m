    function [unifiedDerivedTimes,...
    timeDomainData, timeDomainData_onlyTimeVariables, timeDomain_timeVariableNames,...
    ...
    AccelData, AccelData_onlyTimeVariables, Accel_timeVariableNames,...
    ...
    PowerData, PowerData_onlyTimeVariables, Power_timeVariableNames,...
    ...
    FFTData, FFTData_onlyTimeVariables, FFT_timeVariableNames,...
    ...
    AdaptiveData, AdaptiveData_onlyTimeVariables, Adaptive_timeVariableNames,...
    ...
    timeDomainSettings, powerSettings, fftSettings, eventLogTable,...
    ...
    metaData, stimSettingsOut, stimMetaData, stimLogSettings,...
    ...
    DetectorSettings, AdaptiveStimSettings, AdaptiveEmbeddedRuns_StimSettings,...
    versionInfo] = ProcessRCS(cfg, varargin)
%%
% Demo wrapper script for importing raw .JSON files from RC+S, parsing
% into Matlab table format, and handling missing packets / harmonizing
% timestamps across data streams. Currently assuming you always have
% timeDomain data
%
% Depedencies:
% https://github.com/JimHokanson/turtle_json
% in the a folder called "toolboxes" in the same directory as the processing scripts
%
% Input =
% (1) RC+S Device folder, containing raw JSON files
% (2) Flag indicating if data should be saved (or read if already created):
%       1 = Process and save (overwrite if processed file already exist) -- DEFAULT
%       2 = Process and do not save
%       3 = If processed file already exists, then load. If it does not
%       exist, process and save
%       4 = If processed file already exists, then load. If it does not
%       exist, process but do not save
% (3) Advanced user: input '1' to use systemTick for time alignment of chunks if gaps
% < 6 seconds; default (input '0') will use PacketGenTime to align the start of all
% chunks regardless of preceeding gap length
%
% The raw data directory indicated or selected will be checked for the
% processed data file
%
% JSON files:
% AdaptiveLog.json
% DeviceSettings.json
% DiagnosticLogs.json
% ErrorLog.json
% EventLog.json
% RawDataAccel.json
% RawDataFFT.json
% RawDataPower.json
% RawDataTD.json
% StimLog.json
% TimeSync.json
%%
% Parse input variables, indicating folderPath and/or processFlag
switch nargin -1
    case 0
        folderPath = uigetdir();
        processFlag = 1;
        shortGaps_systemTick = 0;
    case 1
        if length(varargin{1}) == 1 % this indicates processFlag was input
            folderPath = uigetdir();
            processFlag = varargin{1};
        else
            folderPath  = varargin{1};
            processFlag = 1;
        end
        shortGaps_systemTick = 0;
    case 2
        folderPath  = varargin{1};
        processFlag = varargin{2};
        shortGaps_systemTick = 0;
    case 3
        if isempty(varargin{1})
            folderPath = uigetdir();
        else
            folderPath  = varargin{1};
        end
        processFlag = varargin{2};
        shortGaps_systemTick = varargin{3};
end

% Check if processed file exists
outputFileName = fullfile(folderPath,'AllDataTables.mat');
if processFlag == 3
    if isfile(outputFileName)
        disp('Loading previously processed file');
        load(outputFileName);
        processFlag = 0; % Do nothing else
    else
        processFlag = 1; % If no processed file, switch to process flag 1
    end
elseif processFlag == 4
    if isfile(outputFileName)
        disp('Loading previously processed file');
        load(outputFileName);
        processFlag = 0; % Do nothing else
    else
        processFlag = 2; % If no processed file, switch to process flag 12
    end
end

%%
if processFlag == 1 || processFlag == 2
    % DeviceSettings data
    if cfg.textoutputs;   disp('Collecting Device Settings data');   end
    DeviceSettings_fileToLoad = [folderPath filesep 'DeviceSettings.json'];
    if isfile(DeviceSettings_fileToLoad)
        [timeDomainSettings, powerSettings, fftSettings, metaData] = createDeviceSettingsTable(folderPath);
    else
        error('No DeviceSettings.json file')
    end
    %%
    % Stimulation settings
    if cfg.textoutputs;    disp('Collecting Stimulation Settings from Device Settings file');     end
    if isfile(DeviceSettings_fileToLoad)
        [stimSettingsOut, stimMetaData] = createStimSettingsFromDeviceSettings(folderPath);
    else
        warning('No DeviceSettings.json file - could not extract stimulation settings')
    end
    
    if cfg.textoutputs;        disp('Collecting Stimulation Settings from Stim Log file');        end
    StimLog_fileToLoad = [folderPath filesep 'StimLog.json'];
    if isfile(StimLog_fileToLoad)
        [stimLogSettings] = createStimSettingsTable(folderPath,stimMetaData);
    else
        warning('No StimLog.json file')
    end
    %%
    % Adaptive Settings
    if cfg.textoutputs;      disp('Collecting Adaptive Settings from Device Settings file');     end
    if isfile(DeviceSettings_fileToLoad)
        [DetectorSettings,AdaptiveStimSettings,AdaptiveEmbeddedRuns_StimSettings] = createAdaptiveSettingsfromDeviceSettings(folderPath);
    else
        error('No DeviceSettings.json file - could not extract detector and adaptive stimulation settings')
    end
    %%
    % Event Log
    if cfg.textoutputs;     disp('Collecting Event Information from Event Log file');       end
    EventLog_fileToLoad = [folderPath filesep 'EventLog.json'];
    
    if isfile(EventLog_fileToLoad)
        try
        [eventLogTable] = createEventLogTable(folderPath);
        catch
             eventLogTable = [];
              warning('No EventLog.json failed to load');

        end

    else
        warning('No EventLog.json file');
         eventLogTable = [];
    end
    
    %%
    % TimeDomain data
  
    TD_fileToLoad = [folderPath filesep 'RawDataTD.json'];
    if isfile(TD_fileToLoad)
        jsonobj_TD = deserializeJSON(TD_fileToLoad);
        if isfield(jsonobj_TD,'TimeDomainData') && ~isempty(jsonobj_TD.TimeDomainData)
             
            disp('Loading Time Domain Data');       
            [outtable_TD, srates_TD] = createTimeDomainTable(jsonobj_TD);
            
            disp('Creating derivedTimes for time domain:');     
            timeDomainData = assignTime(outtable_TD, shortGaps_systemTick);
        else
            timeDomainData = [];
        end
    else
        timeDomainData = [];
    end
    
    %%
    % Accelerometer data
    if cfg.textoutputs;    disp('Checking for Accelerometer Data');    end
    Accel_fileToLoad = [folderPath filesep 'RawDataAccel.json'];
    if isfile(Accel_fileToLoad)
        jsonobj_Accel = deserializeJSON(Accel_fileToLoad);
        if isfield(jsonobj_Accel,'AccelData') && ~isempty(jsonobj_Accel.AccelData)
            try
                if cfg.textoutputs;     disp('Loading Accelerometer Data');    end
                [outtable_Accel, srates_Accel] = createAccelTable(jsonobj_Accel);
                disp('Creating derivedTimes for accelerometer:')
                AccelData = assignTime(outtable_Accel, shortGaps_systemTick);
            catch
                warning('Accelerometer data present but error when trying to process. Will be excluded.')
                AccelData = [];
            end
        else
            AccelData = [];
        end
    else
        AccelData = [];
    end
    
    %%
    % Power data
    if cfg.textoutputs;      disp('Checking for Power Data');        end
    Power_fileToLoad = [folderPath filesep 'RawDataPower.json'];
    if isfile(Power_fileToLoad)
        if cfg.textoutputs;      disp('Loading Power Data');     end
        % Checking if power data is empty happens within createPowerTable
        % function
        [outtable_Power] = createPowerTable(cfg, folderPath);
        
        % Calculate power band cutoffs (in Hz) and add column to powerSettings
        if ~isempty(outtable_Power)
            try
                % Add samplerate and packetsizes column to outtable_Power -- samplerate is inverse
                % of fftConfig.interval
                numSettings = size(powerSettings,1);
                % Determine if more than one sampling rate across recording
                for iSetting = 1:numSettings
                    all_powerFs(iSetting) =  1/((powerSettings.fftConfig(iSetting).interval)/1000);
                end
                
                if length(unique(all_powerFs)) > 1
                    % Multiple sample rates for power data in the full file
                    PowerData = createDataTableWithMultipleSamplingRates(all_powerFs,powerSettings,outtable_Power,shortGaps_systemTick);
                else
                    % Same sample rate for power data for the full file
                    powerDomain_sampleRate = unique(all_powerFs);
                    outtable_Power.samplerate(:) = powerDomain_sampleRate;
                    outtable_Power.packetsizes(:) = 1;
                    PowerData = assignTime(outtable_Power, shortGaps_systemTick);
                end
            catch
                warning('Power data present but error when trying to process. Will be excluded.')
                PowerData = [];
            end
        else
            PowerData = [];
        end
    else
        PowerData = [];
    end
    
    %%
    % FFT data
    if cfg.textoutputs;      disp('Checking for FFT Data');      end
    FFT_fileToLoad = [folderPath filesep 'RawDataFFT.json'];
    if isfile(FFT_fileToLoad)
        jsonobj_FFT = deserializeJSON(FFT_fileToLoad);
        if isfield(jsonobj_FFT,'FftData') && ~isempty(jsonobj_FFT.FftData)
            try
                if cfg.textoutputs;     disp('Loading FFT Data');    end
                outtable_FFT = createFFTtable(jsonobj_FFT);
                
                % Add FFT parameter info to fftSettings
                numSettings = size(fftSettings,1);
                for iSetting = 1:numSettings
                    currentFFTconfig = fftSettings.fftConfig(iSetting);
                    currentTDsampleRate = fftSettings.TDsampleRates(iSetting);
                    fftParameters = getFFTparameters(currentFFTconfig,currentTDsampleRate);
                    fftSettings.fftParameters(iSetting) = fftParameters;
                end
                % Add samplerate and packetsizes column to outtable_FFT -- samplerate is inverse
                % of fftConfig.interval; in principle this interval could change
                % over the course of the recording
                
                % Determine if more than one sampling rate across recording
                for iSetting = 1:numSettings
                    all_fftFs(iSetting) =  1/((fftSettings.fftConfig(iSetting).interval)/1000);
                end
                
                if length(unique(all_fftFs)) > 1
                    FFTData = createDataTableWithMultipleSamplingRates(all_fftFs,fftSettings,outtable_FFT,shortGaps_systemTick);
                else
                    % Same sample rate for FFT data for the full file
                    FFT_sampleRate = unique(all_fftFs);
                    outtable_FFT.samplerate(:) = FFT_sampleRate;
                    outtable_FFT.packetsizes(:) = 1;
                    if cfg.textoutputs;      disp('Creating derivedTimes for FFT:');     end
                    FFTData = assignTime(outtable_FFT, shortGaps_systemTick);
                end
            catch
                warning('FFT data present but error when trying to process. Will be excluded.')
                FFTData = [];
            end
        else
            FFTData = [];
        end
    else
        FFTData = [];
    end
    %%
    % Adaptive data
    if cfg.textoutputs;       disp('Checking for Adaptive Data');     end
    Adaptive_fileToLoad = [folderPath filesep 'AdaptiveLog.json'];
    if isfile(Adaptive_fileToLoad)
        jsonobj_Adaptive = deserializeJSON(Adaptive_fileToLoad);
        if isfield(jsonobj_Adaptive,'AdaptiveUpdate') && ~isempty(jsonobj_Adaptive(1).AdaptiveUpdate)
            try
                if cfg.textoutputs;       disp('Loading Adaptive Data');      end
                outtable_Adaptive = createAdaptiveTable(jsonobj_Adaptive);
                % Note: StateTime must still be converted to sec in
                % outtable_Adaptive
                
                % Calculate adaptive_sampleRate - determine if more than one
                if size(fftSettings,1) == 1
                    adaptive_sampleRate =  1/((fftSettings.fftConfig(1).interval)/1000);
                    outtable_Adaptive.samplerate(:) = adaptive_sampleRate;
                    outtable_Adaptive.packetsizes(:) = 1;
                    outtable_Adaptive.StateTime = outtable_Adaptive.StateTime * (fftSettings.fftConfig(1).interval/1000);
                    
                    if cfg.textoutputs;      disp('Creating derivedTimes for Adaptive:');        end
                    AdaptiveData = assignTime(outtable_Adaptive, shortGaps_systemTick);
                else
                    for iSetting = 1:size(fftSettings,1)
                        all_adaptiveFs(iSetting) =  1/((fftSettings.fftConfig(iSetting).interval)/1000);
                    end
                    if length(unique(all_adaptiveFs)) > 1
                        AdaptiveData = createDataTableWithMultipleSamplingRates(all_adaptiveFs,fftSettings,outtable_Adaptive,shortGaps_systemTick);
                    else
                        adaptive_sampleRate = all_adaptiveFs(1);
                        outtable_Adaptive.samplerate(:) = adaptive_sampleRate;
                        outtable_Adaptive.packetsizes(:) = 1;
                        
                        if cfg.textoutputs;      disp('Creating derivedTimes for Adaptive:');        end
                        AdaptiveData = assignTime(outtable_Adaptive, shortGaps_systemTick);
                    end
                end
            catch
                 warning('Adaptive data present but error when trying to process. Will be excluded.')
                AdaptiveData = [];
            end
        else
            AdaptiveData = [];
        end
    else
        AdaptiveData = [];
    end
    
    %%
    % First, need to create unifiedDerivedTimes - which has DerivedTimes
    % filling in the gaps (even when there is no TD data)
    if ~isempty(timeDomainData)
        unifiedDerivedTimes = timeDomainData.DerivedTime(1):1000/srates_TD(1):timeDomainData.DerivedTime(end);
        unifiedDerivedTimes = unifiedDerivedTimes';

        % Time format for human-readable time
        timeFormat = 'America/Los_Angeles';

        % Harmonize Accel with unifiedDerivedTimes
        if ~isempty(AccelData)
            if cfg.textoutputs;      disp('Harmonizing time of Accelerometer data with unifiedDerivedTimes');        end
            derivedTime_Accel = AccelData.DerivedTime;
            [newDerivedTime,newDerivedTimes_Accel] = harmonizeTimeAcrossDataStreams(unifiedDerivedTimes, derivedTime_Accel, srates_TD(1));

            % Update unifiedDerivedTimes with newDerivedTime, as additional times
            % may have been added at the beginning and/or end
            unifiedDerivedTimes = newDerivedTime;

            % Add newDerivedTime to AccelData
            AccelData.newDerivedTime = newDerivedTimes_Accel;
            AccelData = movevars(AccelData, 'newDerivedTime','Before',1);

            % Add human-readable version of newDerivedTime to AccelData
            localTime = datetime(AccelData.newDerivedTime/1000,...
                'ConvertFrom','posixTime','TimeZone',timeFormat,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
            AccelData = addvars(AccelData,localTime,'Before',1);

            % Create sparse matrix with only timing variables
            subsetAccelData = AccelData(:,{'newDerivedTime','DerivedTime','timestamp','systemTick','PacketGenTime'});
            AccelData_onlyTimeVariables = table2array(subsetAccelData);
            AccelData_onlyTimeVariables(isnan(AccelData_onlyTimeVariables)) = 0;
            AccelData_onlyTimeVariables = sparse(AccelData_onlyTimeVariables);
            Accel_timeVariableNames = {'newDerivedTime','DerivedTime','timestamp','systemTick','PacketGenTime'};

            % Remove un-needed variables from AccelData
            AccelData = removevars(AccelData, {'DerivedTime','timestamp','systemTick','PacketGenTime',...
                'PacketRxUnixTime','dataTypeSequence','packetsizes'});
        else
            AccelData_onlyTimeVariables = [];
            Accel_timeVariableNames = [];
        end

        % Harmonize Power with unifiedDerivedTimes
        if ~isempty(PowerData)
            if cfg.textoutputs;      disp('Harmonizing time of Power data with unifiedDerivedTimes');        end
            derivedTime_Power = PowerData.DerivedTime;
            [newDerivedTime,newDerivedTimes_Power] = harmonizeTimeAcrossDataStreams(unifiedDerivedTimes, derivedTime_Power, srates_TD(1));

            % Update unifiedDerivedTimes with newDerivedTime, as additional times
            % may have been added at the beginning and/or end
            unifiedDerivedTimes = newDerivedTime;

            % Add newDerivedTime to PowerData
            PowerData.newDerivedTime = newDerivedTimes_Power;
            PowerData = movevars(PowerData, 'newDerivedTime','Before',1);

            % Add human-readable version of newDerivedTime to PowerData
            localTime = datetime(PowerData.newDerivedTime/1000,...
                'ConvertFrom','posixTime','TimeZone',timeFormat,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
            PowerData = addvars(PowerData,localTime,'Before',1);

            % Create sparse matrix with only timing variables
            subsetPowerData = PowerData(:,{'newDerivedTime','DerivedTime','timestamp','systemTick','PacketGenTime'});
            PowerData_onlyTimeVariables = table2array(subsetPowerData);
            PowerData_onlyTimeVariables(isnan(PowerData_onlyTimeVariables)) = 0;
            PowerData_onlyTimeVariables = sparse(PowerData_onlyTimeVariables);
            Power_timeVariableNames = {'newDerivedTime','DerivedTime','timestamp','systemTick','PacketGenTime'};

            % Remove un-needed variables from PowerData
            PowerData = removevars(PowerData, {'DerivedTime','timestamp','systemTick','PacketGenTime',...
                'PacketRxUnixTime','dataTypeSequence','packetsizes'});
        else
            PowerData_onlyTimeVariables = [];
            Power_timeVariableNames = [];
        end

        % Harmonize FFT with unifiedDerivedTimes
        if ~isempty(FFTData)
           if cfg.textoutputs;        disp('Harmonizing time of FFT data with unifiedDerivedTimes');     end
            derivedTime_FFT = FFTData.DerivedTime;
            [newDerivedTime,newDerivedTimes_FFT] = harmonizeTimeAcrossDataStreams(unifiedDerivedTimes, derivedTime_FFT, srates_TD(1));

            % Update unifiedDerivedTimes with newDerivedTime, as additional times
            % may have been added at the beginning and/or end
            unifiedDerivedTimes = newDerivedTime;

            FFTData.newDerivedTime = newDerivedTimes_FFT;
            FFTData = movevars(FFTData, 'newDerivedTime','Before',1);

            % Add human-readable version of newDerivedTime to FFTData
            localTime = datetime(FFTData.newDerivedTime/1000,...
                'ConvertFrom','posixTime','TimeZone',timeFormat,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
            FFTData = addvars(FFTData,localTime,'Before',1);

            % Create sparse matrix with only timing variables
            subsetFFTData = FFTData(:,{'newDerivedTime','DerivedTime','timestamp','systemTick','PacketGenTime'});
            FFTData_onlyTimeVariables = table2array(subsetFFTData);
            FFTData_onlyTimeVariables(isnan(FFTData_onlyTimeVariables)) = 0;
            FFTData_onlyTimeVariables = sparse(FFTData_onlyTimeVariables);
            FFT_timeVariableNames = {'newDerivedTime','DerivedTime','timestamp','systemTick','PacketGenTime'};

            % Remove un-needed variables from FFTData
            FFTData = removevars(FFTData, {'DerivedTime','timestamp','systemTick','PacketGenTime',...
                'PacketRxUnixTime','dataTypeSequence','packetsizes'});
        else
            FFTData_onlyTimeVariables = [];
            FFT_timeVariableNames = [];
        end

        % Harmonize Adaptive with unifiedDerivedTimes
        if ~isempty(AdaptiveData)
            if cfg.textoutputs;      disp('Harmonizing time of Adaptive data with unifiedDerivedTimes');     end
            derivedTime_Adaptive = AdaptiveData.DerivedTime;
            [newDerivedTime,newDerivedTimes_Adaptive] = harmonizeTimeAcrossDataStreams(unifiedDerivedTimes, derivedTime_Adaptive, srates_TD(1));

            % Update unifiedDerivedTimes with newDerivedTime, as additional times
            % may have been added at the beginning and/or end
            unifiedDerivedTimes = newDerivedTime;

            AdaptiveData.newDerivedTime = newDerivedTimes_Adaptive;
            AdaptiveData = movevars(AdaptiveData, 'newDerivedTime','Before',1);

            % Add human-readable version of newDerivedTime to AdaptiveData
            localTime = datetime(AdaptiveData.newDerivedTime/1000,...
                'ConvertFrom','posixTime','TimeZone',timeFormat,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
            AdaptiveData = addvars(AdaptiveData,localTime,'Before',1);

            % Create sparse matrix with only timing variables
            subsetAdaptiveData = AdaptiveData(:,{'newDerivedTime','DerivedTime','timestamp','systemTick','PacketGenTime'});
            AdaptiveData_onlyTimeVariables = table2array(subsetAdaptiveData);
            AdaptiveData_onlyTimeVariables(isnan(AdaptiveData_onlyTimeVariables)) = 0;
            AdaptiveData_onlyTimeVariables = sparse(AdaptiveData_onlyTimeVariables);
            Adaptive_timeVariableNames = {'newDerivedTime','DerivedTime','timestamp','systemTick','PacketGenTime'};

            % Remove un-needed variables from AdaptiveData
            AdaptiveData = removevars(AdaptiveData, {'DerivedTime','timestamp','systemTick','PacketGenTime',...
                'PacketRxUnixTime','dataTypeSequence','packetsizes'});
        else
            AdaptiveData_onlyTimeVariables = [];
            Adaptive_timeVariableNames = [];
        end

        %%
        % Separate timeDomainData for saving

        % Add human-readable version of DerivedTime to timeDomainData
        localTime = datetime(timeDomainData.DerivedTime/1000,...
            'ConvertFrom','posixTime','TimeZone',timeFormat,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        timeDomainData = addvars(timeDomainData,localTime,'Before',1);

        % Create sparse matrix with only timing variables
        subsetTimeDomainData = timeDomainData(:,{'DerivedTime','timestamp','systemTick','PacketGenTime'});
        timeDomainData_onlyTimeVariables = table2array(subsetTimeDomainData);
        timeDomainData_onlyTimeVariables(isnan(timeDomainData_onlyTimeVariables)) = 0;
        timeDomainData_onlyTimeVariables = sparse(timeDomainData_onlyTimeVariables);
        timeDomain_timeVariableNames = {'DerivedTime','timestamp','systemTick','PacketGenTime'};

        % Remove un-needed variables from TimeDomainData
        timeDomainData = removevars(timeDomainData, {'timestamp','systemTick','PacketGenTime',...
            'PacketRxUnixTime','dataTypeSequence','packetsizes'});

        %%
        % Load txt file with code version

        currentFunctionPath = mfilename('fullpath');
        [filepath, name, ext] = fileparts(currentFunctionPath);
        toLoad = [filepath filesep 'VersionInfo.txt'];

        versionInfo = importdata(toLoad);

        %%
        % Save output file if indicated
        if processFlag == 1
            disp('Saving output')

            save(outputFileName,'unifiedDerivedTimes',...
                'timeDomainData','timeDomainData_onlyTimeVariables','timeDomain_timeVariableNames',...
                'AccelData','AccelData_onlyTimeVariables','Accel_timeVariableNames',...
                'PowerData','PowerData_onlyTimeVariables','Power_timeVariableNames',...
                'FFTData','FFTData_onlyTimeVariables','FFT_timeVariableNames',...
                'AdaptiveData','AdaptiveData_onlyTimeVariables','Adaptive_timeVariableNames',...
                'timeDomainSettings','powerSettings','fftSettings','eventLogTable','metaData',...
                'stimSettingsOut','stimMetaData','stimLogSettings','DetectorSettings','AdaptiveStimSettings',...
                'AdaptiveEmbeddedRuns_StimSettings','versionInfo','-v7.3');
        end
    else
        warning('Time domain file is empty')
        unifiedDerivedTimes=[];
        timeDomainData_onlyTimeVariables=[];
        timeDomain_timeVariableNames=[];
        AccelData_onlyTimeVariables=[];
        Accel_timeVariableNames=[];
        PowerData_onlyTimeVariables=[];
        Power_timeVariableNames=[];
        FFTData_onlyTimeVariables=[];
        FFT_timeVariableNames=[];
        AdaptiveData_onlyTimeVariables=[];
        Adaptive_timeVariableNames=[];
        versionInfo=[];
    end
end


end

