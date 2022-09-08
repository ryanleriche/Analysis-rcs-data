function [stimSettingsOut, stimMetaData] = createStimSettingsFromDeviceSettings(folderPath)
%%
% Extract information pertaining to stimulation from DeviceSettings.json
%
% Input: Folder path to *Device folder containing json files
%
% Output: stimSettingsOut table, stimMetaData containing info about
% groups/programs which were active and which contacts were used for stim
%
%%
stimSettingsOut = table();

%%
DeviceSettings = deserializeJSON([folderPath filesep 'DeviceSettings.json']);
%%
% Fix format - Sometimes device settings is a struct or cell array
if isstruct(DeviceSettings)
    DeviceSettings = {DeviceSettings};
end

%%
% Get enabled programs from first record; isEnabled = 0 means program is
% enabled; isEnabled = 131 means program is disabled. These are static for
% the full recording. Get contact information for these programs

currentSettings = DeviceSettings{1};

HostUnixTime = currentSettings.RecordInfo.HostUnixTime;

% initialize cathodes and anodes for all programs
stimMetaData.anodes = cell(4,4);
stimMetaData.cathodes = cell(4,4);

counter = 1;
for iGroup = 0:3
    printGroupName = sprintf('TherapyConfigGroup%d',iGroup);
    for iProgram = 1:4
        temp = currentSettings.(printGroupName).programs(iProgram).isEnabled;
        if temp == 0
            stimMetaData.validPrograms(iGroup + 1,iProgram) = 1;
            
            switch iGroup
                case 0
                    currentGroupName = 'GroupA';
                case 1
                    currentGroupName = 'GroupB';
                case 2
                    currentGroupName = 'GroupC';
                case 3
                    currentGroupName = 'GroupD';
            end
            stimMetaData.validProgramNames{counter,1} = [currentGroupName '_program' num2str(iProgram)];
            
            rawElectrodeTable = currentSettings.(printGroupName).programs(iProgram).electrodes.electrodes;
            % Find electrode(s) which are enabled; record if they are anode
            % (1) or cathode (0)
            temp_anode = [];
            temp_cathode = [];
            for iElectrode = 1:length(rawElectrodeTable)
                isOff = rawElectrodeTable(iElectrode).isOff;
                if isOff == 0 % Indicates channels is active
                    
                    % Subtract one to get electrode contact, because zero indexed
                    if rawElectrodeTable(iElectrode).electrodeType == 1 % anode
                        temp_anode = [temp_anode iElectrode - 1];
                    elseif rawElectrodeTable(iElectrode).electrodeType == 0 % cathode
                        temp_cathode = [temp_cathode iElectrode - 1];
                    end
                    
                end
            end
            
            % Contact 16 indicates can
            stimMetaData.anodes{iGroup + 1,iProgram} = temp_anode;
            stimMetaData.cathodes{iGroup + 1,iProgram} = temp_cathode;
            
            counter = counter + 1;
        else
            stimMetaData.validPrograms(iGroup+1,iProgram) = 0;
        end
    end
end

%%
% Set up output table (stimSettingsOut) with initial settings
entryNumber = 1;

switch currentSettings.GeneralData.therapyStatusData.activeGroup
    case 0
        activeGroup = 'A';
    case 1
        activeGroup = 'B';
    case 2
        activeGroup = 'C';
    case 3
        activeGroup = 'D';
end
therapyStatus = currentSettings.GeneralData.therapyStatusData.therapyStatus;

stimSettingsOut.HostUnixTime(entryNumber)             = HostUnixTime;
stimSettingsOut.activeGroup{entryNumber}              = activeGroup;
stimSettingsOut.therapyStatus(entryNumber)            = therapyStatus;
stimSettingsOut.therapyStatusDescription{entryNumber} = convertTherapyStatus(therapyStatus);


stimSettingsOut.cyclingEnabled{entryNumber}           = currentSettings.TherapyConfigGroup0.cyclingEnabled;


%% see'Medtronic.NeuroStim.Olympus.DataTypes.Therapy.CyclingUnits' w/n Summit API HTML 

% cycling unit "0" means in 100 ms
if currentSettings.TherapyConfigGroup0.cycleOnTime.units == 0

    stimSettingsOut.cycleOnTime{entryNumber} = currentSettings.TherapyConfigGroup0.cycleOnTime.time /10;
    

% cycling unit "1" means in 1s
elseif currentSettings.TherapyConfigGroup0.cycleOnTime.units == 1

    stimSettingsOut.cycleOnTime{entryNumber} = currentSettings.TherapyConfigGroup0.cycleOnTime.time;

% cycling unit "2" means in 10s
elseif currentSettings.TherapyConfigGroup0.cycleOnTime.units == 2

    stimSettingsOut.cycleOnTime{entryNumber} = currentSettings.TherapyConfigGroup0.cycleOnTime.time * 10;
end


% repeat for cycleOff
if currentSettings.TherapyConfigGroup0.cycleOffTime.units == 0

    stimSettingsOut.cycleOffTime{entryNumber} = currentSettings.TherapyConfigGroup0.cycleOffTime.time /10;
    

% cycling unit "1" means in 1s
elseif currentSettings.TherapyConfigGroup0.cycleOffTime.units == 1

    stimSettingsOut.cycleOffTime{entryNumber} = currentSettings.TherapyConfigGroup0.cycleOffTime.time;

% cycling unit "2" means in 10s
elseif currentSettings.TherapyConfigGroup0.cycleOffTime.units == 2

    stimSettingsOut.cycleOffTime{entryNumber} = currentSettings.TherapyConfigGroup0.cycleOffTime.time * 10;
end
 
% rampTime LSB is 100 ms
% (Medtronic.NeuroStim.Olympus.DataTypes.Therapy.TherapyGroup.RampTime)
stimSettingsOut.rampTime{entryNumber} = currentSettings.TherapyConfigGroup0.rampTime /10;

%%
previousActiveGroup = activeGroup;
previousTherapyStatus = therapyStatus;

updateActiveGroup = 0;
updateTherapyStatus = 0;
% Determine if activeGroup and/or therapyStatus has changed
for iRecord = 1 : length(DeviceSettings)
    
    currentSettings = DeviceSettings{iRecord};
    HostUnixTime = currentSettings.RecordInfo.HostUnixTime;
    
    % Check if activeGroup has changed
    if isfield(currentSettings,'GeneralData') && isfield(currentSettings.GeneralData, 'therapyStatusData') &&...
            isfield(currentSettings.GeneralData.therapyStatusData, 'activeGroup')
        switch currentSettings.GeneralData.therapyStatusData.activeGroup
            case 0
                activeGroup = 'A';
            case 1
                activeGroup = 'B';
            case 2
                activeGroup = 'C';
            case 3
                activeGroup = 'D';
        end
        if ~isequal(activeGroup, previousActiveGroup)
            updateActiveGroup = 1;
        end
    end
    
    % Check if therapyStatus has changed (turned on/off)
    if isfield(currentSettings,'GeneralData') && isfield(currentSettings.GeneralData, 'therapyStatusData') &&...
            isfield(currentSettings.GeneralData.therapyStatusData, 'therapyStatus')
        
        therapyStatus = currentSettings.GeneralData.therapyStatusData.therapyStatus;
        if ~isequal(therapyStatus, previousTherapyStatus)
            updateTherapyStatus = 1;
        end
    end
    
    % If either activeGroup or therapyStatus has changed, add row to
    % output table
    if updateActiveGroup || updateTherapyStatus
        % Update table if either activeGroup or therapyStatus has changed
        toAdd.HostUnixTime = HostUnixTime;
        toAdd.activeGroup = activeGroup;
        toAdd.therapyStatus = therapyStatus;
        toAdd.therapyStatusDescription = convertTherapyStatus(therapyStatus);

        toAdd.cyclingEnabled = currentSettings.TherapyConfigGroup0.cyclingEnabled;
        toAdd.cycleOnOff    =[currentSettings.TherapyConfigGroup0.cycleOnTime.time,...
                                  currentSettings.TherapyConfigGroup0.cycleOffTime.time];



        stimSettingsOut = [stimSettingsOut; struct2table(toAdd)];
        
        clear toAdd
        % Update for next loop
        previousActiveGroup = activeGroup;
        previousTherapyStatus = therapyStatus;
        
        % Reset flags
        updateActiveGroup = 0;
        updateTherapyStatus = 0;



      
    end
end
end
