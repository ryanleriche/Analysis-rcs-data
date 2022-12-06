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
warning('off','all');

stimSettingsOut = table();


DeviceSettings = deserializeJSON([folderPath filesep 'DeviceSettings.json']);
%%
% Fix format - Sometimes device settings is a struct or cell array
if isstruct(DeviceSettings)
    DeviceSettings = {DeviceSettings};
end

currentSettings = DeviceSettings{1};

HostUnixTime    = currentSettings.RecordInfo.HostUnixTime;

%%
% Set up output table (stimSettingsOut) with initial settings
entryNumber       = 1;
activeGroup_num   = currentSettings.GeneralData.therapyStatusData.activeGroup;

switch activeGroup_num 
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

%%
% Get enabled programs from first record; isEnabled = 0 means program is
% enabled; isEnabled = 131 means program is disabled. These are static for
% the full recording. Get contact information for these programs

% initialize cathodes and anodes for all programs
stimMetaData.anodes = cell(4,4);
stimMetaData.cathodes = cell(4,4);

counter = 1;
for iGroup = 0:3

    printGroupName = ['TherapyConfigGroup', num2str(iGroup)];

    temp_group                = currentSettings.(printGroupName);

    group                     = struct();

    group.rateInHz            = temp_group.rateInHz;
    group.cyclingEnabled      = currentSettings.(printGroupName).cyclingEnabled;

    % four programs per group
    group.ampInMilliamps                 = [temp_group.programs.amplitudeInMilliamps];
    group.pulseWidthInMicroseconds       = [temp_group.programs.pulseWidthInMicroseconds]; 

    temp                      = [temp_group.programs.miscSettings];
    group.actRechRatio        = [temp.activeRechargeRatio];


    % see'Medtronic.NeuroStim.Olympus.DataTypes.Therapy.CyclingUnits' w/n Summit API HTML 
    if group.cyclingEnabled == 1 
        cycleOnTime = temp_group.cycleOnTime.time;
    
        switch temp_group.cycleOnTime.units
            % cycling unit "0" means LBS of 100 ms
            case 0
                group.cycleOnInSecs =  cycleOnTime/10;
        
            % cycling unit "1" means LBS of 1s
            case 1
                group.cycleOnInSecs = cycleOnTime;
        
            % cycling unit "2" means LBS of 10s
            case 2
                group.cycleOnInSecs = cycleOnTime * 10;
        end
    
        cycleOffTime = temp_group.cycleOffTime.time;
    
        switch temp_group.cycleOffTime.units
            % cycling unit "0" means in 100 ms
            case 0
                group.cycleOffInSecs = cycleOffTime /10;
        
            % cycling unit "1" means in 1s
            case 1
                group.cycleOffInSecs = cycleOffTime;
        
            % cycling unit "2" means in 10s
            case 2
                group.cycleOffInSecs = cycleOffTime * 10;
        end
    
    else
        group.cycleOnInSecs  = NaN;
        group.cycleOffInSecs = NaN;
    end
    
    % rampTime LSB is 100 ms
    % (Medtronic.NeuroStim.Olympus.DataTypes.Therapy.TherapyGroup.RampTime)
    group.rampInSecs  = temp_group.rampTime /10;

    group.rampRepeat = temp_group.rampRepeat;

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

    group.validPrograms     = [];

    group.validProgramNames  = cell(4,1);
    group.contacts           = table();

    for iProgram = 1:4

        temp = currentSettings.(printGroupName).programs(iProgram).isEnabled;

        if temp == 0

            stimMetaData.validPrograms(iGroup + 1,iProgram) = 1;
            group.validPrograms(iProgram)                   = 1;

            stimMetaData.validProgramNames{counter,1} = [currentGroupName '_program' num2str(iProgram)];
            group.validProgramNames{iProgram}         = [currentGroupName '_program' num2str(iProgram)];

            rawElectrodeTable = currentSettings.(printGroupName).programs(iProgram).electrodes.electrodes;
            % Find electrode(s) which are enabled; record if they are anode
            % (1) or cathode (0)
            temp_anode        = [];
            temp_cathode      = [];

            for iElectrode = 1 : length(rawElectrodeTable)
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
            stimMetaData.anodes{iGroup + 1,iProgram}   = temp_anode;
            stimMetaData.cathodes{iGroup + 1,iProgram} = temp_cathode;
            
            group.contacts.anodes{iProgram}            = temp_anode;
            group.contacts.cathodes{iProgram}          = temp_cathode;

            counter = counter + 1;
        else
            stimMetaData.validPrograms(iGroup+1,iProgram) = 0;
            group.validPrograms(iProgram)                 = 0;


            group.validProgramNames{iProgram}       = '';
            group.contacts.anodes{iProgram}         = '';
            group.contacts.cathodes{iProgram}       = '';

        end
    end

    stimSettingsOut.(currentGroupName)(entryNumber) = group;
end



%%
previousActiveGroup   = activeGroup;
previousTherapyStatus = therapyStatus;

updateActiveGroup     = 0;
updateTherapyStatus   = 0;
% Determine if activeGroup and/or therapyStatus has changed
for iRecord = 1 : length(DeviceSettings)
    
    currentSettings = DeviceSettings{iRecord};
    HostUnixTime    = currentSettings.RecordInfo.HostUnixTime;
    
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

        
        for iGroup = 0:3
    
            printGroupName = ['TherapyConfigGroup', num2str(iGroup)];
        
            temp_group                = currentSettings.(printGroupName);
        
            group                     = struct();
        
            group.rateInHz             = temp_group.rateInHz;
            group.cyclingEnabled       = temp_group.cyclingEnabled;
        
            % four programs per group
            group.ampInMilliamps             = [temp_group.programs.amplitudeInMilliamps];
            group.pulseWidthInMicroseconds   = [temp_group.programs.pulseWidthInMicroseconds]; 
        
            temp                      = [temp_group.programs.miscSettings];
            group.actRechRatio        = [temp.activeRechargeRatio];
        
        
            % see'Medtronic.NeuroStim.Olympus.DataTypes.Therapy.CyclingUnits' w/n Summit API HTML 
            if group.cyclingEnabled == 1 
                cycleOnTime = temp_group.cycleOnTime.time;
            
                switch temp_group.cycleOnTime.units
                    % cycling unit "0" means LBS of 100 ms
                    case 0
                        group.cycleOnInSecs =  cycleOnTime/10;
                
                    % cycling unit "1" means LBS of 1s
                    case 1
                        group.cycleOnInSecs = cycleOnTime;
                
                    % cycling unit "2" means LBS of 10s
                    case 2
                        group.cycleOnInSecs = cycleOnTime * 10;
                end
            
                cycleOffTime = temp_group.cycleOffTime.time;
            
                switch temp_group.cycleOffTime.units
                    % cycling unit "0" means in 100 ms
                    case 0
                        group.cycleOffInSecs = cycleOffTime /10;
                
                    % cycling unit "1" means in 1s
                    case 1
                        group.cycleOffInSecs = cycleOffTime;
                
                    % cycling unit "2" means in 10s
                    case 2
                        group.cycleOffInSecs = cycleOffTime * 10;
                end
            
            else
                group.cycleOnInSecs  = NaN;
                group.cycleOffInSecs = NaN;
            end
            
            % rampTime LSB is 100 ms
            % (Medtronic.NeuroStim.Olympus.DataTypes.Therapy.TherapyGroup.RampTime)
            group.rampInSecs  = temp_group.rampTime /10;
        
            group.rampRepeat = temp_group.rampRepeat;

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
        
            group.validPrograms     = [];
        
            group.validProgramNames  = cell(4,1);
            group.contacts           = table();
        
            for iProgram = 1:4
        
                temp = temp_group.programs(iProgram).isEnabled;
        
                if temp == 0
        
                    group.validPrograms(iProgram)                   = 1;
                    group.validProgramNames{iProgram}         = [currentGroupName '_program' num2str(iProgram)];
        
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
                    
                    group.contacts.anodes{iProgram}            = temp_anode;
                    group.contacts.cathodes{iProgram}          = temp_cathode;
        
                else
                    group.validPrograms(iProgram)                 = 0;
        
        
                    group.validProgramNames{iProgram}       = '';
                    group.contacts.anodes{iProgram}         = '';
                    group.contacts.cathodes{iProgram}       = '';
        
                end
            end
        
            toAdd.(currentGroupName)(entryNumber) = group;
        end

        stimSettingsOut = [stimSettingsOut; struct2table(toAdd)];
        
        clear toAdd
        % Update for next loop
        previousActiveGroup   = activeGroup;
        previousTherapyStatus = therapyStatus;
        
        % Reset flags
        updateActiveGroup = 0;
        updateTherapyStatus = 0;
      
    end
end
end
