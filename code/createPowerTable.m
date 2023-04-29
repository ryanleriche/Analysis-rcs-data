function [powerTable]  = createPowerTable(cfg, folderPath)
%%
% Function to unravel Power data
% Input:
%   folderPath: path to the Device* folder, which contains all the .json
%   files for this recording
%
% Output: powerTable (power data in table format, without band limits
% decoded)
%%
try
    % Load power data
    rawPowerData = jsondecode(fixMalformedJson(fileread([folderPath filesep 'RawDataPower.json']),'EventLog'));
    
    % Initalize output
    powerTable = table();
    
    % If no power data, return empty tables, otherwise start parsing
    if isempty(rawPowerData) || isempty(rawPowerData.PowerDomainData)
         if cfg.textoutputs;   fprintf('Power data  is empty\n');   end
        fprintf('Creating dummy event table\n');
        powerTable  = [];
    else
        % Parsing data contained in headers
        Header = [rawPowerData.PowerDomainData.Header];
        variableNames = {'dataTypeSequence','systemTick'};
        powerData = struct();
        for iVariable = 1:length(variableNames)
            powerData.(variableNames{iVariable}) = [Header.(variableNames{iVariable})]';
        end
        
        % Parsing timestamp (stored inside struct)
        timestamps = [Header.timestamp];
        powerData.timestamp =  struct2array(timestamps)';
        
        % Converting TD samplerate
        powerData.TDsamplerate = getSampleRate([rawPowerData.PowerDomainData.SampleRate]');
        
        % Parsing data conatined in rawPowerData.PowerDomainData
        PowerDomainData = [rawPowerData.PowerDomainData];
        variableNames = {'PacketGenTime','PacketRxUnixTime','IsPowerChannelOverrange'};
        for iVariable = 1:length(variableNames)
            powerData.(variableNames{iVariable}) = [PowerDomainData.(variableNames{iVariable})]';
        end
        
        % Convert FftSize to human-readable
        temp_FftSize = [PowerDomainData.FftSize]';
        temp_FftSize(temp_FftSize == 0) = 64;
        temp_FftSize(temp_FftSize == 1) = 256;
        temp_FftSize(temp_FftSize == 3) = 1024;     
        powerData.('FftSize') = temp_FftSize;
        
        % Convert ValidDataMask to bit string
        temp_ValidDataMask = [PowerDomainData.ValidDataMask]';
        powerData.('ValidDataMask') = cellstr(dec2bin(temp_ValidDataMask,8));
        
         % Convert ExternalValuesMask to bit string
        temp_ExternalValuesMask = [PowerDomainData.ExternalValuesMask]';
        powerData.('ExternalValuesMask') = cellstr(dec2bin(temp_ExternalValuesMask,8));
        
        % Parsing data conatined in Bands
        bands = [PowerDomainData.Bands]';
        for iBand = 1:size(bands,2)
            bandName = sprintf('Band%d',iBand);
            powerData.(bandName) = bands(:,iBand);
        end
        
        % If no power data streamed for a band, change zeros to Nans
        if sum(powerData.Band1) == 0
            powerData.Band1 = NaN(size(powerData.Band1));
        end
        if sum(powerData.Band2) == 0
            powerData.Band2 = NaN(size(powerData.Band2));
        end
        if sum(powerData.Band3) == 0
            powerData.Band3 = NaN(size(powerData.Band3));
        end
        if sum(powerData.Band4) == 0
            powerData.Band4 = NaN(size(powerData.Band4));
        end
        if sum(powerData.Band5) == 0
            powerData.Band5 = NaN(size(powerData.Band5));
        end
        if sum(powerData.Band6) == 0
            powerData.Band6 = NaN(size(powerData.Band6));
        end
        if sum(powerData.Band7) == 0
            powerData.Band7 = NaN(size(powerData.Band7));
        end
        if sum(powerData.Band8) == 0
            powerData.Band8 = NaN(size(powerData.Band8));
        end
        
        powerTable = struct2table(powerData);
    end
catch
    powerTable = table();
end
end

