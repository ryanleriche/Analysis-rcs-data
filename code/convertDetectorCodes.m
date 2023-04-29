function [convertedLd] = convertDetectorCodes(inputLd,FFTinterval)
%%
% Takes information from Ld0 and Ld1 fields of DeviceSettings.DetectionConfig
% and converts codes into human readable information
%%
% unchanged fields
convertedLd.biasTerm = inputLd.biasTerm;
convertedLd.features = inputLd.features;
convertedLd.fractionalFixedPointValue = inputLd.fractionalFixedPointValue;

convertedLd.updateRate                      = inputLd.updateRate;
convertedLd.blankingDurationUponStateChange = inputLd.blankingDurationUponStateChange;

convertedLd.onsetDuration                   = inputLd.onsetDuration;
convertedLd.holdoffTime                     = inputLd.holdoffTime;
convertedLd.terminationDuration             = inputLd.terminationDuration;

%% Fields which encode information per bit

input_binary_str =  dec2bin(inputLd.detectionInputs,8);


% Detection inputs
convertedLd.Inputs_BinaryCode = input_binary_str;

% Detection Enable
convertedLd.Enable_BinaryCode = dec2bin(inputLd.detectionEnable,8);

% from binary string return Ch_powerband as cell
i_ch_pb_in_LD    = regexp(reverse(input_binary_str), '1');

convertedLd.powerband_Inputs = cell(length(i_ch_pb_in_LD),1);

for i=1:length(i_ch_pb_in_LD)
    
    switch i_ch_pb_in_LD(i)

        case 1;     pb  = 'Ch0Band0';
        case 2;     pb  = 'Ch0Band1';

        case 3;     pb  = 'Ch1Band0';
        case 4;     pb  = 'Ch1Band1';

        case 5;     pb  = 'Ch2Band0';
        case 6;     pb  = 'Ch2Band1';

        case 7;     pb  = 'Ch3Band0';
        case 8;     pb  = 'C32Band1';   
    end
    convertedLd.powerband_Inputs{i} = pb;
end
end