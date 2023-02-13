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
convertedLd.terminationDuration             =  inputLd.terminationDuration;

%% Fields which encode information per bit
% Detection inputs
convertedLd.detectionInputs_BinaryCode = dec2bin(inputLd.detectionInputs,8);

% Detection Enable
convertedLd.detectionEnable_BinaryCode = dec2bin(inputLd.detectionEnable,8);

convertedLd.powerband_detectionInputs = cell(length(inputLd.detectionInputs),1);

for i=1:length(inputLd.detectionInputs)

    switch inputLd.detectionInputs(i)
        case 0;   pb  = 'None';

        case 1;   pb  = 'Ch0Band0';
        case 2;   pb  = 'Ch0Band1';

        case 4;   pb  = 'Ch1Band0';
        case 8;   pb  = 'Ch1Band1';

        case 16;  pb  = 'Ch2Band0';
        case 32;  pb  = 'Ch2Band1';

        case 64;  pb  = 'Ch3Band0';
        case 128; pb = 'C32Band1';   
    end

    convertedLd.powerband_detectionInputs{i} = pb;

end
end