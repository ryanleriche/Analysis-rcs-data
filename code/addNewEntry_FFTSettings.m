function [newEntry,fftConfig] = addNewEntry_FFTSettings(actionType,recNum,currentSettings,TDsettings,fftConfig)
%%
% Extract FFT settings data in order to add a new row to the
% FFT_SettingsTable; ensures that all table fields are filled for each
% entry (otherwise warning will print)
%
%%
HostUnixTime = currentSettings.RecordInfo.HostUnixTime;

newEntry.action = actionType;
newEntry.recNum = recNum;
newEntry.time   = HostUnixTime;

% parse and save FFT configuration
if isfield(currentSettings,'SensingConfig') 
    
    if isfield(currentSettings.SensingConfig,'fftConfig')
        fftConfig   = convertFFTCodes(currentSettings.SensingConfig.fftConfig);
        
    end

    if isfield(currentSettings,'StreamState')
        % return streamed FFT channel or if it's disabled
        if  currentSettings.StreamState.FftStreamEnabled == 1
            fftConfig.fftStreamChannel      = num2str(currentSettings.SenseState.fftStreamChannel);
        
        elseif currentSettings.StreamState.FftStreamEnabled == 0
            fftConfig.fftStreamChannel      = 'Disabled';
        end

    else
        fftConfig.fftStreamChannel      = 'No Stream State field';

    end
end

newEntry.fftConfig = fftConfig;

% Get sample rate for each TD channel; all TD channels have
% same Fs (or is listed as NaN)
TDsampleRates = NaN(1,4);
for iChan = 1:4
    if isnumeric(TDsettings(iChan).sampleRate)
        TDsampleRates(iChan) = TDsettings(iChan).sampleRate;
    end
end
TDsampleRates = unique(TDsampleRates);
currentTDsampleRate = TDsampleRates(~isnan(TDsampleRates));
newEntry.TDsampleRates = currentTDsampleRate;

end