function [newPower, newSettings] ...
    =...
    calculateNewPower_RBL(...
    ...
    combinedDataTable, fftSettings, powerSettings, metaData, channel, freqBand)

% calculates a new fft band power series using the streamed time domain signal and based on new power and fft settings
%
% Assumption: only one set of fftSettings and powerSettings each time this function is invoked.
% (1) Either the default (initial) settings of the recording session, or (see DEMO_CalculatePowerRCS.m, line~44)
% (2) The fftSettings and powerSettings passed via the user (see DEMO_CalculatePowerRCS.m, line~100)
% 
% Input = 
% (1) combinedDataTable
% (2) fftSettings
% (3) powerSettings 
% (4) metaData (type = output from DEMO_Process)
% (5) channel (type = integer (1..4), eg usage, channel = 1)
% (6) freqBand (type = integer array, eg usage, freqBand = [20 25])
% 
% If you find errors while using this code or want to help further develop
% it, contact juan.ansoromeo@ucsf.edu or juan.anso@gmail.com

% Parse input variables
newSettings.metaData   = metaData;
newSettings.tdChannel  = channel;
newSettings.bandLimits = freqBand;

% initialize with the default sesttings to have access to default settings
% relying only on first set of power and fft settings (removing all other settings changes of the session)
powerSettings(2:end,:) = []; 
newSettings.powerSettings = powerSettings;
fftSettings(2:end,:) = [];
newSettings.fftSettings = fftSettings;

% read actual amplifier gains per channel from metadata
ampGains = newSettings.metaData.ampGains;

% initialize output power table
newPower = table();
newPower.localTime       = combinedDataTable.localTime;
newPower.DerivedTimes    = combinedDataTable.DerivedTime;
newPower.calculatedPower = nan(1,size(combinedDataTable,1))';
newPower.calculatedFFT   = cell(1,size(combinedDataTable,1))';

% initialize a newPowerSettings array
newSettings.powerSettings.powerBands.indices_BandStart_BandStop(:,:) = []; 
newSettings.powerSettings.powerBands.powerBandsInHz = [];
newSettings.powerSettings.powerBands.lowerBound = [];
newSettings.powerSettings.powerBands.upperBound = [];
newSettings.powerSettings.powerBands.powerBinsInHz = [];

% claculate new frequeny bins within band
tempIndecesBinsA = find(powerSettings.powerBands.fftBins>newSettings.bandLimits(1));
tempIndecesBinsB = find(powerSettings.powerBands.fftBins<newSettings.bandLimits(2));
[C,IA,IB]        = intersect(tempIndecesBinsA,tempIndecesBinsB);
binIndecesInBand = tempIndecesBinsA(IA);
binsInBand       = powerSettings.powerBands.fftBins(tempIndecesBinsA(IA));

disp('---')
disp('Calculating equivalent power based on input defined [TD channel (1..4)] and [Power Band]:')
disp(['TD channel = ', num2str(channel)])
disp(['Power Band = ', num2str(newSettings.bandLimits(1)),'Hz-',num2str(newSettings.bandLimits(2)),'Hz']);
disp(['Lower bin = ', num2str(binsInBand(1)), ' Hz']);
disp(['Upper bin = ', num2str(binsInBand(end)), 'Hz']);
disp(['Total bins = ', num2str(length(binsInBand))]);
disp('---')

% add new bins information to new power band in power newSettings
newSettings.powerSettings.powerBands.indices_BandStart_BandStop =...
    [binIndecesInBand(1) binIndecesInBand(end)];

newSettings.powerSettings.powerBands.powerBinsInHz =...
    strcat(sprintf('%0.2f',binsInBand(1)),'Hz-',sprintf('%0.2f',binsInBand(end)),'Hz');

% extract input parameters
tch       = combinedDataTable.DerivedTime;
t         = seconds(tch-tch(1))/1000;

[interval, binStart, binEnd, fftSize, msb_dropped, samp_rate] ...
    = readFFTsettings(newSettings.powerSettings, newSettings.fftSettings);


switch fftSize % actual fft size of device for 64, 250, 1024 fftpoints
    case 64,   fftSizeActual = 62;   n_pad   = 4;
    case 256,  fftSizeActual = 250;  n_pad   = 6;
    case 1024, fftSizeActual = 1000; n_pad   = 24;
end



% Extract neural channel
keych       = combinedDataTable.(['TD_key',num2str(newSettings.tdChannel-1)]);
ampGain     = ampGains.(['Amp',num2str(newSettings.tdChannel)]);

% Transform to RCS units
td_rcs    = transformTDtoRCS(keych, ampGain);


% number of samples since last FFT
samp_elapsed  = (samp_rate*interval/1e3);
% overlap of current window w/ last window
nonoverlap    = (samp_elapsed/fftSizeActual);
overlap       = 1-nonoverlap;

% create Hann window points
hann_win     = hannWindow(fftSize, fftSettings.fftConfig.windowLoad);

% sample 1 of data set where window starts
stime     = 1; 
% caculate an approximate of the total windows over the entire data set
totalTimeWindows = ceil(length(td_rcs)/fftSize/(nonoverlap)); 

i_meas_pb = find(~isnan(combinedDataTable.Power_Band1));


pb_wo_td  = i_meas_pb < fftSize;

if sum(pb_wo_td) > 0
    fprintf(['Not enough TD data based on FFT size '...
             'for power-band indices %d \n'], find(pb_wo_td));

end

% excluding measured power-bands w/o their corresponding TD data
i_meas_pb = i_meas_pb(~pb_wo_td);


for i = 1 : length(i_meas_pb)

    % zero pad samples of interest based on fftSize
    i_td                   = i_meas_pb(i) - (fftSizeActual -1) : i_meas_pb(i);
    td                     = [td_rcs(i_td)', zeros(1, n_pad)];

    % Apply fft on Hann windowed signal
    X                      = fft(td .* hann_win, fftSize);

    % From double to single-sided FFT
    % divide bin centered on zero by two--ignore negative freq
    X                      = [X(1)/2, X(2:fftSize/2)];
    % save FFT
    newPower.FFTinmV{stime+fftSize-1} = fft_rcs_to_mv(X, ampGain);

    % get power and normalize by fftSize
    norm_fft    = 64* (abs(X) / (fftSize)).^2;

    % perform bit shift
    fftPower    = ceil(norm_fft./ 2^(8-msb_dropped));

    % New values into output power variable (vector indexing is important for time alignment)
    newPower.calculatedPower(i_meas_pb(i)) = sum(fftPower(binStart:binEnd));
    
end
end

%% local functions used
function [interval, binStart, binEnd, fftSize, msb_dropped, samp_rate] ...
    ...
    = readFFTsettings(powerSettings, fftSettings)

    interval = powerSettings.fftConfig.interval; % is given in ms
    binStart = powerSettings.powerBands.indices_BandStart_BandStop(1,1);
    binEnd   = powerSettings.powerBands.indices_BandStart_BandStop(1,2);
    fftSize  = powerSettings.fftConfig.size;    % is given in sample points

    % most signifigant bits dropped (number of bits shifted right)
    msb_dropped = str2double(fftSettings.fftConfig.bandFormationConfig(6));

    samp_rate   = fftSettings.TDsampleRates;
end

% transform to rcs units (equation from manufacturer - hardware specific - same in all RC+S devices)
function td_rcs = transformTDtoRCS(keych, ampGain)
    FP_READ_UNITS_VALUE   = 48644.8683623726;    % constant
    lfp_mv                = nan(1,length(keych))';

    lfp_mv(~isnan(keych)) = keych(~isnan(keych))-mean(keych(~isnan(keych))); % remove mean
    lfp_gain              = 250*(ampGain/255);  % actual amplifier gain ch
    td_rcs                = lfp_mv * (lfp_gain*FP_READ_UNITS_VALUE) / (1000*1.2); 
end

function fft_mv = fft_rcs_to_mv(fft_rcs, ampGain)

    FP_READ_UNITS_VALUE   = 48644.8683623726;    % constant
    lfpGain_ch            = 250*(ampGain/255);  % actual amplifier gain ch

    fft_mv                = (1000*1.2)*fft_rcs ./(lfpGain_ch*FP_READ_UNITS_VALUE); 

end
