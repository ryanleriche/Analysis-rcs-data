%% Function Name: td_fft_pd_ld_state_stim_simulation()
%{
Description: Converts .json files from RC+S recordings into .csv files
and MATLAB tables that can be used with the rcssim module. 

The Analysis-rcs-data repo must be included in the path to run 
(https://github.com/openmind-consortium/Analysis-rcs-data/tree/master)

Inputs:

Outputs:


Author: 
    Ryan Leriche, leriche.ryan@gmail.com
%}

%function [data, settings] = extract_data(data_dir)


% sanity check w/ benchtop INS data w/ 20 Hz sine wave
%{
data_dir = '/Users/Leriche/Github/Analysis-rcs-data/testDataSets/Benchtop/ToTestPowerCalc/DeviceNPC700378H';
% data using Analysis-rcs-data

[unifiedDerivedTimes, timeDomainData, ~, ~, ~, ~, ~, PowerData, ~, ~, ~,...
    ~, ~, AdaptiveData, ~, ~, timeDomainSettings, powerSettings,...
    fftSettings, eventLogTable, metaData, stimSettingsOut, stimMetaData,...
    stimLogSettings, DetectorSettings, AdaptiveStimSettings, ...
    AdaptiveEmbeddedRuns_StimSettings, ~] ...
    ...
    = ProcessRCS(data_dir, 2);

dataStreams = {timeDomainData, PowerData, AdaptiveData};
[combinedDataTable] = createCombinedTable(dataStreams, ...
                                          unifiedDerivedTimes, metaData);
%%
pwr_bins = [powerSettings.powerBands.lowerBound, powerSettings.powerBands.upperBound];

%%
calc_power = table();            h = 0;
figure('Units', 'Inches', 'Position', [0, 0, 10, 7])
for i = 1
    for j = 0
        h = h + 1;
        [power, ~]  = calculateNewPower_RBL(combinedDataTable, fftSettings,...
                            powerSettings, metaData, i, pwr_bins(h, :));

        pb_oi        = ['Power_Band', num2str(h)];

        if h == 1
            calc_power.localTime    = power.localTime;
            calc_power.DerivedTimes = power.DerivedTimes; 
        end

        calc_power.(pb_oi) = power.calculatedPower;

        i_pb_val  = ~isnan(calc_power.(pb_oi));

        q = plot(calc_power.localTime(i_pb_val),...
            calc_power.(pb_oi)(i_pb_val),...
            'LineWidth',4);    q.Color(4)=0.5;

        hold on
        i_pb_val  = ~isnan(combinedDataTable.(pb_oi));
        q = plot(combinedDataTable.localTime(i_pb_val), ...
            combinedDataTable.(pb_oi)(i_pb_val),...
            'LineWidth',1.5);

        grid on
        title(['BenchTop INS data from "Analysis-rcs-data" repo', newline,...
            'ch', num2str(i),newline,...
            powerSettings.powerBands.powerBandsInHz{h}])

        if h == 1
            legend({'Simulated'; 'Measured'}); 
        end
    end
end
%}


%% troubleshooting RCS02R as script
data_dir = db_RCS02R.path{i_sess};

% data using Analysis-rcs-data

[unifiedDerivedTimes,...
    timeDomainData, ~, ~,...
    ~, ~, ~, ...
    PowerData, ~, ~,...
    FFTData, ~, ~,...
    ...
    AdaptiveData, ~, ~, timeDomainSettings, powerSettings,...
    fftSettings, eventLogTable, metaData, stimSettingsOut, stimMetaData,...
    stimLogSettings, DetectorSettings, AdaptiveStimSettings, ...
    AdaptiveEmbeddedRuns_StimSettings, ~] ...
    ...
    = ProcessRCS(data_dir, 2);

dataStreams         = {timeDomainData, PowerData, AdaptiveData, FFTData};
[comb_dt]  =  createCombinedTable(dataStreams, ...
                                          unifiedDerivedTimes, metaData);


plt_meta_data ...
    = sprintf('\n%s | %s | duration %s\nTD  %.3f%% of TD samp dropped | %s', ...
    db_RCS02R.sess_name{i_sess}, db_RCS02R.timeStart(i_sess), db_RCS02R.duration(i_sess),...
    per_TD_lost, fftSettings.fftConfig.bandFormationConfig...
    ); 

% calc power-band (PB) time series based on time domain data


tbl_vars   = comb_dt.Properties.VariableNames;
i_pb       = find(contains(tbl_vars, 'Power_Band'));

empty_i_pb = all(isnan(comb_dt(:, i_pb).Variables));

% continue w/ power bands of interest
i_pb       = i_pb(~empty_i_pb);
pb_nums    = cellfun(@(x) str2double(x(end)), tbl_vars(i_pb));


pwr_bins   = [powerSettings.powerBands.lowerBound(pb_nums),...
              powerSettings.powerBands.upperBound(pb_nums)];

sim_tbl    = comb_dt(:, {'localTime', 'DerivedTime',...
                               'TD_key0','TD_key1','TD_key2','TD_key3',...
                               'Power_Band1','Power_Band2','Power_Band3',...
                               'Power_Band4','Power_Band5','Power_Band6',...
                               'Power_Band7','Power_Band8'});            
n_pb       = length(pb_nums);

% generate PB series based on streaming settings
for i = 1 : n_pb

    switch pb_nums(i)
        case 1, TD_ch = 1;      case 2, TD_ch = 1;
        case 3, TD_ch = 2;      case 4, TD_ch = 2;
        case 5, TD_ch = 3;      case 6, TD_ch = 3;
        case 7, TD_ch = 4;      case 8, TD_ch = 4;
    end

    [power, ~]  = calculateNewPower_RBL(comb_dt, fftSettings,...
                        powerSettings, metaData, TD_ch, pwr_bins(i, :));
    if i == 1
        sim_tbl.localTime    = power.localTime;
        sim_tbl.DerivedTimes = power.DerivedTimes;  
    end

    pb_oi              = ['Power_Band', num2str(pb_nums(i))];

    sim_tbl.(['sim_', pb_oi])    = power.calculatedPower;

    sim_tbl.(['fftinmV_ch', num2str(TD_ch)]) = power.FFTinmV;

end


sim_tbl = sim_tbl(~isnan(sim_tbl.Power_Band1),:);



% time-series plot of measured and simulated power-bands:
figure('Units', 'Inches', 'Position', [1, 1, n_pb * 1.5, n_pb*5])

sgtitle(['Power-bands across channels', plt_meta_data], 'Fontsize', 16);


for i = 1 : n_pb
    switch pb_nums(i)
        case 1, TD_ch = 1;      case 2, TD_ch = 1;
        case 3, TD_ch = 2;      case 4, TD_ch = 2;
        case 5, TD_ch = 3;      case 6, TD_ch = 3;
        case 7, TD_ch = 4;      case 8, TD_ch = 4;
    end

    pb_oi     = ['Power_Band', num2str(pb_nums(i))];
    % overlay simulated and measured PB time series
    subplot(n_pb/2, 2, i)
    
    i_pb_val  = find(~isnan(sim_tbl.(['sim_',pb_oi])));
    j_pb_val  = find(~isnan(comb_dt.(pb_oi)));


    
    q = plot(sim_tbl.localTime(i_pb_val),sim_tbl.(['sim_',pb_oi])(i_pb_val),...
            'LineWidth',4);     q.Color(4)=0.6;     hold on   
    

    plot(comb_dt.localTime(j_pb_val), comb_dt.(pb_oi)(j_pb_val),...
            'LineWidth',1.5);   grid on
     
    title(['ch', num2str(TD_ch-1), ' | ',...
        powerSettings.powerBands.powerBandsInHz{pb_nums(i)}], 'FontSize', 14);
    
    if i == 1
        legend({'Simulated'; 'Measured'}); 
    end

    % x and y labe l only on specific subplots
    if mod(i,2) ~=0
        ylabel('RCS Units','FontSize', 14);
    end

    if i == n_pb || i == n_pb - 1
        xlabel('Time','FontSize', 14); 
    end

    Med =  median([sim_tbl.(['sim_',pb_oi])(i_pb_val); comb_dt.(pb_oi)(j_pb_val)]);

    if Med == 0
        Med = 1;
    end

    ylim([Med/30, Med*100])
end

% scatter plot and quanfiying measured and simulated power-band differences
figure('Units', 'Inches', 'Position', [1, 1, n_pb * 2, n_pb*2])

sgtitle(['Power-bands across channels', plt_meta_data], 'Fontsize', 12);

for i = 1 : n_pb
    switch pb_nums(i)
        case 1, TD_ch = 1;      case 2, TD_ch = 1;
        case 3, TD_ch = 2;      case 4, TD_ch = 2;
        case 5, TD_ch = 3;      case 6, TD_ch = 3;
        case 7, TD_ch = 4;      case 8, TD_ch = 4;
    end

    pb_oi     = ['Power_Band', num2str(pb_nums(i))];
    % overlay simulated and measured PB time series
    subplot(n_pb/2, 2, i)

    stim_pwr  = sim_tbl.(['sim_',pb_oi]);
    mea_pwr   = sim_tbl.(pb_oi);


    scatter(stim_pwr, mea_pwr);  hold on; 

    plot(mea_pwr(~isnan(mea_pwr)), mea_pwr(~isnan(mea_pwr)),...
        'k', 'LineWidth', 1.5);


    [r_val, p_val] = corr(stim_pwr, mea_pwr, 'Rows', 'pairwise');

    rmse           = mean(abs(stim_pwr-mea_pwr), 'omitnan');

    per_rmse       = 100 * mean(abs(stim_pwr-mea_pwr)./ abs(mea_pwr), 'omitnan');

    med_pwr        = median([stim_pwr; mea_pwr],'omitnan');

    if med_pwr == 0
        med_pwr = 1;
    end

    text(5*med_pwr,2* med_pwr,...
        [sprintf('r = %.2f',r_val), newline...
         sprintf('pval = %.3g',p_val), newline,...
         sprintf('RMSE = %.3f',rmse),newline,...
         sprintf('percent RMSE = %.3f ',per_rmse)]); 

    % x and y labe l only on specific subplots
    if mod(i,2) ~=0
        ylabel('Measured (RCS Units)','FontSize', 14);
    end

    if i == n_pb || i == n_pb - 1
        xlabel('Simulated (RCS Units)'); 
    end
    xlim([0, 10*med_pwr]);        ylim([0,  10*med_pwr])

    title(['ch', num2str(TD_ch-1),' | ',...
        powerSettings.powerBands.powerBandsInHz{pb_nums(i)}]);
    
    set(gca, 'Fontsize', 12)
end

disp(['***> Only Power Bands ' num2str(pb_nums),' contain real power values']);
disp(plt_meta_data);
%% visualize simulated FFT
fftConfig       = fftSettings.fftConfig;
fftSize         = fftConfig.size;
samp_rate       = timeDomainSettings.samplingRate;

msb_shift       = num2str(fftSettings.fftConfig.bandFormationConfig(6));


i_fft           = find(cellfun(@(x) ~isempty(x), sim_tbl.fftinmV_ch2));
X               = vertcat(sim_tbl.fftinmV_ch2{i_fft})';

% normalized by signal length (not sure why doubling is required)
norm_fft        = abs(X)./(fftSize)/4;

i_real          = find(any(~isnan(norm_fft)));
sim_fftoutput   = norm_fft(:,i_real);


meas_fftoutput  = [FFTData.FftOutput{:}];
med_fft         = median([sim_fftoutput(:); meas_fftoutput(:)],...
                    'all','omitnan');


cb_lim          = [med_fft/20, med_fft*20];
time_range      = [min(FFTData.localTime), max(FFTData.localTime)];
freq_range      = [1, 100];

cfg             = [];
cfg.fftTime     = sim_tbl.localTime(i_fft(i_real));
cfg.fftBins     = powerSettings.powerBands.fftBins;
cfg.fftoutput   = sim_fftoutput;

cfg.cb_lim      = cb_lim;
cfg.x_lim       = time_range;
cfg.y_lim       = freq_range;

cfg.log10       = false;
cfg.title       = 'Simulated FFT';

figure('Units', 'Inches', 'Position', [1, 1, 15, 15])
sgtitle(sprintf('Ch %d \n %s \n %s \n duration: %s',...
                FFTData.Channel(1),...
                db_RCS02R.sess_name{i_sess}, ...
                db_RCS02R.timeStart(i_sess),...
                db_RCS02R.duration(i_sess)));
subplot(211)

sim_spectro  = plot_spectro(cfg);
% visualize Onboard FFT
cfg.fftTime     = FFTData.localTime;
cfg.fftBins     = powerSettings.powerBands.fftBins(1:fftSettings.fftConfig.streamSizeBins);
cfg.fftoutput   = meas_fftoutput;

cfg.title       = 'Onboard FFT';

subplot(212);       
mear_spectro = plot_spectro(cfg);

%% generate LD based on PB series and detector settings
%{
determine LD value, LD state, and stimulation parameters


11/30/22 RBL
completed:
    bringing in LD weights, multiplication, and subtraction vectos
    pulling PB of interest programmically
    ff val scaling automically done

    some of the PB to LD conversion + plotting code

next steps:
    cont. reading RDK manual and slides to see how PB to LD conversion is
    done
        especially the update rate parameter and hold off time
        -> determine everything based off of the dector settings AND PB
        themselves
            -> case for when using novel LD detector settings

        pay attn to time alignment when calling the previous N PBs
            -> could be low-level alginment bug on your part

___________________________________________________________________________
* validate LD sim (based on measured PB) to LD meas

    * time-series 
    * scatter plot w/ RMSE

* validate LD state and stim simulation 

incrp blanking
    - What steps does blanking occur at (FFT -> PB -> LD)?

For detection inputs, what do multiple PB inputs look like?

%}

tbl_vars   = sim_tbl.Properties.VariableNames;
i_pb       = find(contains(tbl_vars, 'Power_Band'));

% continue w/ power bands of interest
pb_nums    = cellfun(@(x) str2double(x(end)), tbl_vars(i_pb));

disp(['For LD, power-bands [' num2str(pb_nums), '] are available.']);

% initalize time-series of LDs as well as settings struct
ld_tbl  = sim_tbl;
LD      = DetectorSettings.Ld0;

% see "Medtronic.NeuroStim.Olympus.DataTypes.Sensing.DetectionEnables"
switch bin2dec(LD.detectionEnable_BinaryCode)
    % No enables are set. This means that single threshold detect mode is used  
    case 0,     detect = 'None';

    % Enable dual threshold (band) detect. If this bit is not set, single detect mode is used.  
    case 1,     detect = 'DualThresholdEnabled'; 
   
    % Blank both LDs based on a state change from this LD. The blank 
    % duration specified for each LD is used (if 0 duration is set, 
    % then no blanking will occur). If this bit is not set, only this LD is blanked.
    case 2,     detect = 'BlankBoth';
end


n_feat      = length(LD.detectionInputs_BinaryCode);
features    = cell(n_feat,1);

% from LD settings, **explicilty** pull out enabled PBs
for i = 1 : n_feat

    % see "Medtronic.NeuroStim.Olympus.DataTypes.Sensing.DetectionInputs"
    switch bin2dec(LD.detectionInputs_BinaryCode)
        case 0,     pb_detect = 'None';
            
        % Power channel 0 band 0 and 1 
        case 1,     pb_detect = 'Power_Band1';    case 2,   pb_detect = 'Power_Band2';
    
        % Power channel 1 band 0 and 1  
        case 4,     pb_detect = 'Power_Band3';    case 8,    pb_detect = 'Power_Band4';
    
        % Power channel 2 band 0 and 1 
        case 16,    pb_detect = 'Power_Band5';    case 32,   pb_detect = 'Power_Band6';
    
        % Power channel 3 band 0 and 1  
        case 64,    pb_detect = 'Power_Band7';    case 128,  pb_detect = 'Power_Band8';
    end
    features{i} = pb_detect;
end

ff_val  = LD.fractionalFixedPointValue;

% weight, multiply, subtract, and BiasTerm vectors in LD units
W     = [LD.features.weightVector] ./ 2^ff_val;
M_vec     = [LD.features.normalizationMultiplyVector] ./ 2^ff_val;

S     = [LD.features.normalizationSubtractVector]; % subtract vector lacks 
                                                 % fractional fixed point 
                                                 % value scaling
thres = LD.biasTerm ./ 2^ff_val;

%temp_tbl = ld_tbl(~isnan(ld_tbl.Power_FftSize),:);


ld_features  = unique(features);
n_feat_in_LD = length(ld_features);


%%
state_tbl    = table();

for i = 1 : n_feat_in_LD
    i_rec     = find(~isnan(ld_tbl.(ld_features{i})));
    feat_pwr  = ld_tbl.(features{i})(i_rec);

    state_tbl.(ld_features{i}) = feat_pwr;

end

state_tbl.ld_pwr = nan(height(state_tbl),1);

% using eariliest streamed LD state to inform starting
% --> could be better to use INS to check for missed streamed data

i_state            = find(cellfun(@(x) any(~isnan(x)),...
                        comb_dt.Adaptive_CurrentAdaptiveState));

state_tbl.desc(1)  = comb_dt.Adaptive_CurrentAdaptiveState(i_state(1));

pb_count      = 1;      abv_count     = 0;      blw_count     = 0;


fprintf('LD updateRate = %d (i.e., mean of previous %d power-bands, and current power-band). \n',...
    LD.updateRate, LD.updateRate -1);

for i = 2 : height(state_tbl)

    if pb_count == LD.updateRate
        pb_count = 1;
        for j = 1 : n_feat_in_LD

            avg_pwr = mean(state_tbl.(ld_features{j})...
                                   (i + 1 - LD.updateRate : i));


            temp_feat(j)  = W(j) * M_vec(j) * (avg_pwr - S(j)); %#ok<SAGROW> 

        end 

        state_tbl.ld_pwr(i) =  sum(temp_feat, 2);

        % evaluate state counters based on LD power update
        if state_tbl.ld_pwr(i) > thres(1)

            abv_count = abv_count + 1;
            blw_count = 0;

        elseif state_tbl.ld_pwr(i) < thres(1)
    
            blw_count = blw_count + 1;
            abv_count = 0;
    
        end

        % determine state from state counters
        if abv_count >= LD.onsetDuration
    
            state_tbl.desc{i} = 'State 1';
    
        elseif blw_count >= LD.terminationDuration
    
            state_tbl.desc{i} = 'State 0';
    
        else
    
            state_tbl.desc{i} = state_tbl.desc{i- LD.updateRate};
        end

        % from states, see if state changed
        if ~strcmp(state_tbl.desc{i}, state_tbl.desc{i- LD.updateRate})
            fprintf('State change at %d \n', i)
        end

    else
        pb_count = pb_count + 1;
    end
end




%% on RCS, new LD is calculated only every LD.updateRate
state                 = table();
state.ld_pwr          = sum_pwr(1 : LD.updateRate  : length(sum_pwr));
    
% holdoffTime is N of updateRate before starting LD 




%%
state.above_thres     = state.ld_pwr  > thres(1);
state.thres_change    = [state.above_thres(1); diff(state.above_thres)];

below_to_above        = find(state.thres_change   == 1);
above_to_below        = find(state.thres_change   == -1);

% edge-case where session ends before LD falls below threshold
% -> have last detection be end
if length(below_to_above) > length(above_to_below)

    above_to_below = [above_to_below; height(state)];

end

state.conseq_above    = zeros(height(state),1);

for i = 1 : length(below_to_above)
    i_above_thres_chunk = below_to_above(i) : above_to_below(i) -1;

    if length(i_above_thres_chunk) > LD.onsetDuration
        state.desc(i_above_thres_chunk) = repmat({'State1'}, length(i_above_thres_chunk),1);
    end
end

%%

%{
thres_change(thres_change==-1) = 0;

span_locs               = cumsum(thres_change);
span_locs(~supra_thres) = 0;
aboveThresholdIndex = find(supra_thres ==1);
notConsecutiveIndex = [true diff(aboveThresholdIndex) ~= 1];

sumIndex   = cumsum(notConsecutiveIndex);
spanLength = histc(sumIndex, 1:sumIndex(end));

goodSpans = find(spanLength >= LD.onsetDuration);
allInSpans = find(ismember(span_locs, goodSpans));

samp_LD_blank   = ceil(LD.blankingDurationUponStateChange * (fftConfig.interval/2000)* samp_rate);
%}

%
% sum across features of scaled_feat matrix to arrive at LD1 output 

ld_tbl.('sim_ld0_PowerOutput')          = nan(height(ld_tbl),1);
ld_tbl.('sim_ld0_PowerOutput')(i_rec(1 : LD.updateRate: length(sum_pwr)))...
    = state.ld_pwr;

% % holdoffTime is N of updateRate before starting LD 
ld_holdoff_samp    = LD.holdoffTime * LD.updateRate * fftConfig.interval/2250 *samp_rate;

ld_tbl.('sim_ld0_PowerOutput')(1:ld_holdoff_samp) = 0;

    

% LD measured vs. LD simulated | time series
figure('Units', 'Inches', 'Position', [1, 1, 16, 8])

sgtitle(sprintf(...
    'LD simulation performance \n  %s  |  %s  |  duration of %s \n', ...
    db_RCS02R.sess_name{i_sess}, ...
    db_RCS02R.timeStart(i_sess),...
    db_RCS02R.duration(i_sess)), 'Fontsize', 14);



% overlay simulated and measured LD series
subplot(211)

%title(ld_features,'Fontsize', 12);

i_ld_val  = ~isnan(ld_tbl.Adaptive_Ld0_output);
j_ld_val  = ~isnan(ld_tbl.sim_ld0_PowerOutput);

q = stairs(ld_tbl.localTime(i_ld_val), ld_tbl.Adaptive_Ld0_output(i_ld_val),...
        'LineWidth',3);     q.Color(4)=0.7;     hold on   

stairs(ld_tbl.localTime(j_ld_val), ld_tbl.sim_ld0_PowerOutput(j_ld_val),...
        'LineWidth',2);   grid on


if ~contains(detect, 'DualThresholdEnabled')
    yline(LD.biasTerm(1), 'LineWidth', 3)

    legend({'Measured', 'Simulated', 'LD threshold'}, ...
        'Fontsize', 12, 'Location','bestoutside');
end

box off; legend box off; grid on; grid minor;

title(['ch', num2str(TD_ch-1), ' | ',...
    powerSettings.powerBands.powerBandsInHz{pb_nums(i)}], 'FontSize', 14);


% x and y labe l only on specific subplots
if mod(i,2) ~=0
    ylabel('RCS Units','FontSize', 14);
end

if i == n_pb || i == n_pb - 1
    xlabel('Time','FontSize', 14); 
end

Med =  median([sim_tbl.(pb_oi)(i_pb_val); comb_dt.(pb_oi)(j_pb_val)]);

if Med == 0
    Med = 1;
end

ylim([0, Med*50])


%%










%% Create a simplified data table

% timestamp
timestamp = comb_dt.DerivedTime/1000;

% Time-Domain data
td1 = comb_dt.TD_key0;
td2 = comb_dt.TD_key1;
td3 = comb_dt.TD_key2;
td4 = comb_dt.TD_key3;

% Power-Band data
ch1_pb1 = comb_dt.Power_Band1;
ch1_pb2 = comb_dt.Power_Band2;

ch2_pb1 = comb_dt.Power_Band3;
ch2_pb2 = comb_dt.Power_Band4;

ch3_pb1 = comb_dt.Power_Band5;
ch3_pb2 = comb_dt.Power_Band6;

ch4_pb1 = comb_dt.Power_Band7;
ch4_pb2 = comb_dt.Power_Band8;



% incrop measured LD outputs into sim_table
if any(...
        contains(comb_dt.Properties.VariableNames, 'Adaptive_Ld0_output'))

    ld1 = correct_ld(...
            comb_dt.Adaptive_Ld0_output);
else
    ld1 = nan(height(comb_dt),1);
end

if any(...
        contains(comb_dt.Properties.VariableNames, 'Adaptive_Ld0_output'))

    ld2 = correct_ld(...
                comb_dt.Adaptive_Ld1_output);
else
    ld2 = nan(height(comb_dt),1);
end

if any(...
                contains(comb_dt.Properties.VariableNames, 'Adaptive_CurrentAdaptiveState'))
    % state
    state = isolate_state_vector(comb_dt);
    
    % stim
    stim = isolate_stim_vector(comb_dt);

else

    state = nan(height(comb_dt),1);
    stim  = nan(height(comb_dt),1);
end

% combine into a single data table and write to file
data = table(timestamp, ...
             td1, td2, td3, td4, ...
             ch1_pb1, ch1_pb2, ch2_pb1, ch2_pb2, ch3_pb1, ch3_pb2, ch4_pb1, ch4_pb2, ...
             ld1, ld2, ...
             state, stim);

writetable(data, fullfile(data_dir, 'td_pb_ld_aligned.csv'))

%% Create a settings file
% If settings were changed in the middle of the selected session, this may
% not work appropriately.
fs_td           = timeDomainSettings.TDsettings{1}(1).sampleRate;
fftSize        = fftSettings.fftConfig.size;
interval        = fftSettings.fftConfig.interval;
bit_shift       = str2double(fftSettings.fftConfig.bandFormationConfig(6));

percent_han_win = str2double(fftSettings.fftConfig.windowLoad(1:end-6));

band_edges_hz = '[';
for i = 1:8
    band = strcat('[',num2str(powerSettings.powerBands(1).lowerBound(i)), ...
                  ',', ...
                  num2str(powerSettings.powerBands(1).upperBound(i)), ...
                  '],');
    band_edges_hz = strcat(band_edges_hz, band);
end
band_edges_hz(end) = ']';      band_edges_hz      = {band_edges_hz};

%{
pwr = powerSettings.powerBands;
band_bins_hz = '[';
for i = 1 : 8

    band   =  strcat('[',num2str(pwr.fftBins(pwr.indices_BandStart_BandStop(i,1))),...
              ',',...
              num2str(pwr.fftBins(pwr.indices_BandStart_BandStop(i,2))), '],');

    band_bins_hz = [band_bins_hz, band]; %#ok<AGROW> 

end
band_bins_hz(end) = ']';         band_bins_hz       = {band_bins_hz};
%}

amp_gain           = {metaData.ampGains.Amp1, metaData.ampGains.Amp2, metaData.ampGains.Amp3, metaData.ampGains.Amp4};


% LD settings
ld1_subtract        = {DetectorSettings.Ld0.features.normalizationSubtractVector};
ld2_subtract        = {DetectorSettings.Ld1.features.normalizationSubtractVector};

ld1_multiply        = {DetectorSettings.Ld0.features.normalizationMultiplyVector};
ld2_multiply        = {DetectorSettings.Ld1.features.normalizationMultiplyVector};

ld1_weights         = {DetectorSettings.Ld0.features.weightVector};
ld2_weights         = {DetectorSettings.Ld1.features.weightVector};

ld1_update_rate     = {DetectorSettings.Ld0.updateRate};
ld2_update_rate     = {DetectorSettings.Ld1.updateRate};

ld1_frac_fix_vals   = {DetectorSettings.Ld0.fractionalFixedPointValue};
ld2_frac_fix_vals   = {DetectorSettings.Ld1.fractionalFixedPointValue};

ld1_blank_dur       = {DetectorSettings.Ld0.blankingDurationUponStateChange};
ld2_blank_dur       = {DetectorSettings.Ld1.blankingDurationUponStateChange};

ld1_onset           = {DetectorSettings.Ld0.onsetDuration};
ld2_onset           = {DetectorSettings.Ld1.onsetDuration};

ld1_termination     = {DetectorSettings.Ld0.terminationDuration};
ld2_termination     = {DetectorSettings.Ld1.terminationDuration};

ld1_dual_threshold  = {DetectorSettings.Ld0.biasTerm};
ld2_dual_threshold  = {DetectorSettings.Ld1.biasTerm};

ld1_threshold       = {DetectorSettings.Ld0.biasTerm(1)};
ld2_threshold       = {DetectorSettings.Ld1.biasTerm(1)};

% now mA adminstered at given states
state1_mA     = AdaptiveEmbeddedRuns_StimSettings.states.state0_AmpInMilliamps;
state2_mA     = AdaptiveEmbeddedRuns_StimSettings.states.state1_AmpInMilliamps;
state3_mA     = AdaptiveEmbeddedRuns_StimSettings.states.state2_AmpInMilliamps;

state4_mA     = AdaptiveEmbeddedRuns_StimSettings.states.state3_AmpInMilliamps;
state5_mA     = AdaptiveEmbeddedRuns_StimSettings.states.state4_AmpInMilliamps;
state6_mA     = AdaptiveEmbeddedRuns_StimSettings.states.state5_AmpInMilliamps;

state7_mA     = AdaptiveEmbeddedRuns_StimSettings.states.state6_AmpInMilliamps;
state8_mA     = AdaptiveEmbeddedRuns_StimSettings.states.state7_AmpInMilliamps;
state9_mA     = AdaptiveEmbeddedRuns_StimSettings.states.state8_AmpInMilliamps;

rise_mA_sec    = {AdaptiveEmbeddedRuns_StimSettings.deltas{1}.rise};
fall_mA_sec    = {AdaptiveEmbeddedRuns_StimSettings.deltas{1}.fall};

pwr_bins       = {pwr_bins};
% combine into a single data table and write to file
settings = table(fs_td,...
                 fftSize, ...
                 interval, ...
                 bit_shift, ...
                 band_edges_hz,...
                 percent_han_win,...
                 amp_gain,...
                 ... % LD1 settings
                 ld1_subtract, ld1_multiply, ld1_weights,...
                 ld1_update_rate,  ld1_onset, ld1_termination,...
                 ld1_blank_dur,...
                 ld1_dual_threshold, ld1_threshold,...
                 ld1_frac_fix_vals,...
                 ...
                 ld2_subtract, ld2_multiply, ld2_weights,...
                 ld2_update_rate,  ld2_onset, ld2_termination,...
                 ld2_blank_dur,...
                 ld2_dual_threshold, ld2_threshold,...
                 ld2_frac_fix_vals,...
                 ...
                 state1_mA, state2_mA, state3_mA,...
                 state4_mA, state5_mA, state6_mA,...
                 state7_mA, state8_mA, state9_mA,...
                 rise_mA_sec, fall_mA_sec );

writetable(settings, fullfile(data_dir, 'stream_settings.csv'))


%
%{
# settings['blank_both']      = settings['blank_both'].apply(literal_eval)
%}


%%
%end
