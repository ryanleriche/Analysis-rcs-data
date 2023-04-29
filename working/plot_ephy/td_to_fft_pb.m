function sim_tbl     = td_to_fft_pb(i_sess, db_RCS02R)
%{
Description: Converts .json files from RC+S recordings into .csv files
and MATLAB tables that can be used with the rcssim module. 


Ryan Leriche, ryan.leriche@ucsf.edu (leriche.ryan@gmail.com) December 2022
%}

% can troubleshoot as script by commenting out 'function' and 'end' and uncommenting the line below
%i_sess = 16;

% pull in RCS data via ProcessRCS
data_dir = db_RCS02R.path{i_sess};

[unifiedDerivedTimes,...
    timeDomainData, ~, ~,...
    ~, ~, ~, ...
    PowerData, ~, ~,...
    FFTData, ~, ~,...
    ...
    AdaptiveData, ~, ~, ~, powerSettings,...
    fftSettings, ~, metaData, ~, ~,...
    ~, ~, ~, ...
    ~, ~] ...
    ...
    = ProcessRCS(data_dir, 3);

dataStreams   = {timeDomainData, PowerData, AdaptiveData, FFTData};
comb_dt       =  createCombinedTable(dataStreams, unifiedDerivedTimes, metaData);


plt_meta_data ...
    = sprintf(['\n%s | %s | duration %s\nTD  %.3f%% of TD samp dropped | %s \n', ...
    'stim @ %s %s'],...
    db_RCS02R.sess_name{i_sess}, db_RCS02R.timeStart(i_sess), db_RCS02R.duration(i_sess),...
    db_RCS02R.per_TD_lost(i_sess), db_RCS02R.fftSettings{i_sess}.fftConfig.bandFormationConfig,...
    db_RCS02R.stimReg{i_sess}, db_RCS02R.stimLogSettings{i_sess}.stimParams_prog1{1}); 

disp(plt_meta_data);
% calc power-band (PB) time series based on time domain data


tbl_vars   = comb_dt.Properties.VariableNames;
i_pb       = find(contains(tbl_vars, 'Power_Band'));

empty_i_pb = all(isnan(comb_dt(:, i_pb).Variables));

% continue w/ power bands of interest
i_pb       = i_pb(~empty_i_pb);
pb_nums    = cellfun(@(x) str2double(x(end)), tbl_vars(i_pb));


pwr_bins   = [powerSettings.powerBands.lowerBound(pb_nums),...
              powerSettings.powerBands.upperBound(pb_nums)];

sim_tbl    = comb_dt(:,   {'localTime', 'DerivedTime',...
                           'TD_key0','TD_key1','TD_key2','TD_key3',...
                           'Power_Band1','Power_Band2','Power_Band3',...
                           'Power_Band4','Power_Band5','Power_Band6',...
                           'Power_Band7','Power_Band8'});            
n_pb       = length(pb_nums);

%% generate PB series based on streaming settings
for i = 1 : n_pb

    switch pb_nums(i)
        case 1, TD_ch = 0;      case 2, TD_ch = 0;
        case 3, TD_ch = 1;      case 4, TD_ch = 1;
        case 5, TD_ch = 2;      case 6, TD_ch = 2;
        case 7, TD_ch = 3;      case 8, TD_ch = 3;
    end

    %****where TD to FFT to PB calculation occurs
    [tmp_sim_tbl, ~]  = calculateNewPower_RBL(comb_dt, fftSettings,...
                        powerSettings, metaData, TD_ch, pwr_bins(i, :));

    pb_oi              = ['Power_Band', num2str(pb_nums(i))];

    sim_tbl.(['sim_', pb_oi])    = tmp_sim_tbl.calculatedPower;

    sim_tbl.(['fftinmV_ch', num2str(TD_ch)]) = tmp_sim_tbl.FFTinmV;

end
%% **** rest of fxn is for plotting (saved automatically as .pngs in path) ****
% time-series plot of measured and simulated power-bands:

sim_tbl = sim_tbl(~isnan(sim_tbl.Power_Band1),:);
if n_pb > 1
    figure('Units', 'Inches', 'Position', [1, 1, n_pb * 2.25 , n_pb*1.75])
else
    figure('Units', 'Inches', 'Position', [1, 1, 20, 20])
end

sgtitle(['Power-bands across channels', plt_meta_data], 'Fontsize', 12, 'Interpreter','none');

for i = 1 : n_pb
    switch pb_nums(i)
        case 1, TD_ch = 0;      case 2, TD_ch = 0;
        case 3, TD_ch = 1;      case 4, TD_ch = 1;
        case 5, TD_ch = 2;      case 6, TD_ch = 2;
        case 7, TD_ch = 3;      case 8, TD_ch = 3;
    end

    pb_oi     = ['Power_Band', num2str(pb_nums(i))];
    % overlay simulated and measured PB time series
    if n_pb > 1
        subplot(n_pb/2, 2, i)
    end
    
    i_pb_val  = find(~isnan(sim_tbl.(['sim_',pb_oi])));
    j_pb_val  = find(~isnan(comb_dt.(pb_oi)));

    q = plot(sim_tbl.localTime(i_pb_val), sim_tbl.(['sim_',pb_oi])(i_pb_val),...
            'LineWidth',4);     q.Color(4)=0.6;     hold on   
    

    plot(comb_dt.localTime(j_pb_val), comb_dt.(pb_oi)(j_pb_val),...
            'LineWidth',1.5);   grid on
     
    title_txt = ['ch', num2str(TD_ch),...
        ' (' db_RCS02R.timeDomainSettings{i_sess}.(['chan', num2str(TD_ch +1)]){1}, ')', newline,...
        pb_oi, ' (', powerSettings.powerBands.powerBandsInHz{pb_nums(i)}, ')'];

    title(title_txt, 'FontSize', 10, 'Interpreter','none');
    
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

    mean_pwr =  mean([sim_tbl.(['sim_',pb_oi])(i_pb_val); comb_dt.(pb_oi)(j_pb_val)]);
    std_pwr  = std([sim_tbl.(['sim_',pb_oi])(i_pb_val); comb_dt.(pb_oi)(j_pb_val)]);
    if mean_pwr == 0
        mean_pwr = 1;
    end

    ylim([0, mean_pwr+3*std_pwr])

    set(gca, 'Fontsize', 12);       grid on; grid minor
end

save_dir = [cd,'/plot_ephy/session_by_session/'...
        'SummitContinuousBilateralStreaming/RCS02R/', db_RCS02R.sess_name{i_sess},'/'];

if ~isfolder(save_dir)
    mkdir(save_dir)
end

saveas(gcf, [save_dir,'pb_series_simulated_vs_measured.png'])

% scatter plot
if n_pb > 1
    figure('Units', 'Inches', 'Position', [1, 1, n_pb * 2.25 , n_pb*1.75])
else
    figure('Units', 'Inches', 'Position', [1, 1, 20, 20])
end

sgtitle(['Power-bands across channels', plt_meta_data], 'Fontsize', 12, 'Interpreter','none');

for i = 1 : n_pb
    switch pb_nums(i)
        case 1, TD_ch = 0;      case 2, TD_ch = 0;
        case 3, TD_ch = 1;      case 4, TD_ch = 1;
        case 5, TD_ch = 2;      case 6, TD_ch = 2;
        case 7, TD_ch = 3;      case 8, TD_ch = 3;
    end

    pb_oi     = ['Power_Band', num2str(pb_nums(i))];

    if n_pb > 1
        subplot(n_pb/2, 2, i)
    end
    stim_pwr  = sim_tbl.(['sim_',pb_oi]);
    mea_pwr   = sim_tbl.(pb_oi);


    scatter(stim_pwr, mea_pwr);  hold on; 

    plot(mea_pwr(~isnan(mea_pwr)), mea_pwr(~isnan(mea_pwr)),...
        'k', 'LineWidth', 1.5);

    [r_val, p_val] = corr(stim_pwr, mea_pwr, 'Rows', 'pairwise');

    rmse           = mean(abs(stim_pwr-mea_pwr), 'omitnan');

    mean_pwr        = mean([stim_pwr; mea_pwr],'omitnan');
    std_pwr         = std([stim_pwr; mea_pwr],'omitnan');
    if mean_pwr == 0
        mean_pwr = 1;
    end

    % x and y labe l only on specific subplots
    if mod(i,2) ~=0
        ylabel('Measured (RCS Units)','FontSize', 14);
    end

    if i == n_pb || i == n_pb - 1
        xlabel('Simulated (RCS Units)'); 
    end

    xlim([0, 3*std_pwr+mean_pwr]);  
    ylim([0, 3*std_pwr+mean_pwr]);

    title_txt = ['ch', num2str(TD_ch),...
        ' (' db_RCS02R.timeDomainSettings{i_sess}.(['chan', num2str(TD_ch +1)]){1}, ')', newline,...
        pb_oi, ' (', powerSettings.powerBands.powerBandsInHz{pb_nums(i)}, ')'];

    title(title_txt, 'Interpreter','none');
    
    text(mean_pwr, mean_pwr + 1.5*std_pwr,...
        [sprintf('r^2 = %.3f',r_val^2), newline...
         sprintf('pval = %.3g',p_val), newline,...
         sprintf('RMSE = %.3f',rmse)]);

    set(gca, 'Fontsize', 12);       grid on; grid minor
    
end
saveas(gcf, [save_dir,'pb_scatter_simulated_vs_measured.png'])
end