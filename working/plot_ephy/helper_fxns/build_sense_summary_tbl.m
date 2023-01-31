function ...
    [app_ss_tbl, proc_app_log, td_fft_pm_ld_state_stim_tbl]...
    ...
    = build_sense_summary_tbl(...
    ...
    db_RCSXXX, proc_app_log)
%%
%{
take every sub-setting (i.e., time-domain, FFT, power-band, linear discrinant,
state deicions), from every streaming session into a single-summary table

THEN

find nearest streaming session--per INS log entry--to explicilty have ALL
sensing/stimulation parameters

RETURN

AppLog Streaming Session aligned table (app_ss_tbl) for
subsequent plotting/analysis of aDBS settings

%}


%%% uncomment these to trouble-shoot as script:


% db_RCSXXX    = db.RCS02R;
% proc_app_log = INS_logs.RCS02R.app;


%% Expand DetectorSettings (i.e., the LD parameters)
[~, i_u]        = unique(db_RCSXXX.sess_name);
db_RCSXXX       = db_RCSXXX(i_u,:);

i_detSett        = cellfun(@(x) ~isempty(x), db_RCSXXX.DetectorSettings);


% take LAST entry from as generated from Streaming Session
tmp_lds           = cellfun(@(x) x(end,:), db_RCSXXX.DetectorSettings(i_detSett),...
                             'UniformOutput',false);

tmp_LD_tbl  = vertcat(tmp_lds{:});

% explicitly include Streaming Session API time, name, and path
tmp_LD_tbl.sess_name       = db_RCSXXX.sess_name(i_detSett);


%%% expand LD0 and LD1 settings as own variable for easy plotting/analysis later

LD0_tbl         = struct2table(tmp_LD_tbl.Ld0);
tmp_ld_feat_tbl = table();

for i=1:4
    LD0_feat = struct2table(cellfun(@(x) x(i), LD0_tbl.features));

    LD0_feat = renamevars(LD0_feat,...
               {'normalizationMultiplyVector', 'normalizationSubtractVector','weightVector'},...
               {['normMultiply', num2str(i-1)],['normSubtract', num2str(i-1)],['normWeight', num2str(i-1)]});


    tmp_ld_feat_tbl = [tmp_ld_feat_tbl, LD0_feat]; %#ok<AGROW> 

end

tmp_ld_feat_tbl.biasTerm0 = cellfun(@(x) x(1), LD0_tbl.biasTerm);
tmp_ld_feat_tbl.biasTerm1 = cellfun(@(x) x(2), LD0_tbl.biasTerm);

LD0_tbl = [tmp_ld_feat_tbl,  LD0_tbl];
LD0_tbl = removevars(LD0_tbl, {'biasTerm', 'features'});


%%% repeat for LD1
LD1_tbl = struct2table(tmp_LD_tbl.Ld1);

tmp_ld_feat_tbl = table();

for i=1:4
    LD1_feat = struct2table(cellfun(@(x) x(i), LD1_tbl.features));

    LD1_feat = renamevars(LD1_feat,...
               {'normalizationMultiplyVector', 'normalizationSubtractVector','weightVector'},...
               {['normMultiply', num2str(i-1)],['normSubtract', num2str(i-1)],['normWeight', num2str(i-1)]});


    tmp_ld_feat_tbl = [tmp_ld_feat_tbl, LD1_feat]; %#ok<AGROW> 

end

tmp_ld_feat_tbl.biasTerm0 = cellfun(@(x) x(1), LD1_tbl.biasTerm);
tmp_ld_feat_tbl.biasTerm1 = cellfun(@(x) x(2), LD1_tbl.biasTerm);

LD1_tbl = [tmp_ld_feat_tbl,  LD1_tbl];
LD1_tbl = removevars(LD1_tbl, {'biasTerm', 'features'});

% reintroduce LD0/1 names into table variable names themselves
LD0_tbl = renamevars(LD0_tbl, LD0_tbl.Properties.VariableNames, ...
            cellfun(@(x) ['LD0_',x], LD0_tbl.Properties.VariableNames, 'UniformOutput',false));

LD1_tbl = renamevars(LD1_tbl, LD1_tbl.Properties.VariableNames, ...
            cellfun(@(x) ['LD1_',x], LD1_tbl.Properties.VariableNames, 'UniformOutput',false));

% now concatenate expanded LD0 and LD1 table to Streaming Session name/time
LD_State_tbl = [tmp_LD_tbl(:, 'sess_name'),...
                LD0_tbl,...
                LD1_tbl];

%% expand State Settings
% ensure alignment of streaming session DetectorSettings and adaptiveStimSettings

i_adapStimSett      = cellfun(@(x) ~isempty(x), db_RCSXXX.AdaptiveStimSettings);
tmp_state           = cellfun(@(x) x(end,:), db_RCSXXX.AdaptiveStimSettings(i_adapStimSett),...
                             'UniformOutput',false);

tmp_State_tbl       = vertcat(tmp_state{:});

delta_tbl = table();

for i=1:4
    state_deltas = struct2table(cellfun(@(x) x(i), tmp_State_tbl.deltas));

    delta = renamevars(state_deltas,...
               {'fall', 'rise'},...
               {['fallInMilliAmpsPerSecProg', num2str(i-1)],['riseInMilliAmpsPerSec', num2str(i-1)]});


    delta_tbl = [delta_tbl, delta]; %#ok<AGROW> 
end

state_amps     = struct2table(tmp_State_tbl.states);
state_i_Prog_j = table;

for i=1:8

    tmp_state = state_amps(:, sprintf('state%0.0f_AmpInMilliamps', i-1));

    tmp_sess_by_prog_amp = tmp_state.Variables;
   
    for j=1:4
        state_i_Prog_j.(sprintf('state%0.0f_AmpInMilliampsProg%0.0f', i-1,j-1)) ...
            = tmp_sess_by_prog_amp(:,j);
    end
end

state_delta_tbl = [db_RCSXXX(i_adapStimSett, 'sess_name'), state_i_Prog_j, delta_tbl];

%% expand TD settings
i_tdSett   = cellfun(@(x) ~isempty(x), db_RCSXXX.timeDomainSettings);
tmp_td     = cellfun(@(x) x(end,:), db_RCSXXX.timeDomainSettings(i_tdSett),...
                             'UniformOutput',false);

tmp_td_tbl = removevars(vertcat(tmp_td{:}),...
              {'recNum', 'duration', 'timeStart', 'timeStop', 'samplingRate'});

TD_tbl     = db_RCSXXX(i_tdSett,'sess_name');
for i=1:4

    ch_i_td_tbl = struct2table(cellfun(@(x) x(i), tmp_td_tbl.TDsettings));

    ch_i_td_tbl = renamevars(ch_i_td_tbl, ch_i_td_tbl.Properties.VariableNames,...
                     cellfun(@(x) sprintf('Ch%0.0f_%s', i-1, x),...
                            ch_i_td_tbl.Properties.VariableNames,...
                            'UniformOutput', false));

    TD_tbl      =  [TD_tbl, ch_i_td_tbl]; %#ok<AGROW> 

end

td_vars = TD_tbl.Properties.VariableNames;

% remove sampleRate as every TD channel MUST have same sampleRate and this
% is parsimonioulsy shown in the FFT settings
TD_tbl= removevars(TD_tbl, td_vars(contains(td_vars, 'sampleRate')));

%% expand FFT settings
i_fftSett    = cellfun(@(x) ~isempty(x), db_RCSXXX.fftSettings);
tmp_fft      = cellfun(@(x) x(end,:), db_RCSXXX.fftSettings(i_fftSett),...
                             'UniformOutput',false);
tmp_fft_tbl           = vertcat(tmp_fft{:});
tmp_fft_tbl.sess_name = db_RCSXXX.sess_name(i_fftSett);

% yep, FFT settings are quite simple to expand
tmp_fftcfg_tbl = struct2table(tmp_fft_tbl.fftConfig);

fftSett_tbl = [tmp_fft_tbl(:, {'sess_name', 'TDsampleRates'}),...
               tmp_fftcfg_tbl];

fftSett_tbl = renamevars(fftSett_tbl, {'interval', 'size'}, {'fft_intervalInMilliseconds', 'fft_sizeInSamples'});

%% expand power-band settings
i_pwrSett    = cellfun(@(x) ~isempty(x), db_RCSXXX.powerSettings);
tmp_pwr      = cellfun(@(x) x(end,:), db_RCSXXX.powerSettings(i_pwrSett),...
                             'UniformOutput',false);

tmp_pwr_tbl  = removevars(vertcat(tmp_pwr{:}),...
                  {'recNum', 'duration', 'timeStart',...
                  'timeStop', 'fftConfig', 'TDsampleRates'});

tmp_pwr_tbl.sess_name =db_RCSXXX.sess_name(i_pwrSett);


tmp_pwrband_tbl = struct2table(tmp_pwr_tbl.powerBands);
tmp_band_tbl  = table();
i_ch            = [0 0 1 1 2 2 3 3];
for i =1:8

    powerBand_i = cellfun(@(x) x(i), tmp_pwrband_tbl.powerBandsInHz);
    ch_n_pb_lbl = sprintf('Ch%0.0f_powerBandInHz%0.0f', i_ch(i),i-1);

    tmp_band_tbl(:,ch_n_pb_lbl) = powerBand_i;


    powerBin_i  = cellfun(@(x) x(i), tmp_pwrband_tbl.powerBinsInHz);
    ch_n_pb_lbl = sprintf('Ch%0.0f_powerBinInHz%0.0f', i_ch(i),i-1);

    tmp_band_tbl(:,ch_n_pb_lbl) = powerBin_i;

end
% note power_band_tbl goes:
% Sess Name -> power band setting -> corresponding fft settings

power_band_tbl = [tmp_pwr_tbl(:, 'sess_name'), tmp_band_tbl, ...
                  tmp_pwrband_tbl(:, ...
                        {'fftBins', 'indices_BandStart_BandStop'})];
%% expand StimSettings
i_stimSett   = cellfun(@(x) ~isempty(x), db_RCSXXX.stimSettingsOut);

tmp_stim_tbl = db_RCSXXX.stimSettingsOut(i_stimSett);
tmp_stim_tbl = cellfun(@(x) x(end,:), tmp_stim_tbl, 'UniformOutput', false);



stim_tbl     = [db_RCSXXX(i_stimSett, {'sess_name'}),...
                vertcat(tmp_stim_tbl{:})];

i_entry    = find(~cellfun(@isempty, stim_tbl.GroupA));

Groups              = {'GroupA', 'GroupB', 'GroupC', 'GroupD'};
Progs               = {'Prog0', 'Prog1', 'Prog2', 'Prog3'};

% replace empty Groups w/ previous non-empty group
for i=1:height(stim_tbl)

    if isempty(stim_tbl{i, 'GroupA'}{1})

        emp_diff = (i - i_entry);

        [~, i_prev] = min(emp_diff(emp_diff>0));

        stim_tbl(i,Groups) = stim_tbl(i_entry(i_prev), Groups);

    end
end



for i =1:4

    tmp_group = struct2table([stim_tbl{:, Groups{i}}{:}]);

    group_prog_tbl = table();
    for j = 1:4

        group_prog_tbl.([Groups{i}, Progs{j},'_', 'ampInMilliamps']) ...
            = tmp_group.ampInMilliamps(:,j);

        group_prog_tbl.([Groups{i}, Progs{j},'_', 'pulseWidthInMicroseconds']) ...
            = tmp_group.pulseWidthInMicroseconds(:,j);

        group_prog_tbl.([Groups{i}, Progs{j},'_', 'actRechRatio']) ...
            = tmp_group.actRechRatio(:,j);


        group_prog_tbl.([Groups{i}, Progs{j},'_', 'contacts']) ...
            ...
            = cellfun(@(x) ...
            [num2str(x.anodes{j}),'+', num2str(x.cathodes{j}),'-'], tmp_group.contacts, ...
                       'UniformOutput', false);

        if i ==1
            group_prog_tbl.([Progs{j},'_', 'Enabled']) ...
                = tmp_group.validPrograms(:,j);
        end
    end
    tmp_group = removevars(tmp_group,...
            {'ampInMilliamps','pulseWidthInMicroseconds',...
            'actRechRatio', 'contacts', 'validPrograms', 'validProgramNames'});

    tmp_group = renamevars(tmp_group, tmp_group.Properties.VariableNames,...
        cellfun(@(x) [Groups{i},'_',x], tmp_group.Properties.VariableNames,...
            'UniformOutput', false));

    stim_tbl = [stim_tbl, tmp_group, group_prog_tbl]; %#ok<AGROW> 
end

i_grp_struct = cellfun(@(x)...
                strcmp(x, stim_tbl.Properties.VariableNames), Groups,...
                'UniformOutput', false);

i_grp_struct = any(vertcat(i_grp_struct{:}));


stim_tbl(:, i_grp_struct) = [];
%% combine all expanded settings into summary table based off of Session name

% make sumary table based on sub-table that has GREATEST N sessions
% allows sub-tbls w/ empty entries to "fit in"

setting_tbls    = {TD_tbl,  fftSett_tbl, power_band_tbl, ...
                   LD_State_tbl, state_delta_tbl,...
                   stim_tbl};


sess_name       = unique([TD_tbl.sess_name;...
                          fftSett_tbl.sess_name; ...
                          power_band_tbl.sess_name; ...
                          LD_State_tbl.sess_name; ...
                          state_delta_tbl.sess_name;...
                          stim_tbl.sess_name]);

tmp_tbl          = table(sess_name);

% the first variable of every tbl is sess_name

stim_vars       = stim_tbl.Properties.VariableNames(2:end);
td_vars         = TD_tbl.Properties.VariableNames(2:end);
fft_vars        = fftSett_tbl.Properties.VariableNames(2:end);
pwr_vars        = power_band_tbl.Properties.VariableNames(2:end);
ld_vars         = LD_State_tbl.Properties.VariableNames(2:end);
state_vars      = state_delta_tbl.Properties.VariableNames(2:end);


% from sub-tbls, determine variable names and class (cell, struct, double,  etc)
% to build summary table
comp_sense_state_vars       =  {td_vars, fft_vars, pwr_vars, ld_vars, state_vars, stim_vars};

sense_state_class      = cellfun(@(x) varfun(@class, x(:,2:end), 'OutputFormat','cell'), ...
                                 setting_tbls, 'UniformOutput', false);

sense_state_class      = [sense_state_class{:}];
exp_sense_state_vars         = [comp_sense_state_vars{:}];

% initalize and then fill
td_fft_pm_ld_state_stim_tbl = [tmp_tbl, ...
                      table('Size', [height(tmp_tbl), length(exp_sense_state_vars)],...
                            'VariableTypes', sense_state_class,...
                            'VariableNames',exp_sense_state_vars)];

% from sub-tbl pull all variables that share the sess name w/ summary tbl
for i =1:length(setting_tbls)

    sett_tbl = setting_tbls{i};
    i_entry  = find(ismember(tmp_tbl.sess_name, sett_tbl.sess_name));

    td_fft_pm_ld_state_stim_tbl(i_entry, comp_sense_state_vars{i}) = sett_tbl(:, comp_sense_state_vars{i});

end

%% organize/clean-up summary table
% fill in empty entries explicitly for easier parsing going forward
for i=1:length(exp_sense_state_vars)

    if iscell(td_fft_pm_ld_state_stim_tbl{:,exp_sense_state_vars(i)})

        i_empty = find(cellfun(@isempty,...
                          td_fft_pm_ld_state_stim_tbl{:,exp_sense_state_vars(i)}...
                         ));

       td_fft_pm_ld_state_stim_tbl{i_empty,exp_sense_state_vars(i)} = {'Empty Entry'};

    end
end

% bring in session path, timeStart, timeStop, and duration
i_raw_tbl = ismember(db_RCSXXX.sess_name, td_fft_pm_ld_state_stim_tbl.sess_name);

td_fft_pm_ld_state_stim_tbl = [td_fft_pm_ld_state_stim_tbl(:, 'sess_name'),...
                            db_RCSXXX(i_raw_tbl, {'path', 'timeStart', 'timeStop','duration'}),...
                            td_fft_pm_ld_state_stim_tbl(:, 2:end)];


i_emp = cellfun(@isempty,...
                td_fft_pm_ld_state_stim_tbl.timeStart);

td_fft_pm_ld_state_stim_tbl(i_emp,:) = [];

% earliest timeStart minus latest timeStop
td_fft_pm_ld_state_stim_tbl.timeStart = cellfun(@(x) x(1),   td_fft_pm_ld_state_stim_tbl.timeStart);
td_fft_pm_ld_state_stim_tbl.timeStop  = cellfun(@(x) x(end), td_fft_pm_ld_state_stim_tbl.timeStop);

td_fft_pm_ld_state_stim_tbl.duration = td_fft_pm_ld_state_stim_tbl.timeStop - td_fft_pm_ld_state_stim_tbl.timeStart;


td_fft_pm_ld_state_stim_tbl = sortrows(td_fft_pm_ld_state_stim_tbl, 'timeStart');
%%%
%%
%%%%%%%%%%%
%{
 with TD, FFT, power-bands, LD, and aDBS States, aligned and parsed from streamming sessions

   --> flesh-out aDBS offline session parameters, from previous streaming session
%}
%%%%%%%%%%%%
%% from nearest previous streaming session to INS log entry

% more explicit naming
proc_app_log = renamevars(proc_app_log,  'time', 'time_INS');

tmp_ss_tbl   = td_fft_pm_ld_state_stim_tbl(...
                   strcmp(td_fft_pm_ld_state_stim_tbl.activeGroup, 'D'), :);
for i = 1 : height(proc_app_log)

    t_diff   =  proc_app_log.time_INS(i) - tmp_ss_tbl.timeStop;

    % find the smallest negative time difference to get nearest
    % stimLog.json settings BEFORE AppLog.txt settings

    near_t  = min(t_diff(t_diff > 0));

    % need to have streaming session occuring BEFORE INS log time
    if ~isempty(near_t)
        j       = find(t_diff == near_t);
        j       = j(1);

        proc_app_log.prev_sess_name(i)    = tmp_ss_tbl.sess_name(j);
    else
         proc_app_log.prev_sess_name(i) = {'No Entry'};
    end
end

%% explitcly save streaming sessions and offline aDBS sessions according to unique
% sensing parameters

var_oi  = {'chanFullStr', ...
         'bandFormationConfig', 'interval', 'size','streamSizeBins','streamOffsetBins', 'windowLoad',...
         'powerBandinHz', 'powerBinInHz','LD0', 'LD1', 'therapyStatus'};


u_aDBS_Settings = exp_sense_state_vars(...
                      contains(exp_sense_state_vars, var_oi));


sess_w_same_settings  = findgroups(td_fft_pm_ld_state_stim_tbl(:,u_aDBS_Settings));


td_fft_pm_ld_state_stim_tbl = addvars(td_fft_pm_ld_state_stim_tbl, sess_w_same_settings,...
                          'After', 'sess_name');

proc_app_log.sess_w_same_settings =nan(height(proc_app_log),1);

for i = 1:length(unique(td_fft_pm_ld_state_stim_tbl.sess_w_same_settings))

    i_same_sett = find(contains(proc_app_log.prev_sess_name, ...
                  td_fft_pm_ld_state_stim_tbl.sess_name(...
                  td_fft_pm_ld_state_stim_tbl.sess_w_same_settings == i)...
                         ));

    proc_app_log.sess_w_same_settings(i_same_sett) = i;


end

proc_app_log = proc_app_log(~isnan(proc_app_log.sess_w_same_settings),:);
% only save streaming sessions (and their expanded parameters) of aDBS
% offline sessions
i_entry      = cellfun(@(x) ~isempty(x), proc_app_log.prev_sess_name);

i_sess_oi    = find(contains(td_fft_pm_ld_state_stim_tbl.sess_name, ...
                          unique(proc_app_log.prev_sess_name(i_entry))));

app_ss_tbl   = td_fft_pm_ld_state_stim_tbl(i_sess_oi, :);




% app_SS_tbl.RCS02R                   = app_ss_tbl;
% proc_App_log.RCS02R                 = proc_app_log;
% TD_FFT_PB_LD_State_stim_tbl.RCS02R  = td_fft_pm_ld_state_stim_tbl;


end