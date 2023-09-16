function ...
    [par_db_RCSXXX, exp_sense_state_vars]...
    ...
    = makeParsedDatabaseRCS(...
    ...
    cfg, pt_side_id, db)
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


%%% uncomment to trouble-shoot as script:
% cfg                        = [];
% cfg.ignoreold_par_db              = true;
% cfg.raw_dir                = [pia_dir, 'raw/'];
% cfg.proc_dir               = [pia_dir, 'processed/parsed_databases/'];
% 
% cfg.pt_sides         = {'RCS02R', 'RCS04R','RCS04L', 'RCS05L','RCS05R',...
%                         'RCS06R','RCS06L', 'RCS07L', 'RCS07R'};


%%
cfg.proc_dir         = [cfg.proc_dir, '/parsed_databases/'];

if ~exist(cfg.proc_dir, 'dir');    mkdir(cfg.proc_dir);    end
    


db_RCSXXX             = db.(pt_side_id);

%% option to add new streaming sessions rather than building from scratch

outputFileName    = fullfile(cfg.proc_dir,[pt_side_id '_parsed_database.mat']);


if exist(outputFileName, 'file') && ~cfg.ignoreold_par_db

    fprintf('%s | Loading previously saved parsed database\n', pt_side_id);

    old = load(outputFileName,'par_db_RCSXXX');

    % of those databases w/ session names, which of them are already parsed
    i_parsed_sess         = contains(db_RCSXXX.sess_name(...
                                        ~cellfun(@isempty,  db_RCSXXX.sess_name)),...
                                        old.par_db_RCSXXX.sess_name);

    % remove parsed sessions to expedient processing
    db_RCSXXX(i_parsed_sess(~cellfun(@isempty,  db_RCSXXX.sess_name)),:) = [];

    db_RCSXXX(cellfun(@isempty,  db_RCSXXX.timeStart),:) = [];
    if isempty(db_RCSXXX)

        fprintf('%s | No new data since %s!\n        --> Existing database returned\n', pt_side_id, old.par_db_RCSXXX.timeStart(end));

        par_db_RCSXXX     = old.par_db_RCSXXX;
        
        % from sess_name -> duration onwards
        exp_sense_state_vars  = old.par_db_RCSXXX.Properties.VariableNames(5:end);
        return

    else
        fprintf('%s | new RCS folders from %s --> %s\n', pt_side_id, db_RCSXXX.timeStart{1}, db_RCSXXX.timeStart{end});

    end

else

    fprintf('%s | Compiling parsed database from scratch\n', pt_side_id)

end


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
LD0_tbl     = parse_LDtbl(tmp_LD_tbl.Ld0);
LD1_tbl     = parse_LDtbl(tmp_LD_tbl.Ld1);


% introduce LD0/1 names into table variable names themselves
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

for d=1:4
    state_deltas = struct2table(cellfun(@(x) x(d), tmp_State_tbl.deltas));

    delta = renamevars(state_deltas,...
               {'fall', 'rise'},...
               {['fallInMilliAmpsPerSecProg', num2str(d-1)],['riseInMilliAmpsPerSecProg', num2str(d-1)]});


    delta_tbl = [delta_tbl, delta]; %#ok<AGROW> 
end

state_amps     = struct2table(tmp_State_tbl.states);
state_i_Prog_j = table;

for d=1:8

    tmp_state = state_amps(:, sprintf('state%g_AmpInMilliamps', d-1));

    tmp_sess_by_prog_amp = tmp_state.Variables;
   
    for j=1:4
        state_i_Prog_j.(sprintf('state%g_AmpInMilliampsProg%g', d-1,j-1)) ...
            = tmp_sess_by_prog_amp(:,j);
    end
end

state_delta_tbl = [db_RCSXXX(i_adapStimSett, 'sess_name'), state_i_Prog_j, delta_tbl];
%% pull shank anatomical targets, and include in all contact/channel naming
i_meta       = cellfun(@(x) ~isempty(x), db_RCSXXX.metaData);
tmp_meta     = cellfun(@(x) x(end,:), db_RCSXXX.metaData(i_meta),...
                             'UniformOutput',false);


tmp_meta_tbl = struct2table(...
                            vertcat(tmp_meta{:})...
                            );


tmp_lead_tar     = [tmp_meta_tbl.leadTargets{:}];


con_0_3_anatomy  = unique(tmp_lead_tar(1, ...
                               cellfun(@(x) ~strcmp(x, ' '), tmp_lead_tar(1,:))...
                          ), 'stable');

con_8_11_anatomy = unique(tmp_lead_tar(3, ...
                               cellfun(@(x) ~strcmp(x, ' '), tmp_lead_tar(3,:))...
                         ), 'stable');

%%% RCS07 intended to have bilateral ACC, but had bilateral subgenual
%%% cingulate cortex leads instead
if contains(pt_side_id, 'RCS07')
    con_8_11_anatomy = {'SGC'};
end

con_8_11_anatomy = [pt_side_id(end),  con_8_11_anatomy{end}];
con_0_3_anatomy  = [pt_side_id(end), con_0_3_anatomy{end}];



%% expand TD settings
i_tdSett   = cellfun(@(x) ~isempty(x), db_RCSXXX.timeDomainSettings);
tmp_td     = cellfun(@(x) x(end,:), db_RCSXXX.timeDomainSettings(i_tdSett),...
                             'UniformOutput',false);

tmp_td_tbl = removevars(vertcat(tmp_td{:}),...
              {'recNum', 'duration', 'timeStart', 'timeStop', 'samplingRate'});

TD_tbl     = db_RCSXXX(i_tdSett,'sess_name');
for d=1:4

    ch_i_td_tbl = struct2table(cellfun(@(x) x(d), tmp_td_tbl.TDsettings));

    ch_i_td_tbl = renamevars(ch_i_td_tbl, ch_i_td_tbl.Properties.VariableNames,...
                     cellfun(@(x) sprintf('Ch%0.0f_%s', d-1, x),...
                            ch_i_td_tbl.Properties.VariableNames,...
                            'UniformOutput', false));

    TD_tbl      =  [TD_tbl, ch_i_td_tbl]; %#ok<AGROW> 

end
td_vars          = TD_tbl.Properties.VariableNames;
to_add_brain_ROI = td_vars(contains(td_vars, {'chanOut', 'chanFullStr'}));


for i_var = 1 : length(to_add_brain_ROI)

    td_var =  to_add_brain_ROI{i_var};

    if contains(td_var, {'0', '1'}) % channels 0 and 1 MUST be contacts 0->3


            TD_tbl.(td_var) = cellfun(@(x) ...
                                        [con_0_3_anatomy, ' ' x], TD_tbl.(td_var),...
                                        'UniformOutput', false);

    elseif contains(td_var, {'2', '3'}) % channels 2 and 3 MUST be contacts 8->11

        TD_tbl.(td_var) = cellfun(@(x) ...
                                        [con_8_11_anatomy, ' ' x], TD_tbl.(td_var),...
                                        'UniformOutput', false);

    end
end


% remove sampleRate as every TD channel MUST have same sampleRate and this
% is parsimonioulsy shown in the FFT settings
TD_tbl = removevars(TD_tbl, td_vars(contains(td_vars, 'sampleRate')));

TD_tbl = movevars(TD_tbl, {'Ch1_chanFullStr', 'Ch2_chanFullStr', 'Ch3_chanFullStr'}, 'After', 'Ch0_chanFullStr');
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

if size(tmp_pwrband_tbl.fftBins,2) ~= 1

    tmp_pwrband_tbl.fftBins = num2cell(tmp_pwrband_tbl.fftBins,2);

end

tmp_band_tbl  = table();
i_ch            = [0 0 1 1 2 2 3 3];
for d =1:8

    powerBand_i = cellfun(@(x) x(d), tmp_pwrband_tbl.powerBandsInHz);
    ch_n_pb_lbl = sprintf('Ch%0.0f_powerBandInHz%0.0f', i_ch(d),d-1);

    tmp_band_tbl(:,ch_n_pb_lbl) = powerBand_i;


    powerBin_i  = cellfun(@(x) x(d), tmp_pwrband_tbl.powerBinsInHz);
    ch_n_pb_lbl = sprintf('Ch%0.0f_powerBinInHz%0.0f', i_ch(d),d-1);

    tmp_band_tbl(:,ch_n_pb_lbl) = powerBin_i;

end
% note power_band_tbl goes:
% Sess Name -> power band setting -> corresponding fft settings

power_band_tbl = [tmp_pwr_tbl(:, 'sess_name'), tmp_band_tbl];
%% expand StimSettings
i_stimSett   = cellfun(@(x) ~isempty(x), db_RCSXXX.stimSettingsOut);

tmp_stim_tbl = db_RCSXXX.stimSettingsOut(i_stimSett);
tmp_stim_tbl = cellfun(@(x) x(end,:), tmp_stim_tbl, 'UniformOutput', false);


stim_tbl     = [db_RCSXXX(i_stimSett, {'sess_name'}),...
                vertcat(tmp_stim_tbl{:})];

i_entry    = find(~cellfun(@isempty, stim_tbl.GroupA));

Groups              = {'GroupA', 'GroupB', 'GroupC', 'GroupD'};
Progs               = {'Prog0', 'Prog1', 'Prog2', 'Prog3'};

i_rmv = [];
% replace empty Groups w/ previous non-empty group
for d=1:height(stim_tbl)
    if isempty(stim_tbl{d, 'GroupA'}{1})

        emp_diff = (d - i_entry);

        [~, i_prev] = min(emp_diff(emp_diff>0));

        if ~isempty(i_prev)
            stim_tbl(d,Groups) = stim_tbl(i_entry(i_prev), Groups);
        % if there is no stimSettings entry before, then remove row as
        % stim has yet to be defined
        else
             i_rmv = [i_rmv, d]; %#ok<AGROW> 
        end
    end
end

% just removes inital entries without any Group definitions
stim_tbl(i_rmv, :) = [];

%%%

for d =1:4

    tmp_group = struct2table([stim_tbl{:, Groups{d}}{:}]);

    group_prog_tbl = table();
    for j = 1:4

        group_prog_tbl.([Groups{d}, Progs{j},'_', 'ampInMilliamps']) ...
            = tmp_group.ampInMilliamps(:,j);

        group_prog_tbl.([Groups{d}, Progs{j},'_', 'pulseWidthInMicroseconds']) ...
            = tmp_group.pulseWidthInMicroseconds(:,j);

        group_prog_tbl.([Groups{d}, Progs{j},'_', 'contacts']) ...
            ...
            = cellfun(@(x) ...
            ['+',num2str(x.anodes{j}), '-', num2str(x.cathodes{j})], tmp_group.contacts, ...
                       'UniformOutput', false);

        if d ==1
            group_prog_tbl.([Progs{j},'_', 'Enabled']) ...
                = tmp_group.validPrograms(:,j);
        end
    end

    group_prog_tbl.([Groups{d},'_', 'actRechRatio']) ...
            = tmp_group.actRechRatio(:,1);

    tmp_group = removevars(tmp_group,...
            {'actRechRatio','ampInMilliamps','pulseWidthInMicroseconds',...
             'contacts', 'validPrograms', 'validProgramNames'});

    tmp_group = renamevars(tmp_group, tmp_group.Properties.VariableNames,...
        cellfun(@(x) [Groups{d},'_',x], tmp_group.Properties.VariableNames,...
            'UniformOutput', false));

    stim_tbl = [stim_tbl, tmp_group, group_prog_tbl]; %#ok<AGROW> 
end

i_grp_struct = cellfun(@(x)...
                strcmp(x, stim_tbl.Properties.VariableNames), Groups,...
                'UniformOutput', false);

i_grp_struct = any(vertcat(i_grp_struct{:}));


stim_tbl(:, i_grp_struct) = [];

%%% add anatomy
stim_vars          = stim_tbl.Properties.VariableNames;
to_add_brain_ROI   = stim_vars(contains(stim_vars , {'contacts'}));

for i_var = 1 : length(to_add_brain_ROI)

    stim_var =  to_add_brain_ROI{i_var};
    
    % adding in side + region for unambiguous contacts when comparing both sides
    tmp          = stim_tbl.(stim_var);
    
    i_con        = cellfun(@(x) ~strcmp(x, '+-') , tmp);
    
    if ~isempty(tmp(i_con))
        ind_contacts = cellfun(@(x) regexp(x,'\d*','Match'), tmp(i_con), 'UniformOutput', false);
        ind_contacts = cellfun(@(x) x(2), ind_contacts);
        
        i_small      = cellfun(@(x) any(strcmp(x, ["0","1","2","3"])), ind_contacts);
        i_large      = cellfun(@(x) any(strcmp(x, ["8","9","10","11"])), ind_contacts);
        
        i_con        = find(i_con);
        
        
        tmp(i_con(i_small)) = cellfun(@(x) [con_0_3_anatomy, ' ', x], ...
                                            tmp(i_con(i_small)), 'UniformOutput', false);
        
        tmp(i_con(i_large)) = cellfun(@(x) [con_8_11_anatomy, ' ', x], ...
                                            tmp(i_con(i_large)), 'UniformOutput', false);
    end    
        % rather than 16, report monopolar stim explictly as case
        stim_tbl.(stim_var) = replace(tmp, '16', 'c');

end


%% from TimeSync.json, return mean INS latency per session (INS time - API time)

i_TimeSync          = find(cellfun(@isstruct, db_RCSXXX.INS_API_latency));


INS_lat_tbl = [db_RCSXXX(i_TimeSync, {'sess_name'})...
                        ...
              renamevars(struct2table([db_RCSXXX.INS_API_latency{i_TimeSync}]),...
              {'mean', 'std', 'max', 'min'},...
               {'INS_lat_mean', 'INS_lat_std', 'INS_lat_max', 'INS_lat_min'})];


%% combine all expanded settings into summary table based off of Session name

% make sumary table based on sub-table that has GREATEST N sessions
% allows sub-tbls w/ empty entries to "fit in"

setting_tbls    = {TD_tbl,  fftSett_tbl, power_band_tbl, ...
                   LD_State_tbl, state_delta_tbl,...
                   stim_tbl, INS_lat_tbl};


sess_name       = unique([TD_tbl.sess_name;...
                          fftSett_tbl.sess_name; ...
                          power_band_tbl.sess_name; ...
                          LD_State_tbl.sess_name; ...
                          state_delta_tbl.sess_name;...
                          stim_tbl.sess_name;...
                          INS_lat_tbl.sess_name]);

tmp_tbl          = table(sess_name);

% the first variable of every tbl is sess_name

stim_vars       = stim_tbl.Properties.VariableNames(2:end);
td_vars         = TD_tbl.Properties.VariableNames(2:end);
fft_vars        = fftSett_tbl.Properties.VariableNames(2:end);
pwr_vars        = power_band_tbl.Properties.VariableNames(2:end);
ld_vars         = LD_State_tbl.Properties.VariableNames(2:end);
state_vars      = state_delta_tbl.Properties.VariableNames(2:end);

INS_lat_vars    = INS_lat_tbl.Properties.VariableNames(2:end); 

% from sub-tbls, determine variable names and class (cell, struct, double,  etc)
% to build summary table
comp_sense_state_vars       =  {td_vars, fft_vars, pwr_vars, ld_vars, ...
                                state_vars, stim_vars, INS_lat_vars};

sense_state_class      = cellfun(@(x) varfun(@class, x(:,2:end), 'OutputFormat','cell'), ...
                                 setting_tbls, 'UniformOutput', false);

sense_state_class      = [sense_state_class{:}];
exp_sense_state_vars         = [comp_sense_state_vars{:}];

% initalize and then fill
par_db_RCSXXX = [tmp_tbl, ...
                      table('Size', [height(tmp_tbl), length(exp_sense_state_vars)],...
                            'VariableTypes', sense_state_class,...
                            'VariableNames',exp_sense_state_vars)];

% from sub-tbl pull all variables that share the sess name w/ summary tbl
for d =1:length(setting_tbls)

    sett_tbl = setting_tbls{d};
    i_entry  = find(ismember(tmp_tbl.sess_name, sett_tbl.sess_name));

    par_db_RCSXXX(i_entry, comp_sense_state_vars{d}) = sett_tbl(:, comp_sense_state_vars{d});

end

%% organize/clean-up summary table
% fill in empty entries from previous non-empty for easier parsing going forward
for d=1:length(exp_sense_state_vars)

    var_oi = par_db_RCSXXX{:,exp_sense_state_vars(d)};
    if iscell(var_oi)

        tf_empty =  cellfun(@isempty, var_oi);

        i_empty  = find(tf_empty);
        i_fill   = find(~tf_empty);

        for j=1:length(i_empty)

            diffs = i_empty(j) - i_fill;

            [~, i_prev] = min(diffs(diffs>0));

            if ~isempty(i_prev)

                par_db_RCSXXX{i_empty(j),exp_sense_state_vars(d)} ...
                     = ...
                var_oi(i_fill(i_prev));

            end
        end

    elseif isduration(var_oi)

         % for durations of '0:00:00' just put in NaN
         var_oi(eq(var_oi, duration('0:00:00')))   = NaN;

         par_db_RCSXXX{:,exp_sense_state_vars(d)}  = var_oi;

    end
end

%% bring in session path, timeStart, timeStop, and duration
i_raw_tbl = ismember(db_RCSXXX.sess_name, par_db_RCSXXX.sess_name);

par_db_RCSXXX = [par_db_RCSXXX(:, 'sess_name'),...
                            db_RCSXXX(i_raw_tbl, {'path', 'timeStart', 'timeStop','duration'}),...
                            par_db_RCSXXX(:, 2:end)];


i_emp = cellfun(@isempty,...
                par_db_RCSXXX.timeStart);

par_db_RCSXXX(i_emp,:) = [];

% earliest timeStart minus latest timeStop
% --> respecify TimeZone due to low-level assumption of SAME UTC offset which
%     neglects daylight savings
par_db_RCSXXX.timeStart          = cellfun(@(x) x(1),   par_db_RCSXXX.timeStart);
par_db_RCSXXX.timeStart.TimeZone = 'America/Los_Angeles';

par_db_RCSXXX.timeStop           = cellfun(@(x) x(end), par_db_RCSXXX.timeStop);
par_db_RCSXXX.timeStop.TimeZone  = 'America/Los_Angeles';

par_db_RCSXXX.duration           = par_db_RCSXXX.timeStop - par_db_RCSXXX.timeStart;



%% COMBINE WITH OLD DATABASE
% IF the old database existed, recombine with new database and sort it
% but first fix cell/ mat class issues

 if exist('old', 'var')

      if ~isempty(par_db_RCSXXX)

            par_db_RCSXXX = [old.par_db_RCSXXX; par_db_RCSXXX];
      else
            par_db_RCSXXX = old.par_db_RCSXXX;

            disp(['Database as of ', datestr(par_db_RCSXXX.timeStart(1)), '; ',...
            datestr(par_db_RCSXXX.timeStart(end)),'---no new files since ALL were empty'])
      end
 end

par_db_RCSXXX = sortrows(par_db_RCSXXX, 'timeStart');
  
%%
par_db_file = fullfile(cfg.proc_dir,[pt_side_id, '_parsed_database.xlsx']);
                
if isfile(par_db_file);       delete(par_db_file);      end

writetable(par_db_RCSXXX, par_db_file);

save(fullfile(cfg.proc_dir,[pt_side_id '_parsed_database.mat']),...
    'par_db_RCSXXX');

fprintf('%s | parsed database saved to %s%s_parsed_database.mat\n',...
    pt_side_id, cfg.proc_dir, pt_side_id);



end