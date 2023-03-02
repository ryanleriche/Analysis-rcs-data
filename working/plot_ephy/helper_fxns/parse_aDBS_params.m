function  [sense_meta, by_pb_meta, ld0_meta, ld1_meta, state_meta]...
        ...
        = parse_aDBS_params(...
        ...
    plt_ss_tbl_oi)

all_vars = plt_ss_tbl_oi.Properties.VariableNames;


static_across_sess = {'TDsampleRates', ...
                     'fft_intervalInMilliseconds', 'fft_sizeInSamples',...
                     'bandFormationConfig',...
                     'Ch0_chanOut','Ch1_chanOut', 'Ch2_chanOut', 'Ch3_chanOut','GroupDProg0_contacts'};


vars = table2cell(plt_ss_tbl_oi(1,  static_across_sess));

sense_meta = sprintf(['TD samp rate    | %g Hz',newline,...
                      'FFT interval    | %g ms', newline,...
                      'FFT size        | %g samples', newline,...
                      'Bit Shift       | %s', newline,...
                      'TD Ch0             | %s', newline,...
                      'TD Ch1             | %s', newline,...
                      'TD Ch2             | %s', newline,...
                      'TD Ch3             | %s', newline,...
                      'Stim Contacts      | %s'],...
                    vars{:});

%% only return bands used in LDs
% check LD powerband detection inputs, if band is in either return its bin
% in Hz
if all(strcmp(plt_ss_tbl_oi{1, 'LD0_powerband_detectionInputs'}, 'None'))

    nonzero_LD0_pbs = [];

else
    nonzero_LD0_pbs = cellfun(@(x) str2num(x(end)), plt_ss_tbl_oi{1, 'LD0_powerband_detectionInputs'});

end


if all(strcmp(plt_ss_tbl_oi{1, 'LD1_powerband_detectionInputs'}, 'None'))

    nonzero_LD1_pbs = [];

else
    nonzero_LD1_pbs = cellfun(@(x) str2num(x(end)), plt_ss_tbl_oi{1, 'LD1_powerband_detectionInputs'});
end

nonzero_pbs     = unique([nonzero_LD0_pbs, nonzero_LD1_pbs]);

by_chans_pb     =  {'Ch0_powerBinInHz0', 'Ch0_powerBinInHz1',...
                    'Ch1_powerBinInHz2', 'Ch1_powerBinInHz3',...
                    'Ch2_powerBinInHz4', 'Ch2_powerBinInHz5',...
                    'Ch3_powerBinInHz6', 'Ch3_powerBinInHz7'};

plt_pbs         = by_chans_pb(nonzero_pbs+1);

form_pbs        = cellfun(@(x,y) [x,' | ', y, newline], ...
                    plt_pbs,...
                    table2cell(plt_ss_tbl_oi(1,  plt_pbs)), 'UniformOutput',false);

if isempty(form_pbs)
    by_pb_meta = '';
else
    by_pb_meta   = [form_pbs{:}];
end

%% 
ld0_vars   = all_vars(contains(all_vars, {'LD0'}) & ... 
                   strcmp(varfun(@class, plt_ss_tbl_oi(:,all_vars),...
                        'OutputFormat','cell'), 'double')...
                   );
    
if plt_ss_tbl_oi.LD0_biasTerm0 == plt_ss_tbl_oi.LD0_biasTerm1

    ld0_vars = ld0_vars(~contains(ld0_vars, 'LD0_biasTerm1'));

end

plt_ld_vars            = ld0_vars(plt_ss_tbl_oi{:, ld0_vars} > 0 |...
                                  contains(ld0_vars, {'fractionalFixed', 'Duration'}))';

tmp_meta_cell          = table2cell(plt_ss_tbl_oi(1,  plt_ld_vars));

ld_meta_cell          = [tmp_meta_cell; tmp_meta_cell];
by_ld_meta            = cell(length(tmp_meta_cell),1);

% for LD updateRate, blanking, onset duration, and termination duration
% include in seconds rather than N of FFT intervals or N updateRates
% (pg. 88-89 of the 4NR010 Research Lab Programmer Guide M979053A001)

for i=1:length(plt_ld_vars)
    switch plt_ld_vars{i}

        case {'LD0_updateRate', 'LD0_blankingDurationUponStateChange'}
            ld_meta_cell{2,i}  = tmp_meta_cell{i} * (plt_ss_tbl_oi.fft_intervalInMilliseconds / 1000);

        case {'LD0_onsetDuration', 'LD0_terminationDuration'}
            if tmp_meta_cell{i} > 0
                ld_meta_cell{2,i} = tmp_meta_cell{i} * plt_ss_tbl_oi.LD0_updateRate * (plt_ss_tbl_oi.fft_intervalInMilliseconds / 1000);
            else
                ld_meta_cell{2,i} = plt_ss_tbl_oi.LD0_updateRate * (plt_ss_tbl_oi.fft_intervalInMilliseconds / 1000);
            end
    end

    if ld_meta_cell{1,i} ~= ld_meta_cell{2,i}

        by_ld_meta{i} = sprintf('    %s | %g (%g seconds)', plt_ld_vars{i}(5:end), ld_meta_cell{1,i}, ld_meta_cell{2,i});

    else
        by_ld_meta{i} = sprintf('    %s | %g', plt_ld_vars{i}(5:end), ld_meta_cell{1,i});

    end
end

str_form = repmat(['%s', newline], 1, length(by_ld_meta));

if isempty(by_ld_meta)
    ld0_meta = 'LD0 was not used';
else
    ld0_meta = sprintf(['LD0\n',str_form], by_ld_meta{:});
end


%%
ld1_vars  = all_vars(contains(all_vars, {'LD1'}) & ... 
                   strcmp(varfun(@class, plt_ss_tbl_oi(:,all_vars),...
                        'OutputFormat','cell'), 'double')...
                   );

if plt_ss_tbl_oi.LD1_biasTerm0 == plt_ss_tbl_oi.LD1_biasTerm1

    ld1_vars = ld1_vars(~contains(ld1_vars, 'LD1_biasTerm1'));

end
    

plt_ld_vars            = ld1_vars(plt_ss_tbl_oi{:, ld1_vars} > 0 |...
                                  contains(ld1_vars, {'fractionalFixed', 'Duration'}))';

tmp_meta_cell         = table2cell(plt_ss_tbl_oi(1,  plt_ld_vars));

ld_meta_cell          = [tmp_meta_cell; tmp_meta_cell];
by_ld_meta            = cell(length(tmp_meta_cell),1);

% for LD updateRate, blanking, onset duration, and termination duration
% include in seconds rather than N of FFT intervals or N updateRates
% (pg. 88-89 of the 4NR010 Research Lab Programmer Guide M979053A001)

for i=1:length(plt_ld_vars)
    switch plt_ld_vars{i}

        case {'LD1_updateRate', 'LD1_blankingDurationUponStateChange'}
            ld_meta_cell{2,i}  = tmp_meta_cell{i} * (plt_ss_tbl_oi.fft_intervalInMilliseconds / 1000);

        case {'LD1_onsetDuration', 'LD1_terminationDuration'}
            
            if tmp_meta_cell{i} > 0
                ld_meta_cell{2,i} = tmp_meta_cell{i} * plt_ss_tbl_oi.LD1_updateRate * (plt_ss_tbl_oi.fft_intervalInMilliseconds / 1000);
            else
                ld_meta_cell{2,i} = plt_ss_tbl_oi.LD1_updateRate * (plt_ss_tbl_oi.fft_intervalInMilliseconds / 1000);
            end
    end

    if ld_meta_cell{1,i} ~= ld_meta_cell{2,i}

        by_ld_meta{i} = sprintf('    %s | %g (%g seconds)', plt_ld_vars{i}(5:end), ld_meta_cell{1,i}, ld_meta_cell{2,i});

    else
        by_ld_meta{i} = sprintf('    %s | %g', plt_ld_vars{i}(5:end), ld_meta_cell{1,i});

    end
end

str_form = repmat(['%s', newline], 1, length(by_ld_meta));

if isempty(by_ld_meta)
    ld1_meta = 'LD1 was not used';
else
    ld1_meta = sprintf(['LD1\n',str_form], by_ld_meta{:});
end

%% return state meta data

state_vars   = all_vars(contains(all_vars, {'state', 'rise', 'fall', 'GroupD'}) & ~contains(all_vars, {'contacts', 'GroupDProg0_ampInMilliamps'}));
prog_enabled = num2cell(....
                       find(...
                       plt_ss_tbl_oi{:,...
                           {'Prog0_Enabled', 'Prog1_Enabled', ...
                           'Prog2_Enabled', 'Prog3_Enabled'}}) -1);


state_txt     = [cellfun(@(x) ['Prog', num2str(x)], prog_enabled, 'UniformOutput', false),...
                 ];

state_vars   = state_vars(contains(state_vars, state_txt))';

state_vars   = state_vars(...
                    plt_ss_tbl_oi{:, state_vars} >=0);

% remove all programs w/ pulseWidth of 850 --> not really enabled
for i = 0:3

    i_var_prog = contains(state_vars, sprintf('Prog%g', i));

    if any(i_var_prog)

        tmp_var = state_vars(i_var_prog);
        if plt_ss_tbl_oi{:, tmp_var(contains(tmp_var, 'pulseWidth'))} ==850

        state_vars(i_var_prog) = [];

        end
    end
end

by_state_meta = cell(length(state_vars), 1);

for i = 1:length(state_vars)
     by_state_meta{i} = sprintf('    %s | %g', state_vars{i}, plt_ss_tbl_oi{:, state_vars{i}});
end


str_form = repmat(['%s', newline], 1, length(state_vars));

state_meta   = sprintf(['State to Stim\n',str_form],  by_state_meta{:});


end