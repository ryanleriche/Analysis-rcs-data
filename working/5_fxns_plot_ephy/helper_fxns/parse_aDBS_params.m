function  meta = parse_aDBS_params(plt_ss_tbl_oi)

all_vars = plt_ss_tbl_oi.Properties.VariableNames;


static_across_sess = {'TDsampleRates', ...
                     'fft_intervalInMilliseconds', 'fft_sizeInSamples',...
                     'bandFormationConfig',...
                     'Ch0_chanOut','Ch1_chanOut', 'Ch2_chanOut', 'Ch3_chanOut','GroupDProg0_contacts'};


vars = table2cell(plt_ss_tbl_oi(1,  static_across_sess));

meta.sense = sprintf(['TD samp rate    | %g Hz',newline,...
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
ld0_pb_vars   = all_vars(contains(all_vars, 'LD0_PowerBand_input_'));
i_pb_input    = ~strcmp(plt_ss_tbl_oi{1, ld0_pb_vars}, 'Disabled');

ld0_pb_inputs = plt_ss_tbl_oi{1, ld0_pb_vars}(i_pb_input);

plt_pbs       = cell(1,length(ld0_pb_inputs));

for i=1:length(ld0_pb_inputs)

    switch ld0_pb_inputs{i}

        case 'Ch0Band0'; pb ='Ch0_powerBinInHz0';
        case 'Ch0Band1'; pb ='Ch0_powerBinInHz1';

        case 'Ch1Band0'; pb ='Ch1_powerBinInHz2';
        case 'Ch1Band1'; pb ='Ch1_powerBinInHz3';

        case 'Ch2Band0'; pb ='Ch2_powerBinInHz4';
        case 'Ch2Band1'; pb ='Ch2_powerBinInHz5';

        case 'Ch3Band0'; pb ='Ch3_powerBinInHz6';
        case {'Ch3Band1', 'C32Band1'}; pb ='Ch3_powerBinInHz7';

    end
    plt_pbs{i}= pb;
end
if isempty(plt_pbs)

    meta.by_ld0_pb = 'No Power-bands inputted to LD0';

else
    form_pbs        = cellfun(@(x,y) [x,' | ', y, newline], ...
                    plt_pbs,...
                    table2cell(plt_ss_tbl_oi(1,  plt_pbs)), 'UniformOutput',false);

    meta.by_ld0_pb   = ['LD0 Inputs:', newline,form_pbs{:}];
end
%%
% check LD powerband detection inputs, if band is in either return its bin
% in Hz
ld1_pb_vars   = all_vars(contains(all_vars, 'LD1_PowerBand_input_'));
i_pb_input    = ~strcmp(plt_ss_tbl_oi{1, ld1_pb_vars}, 'Disabled');

ld1_pb_inputs = plt_ss_tbl_oi{1, ld1_pb_vars}(i_pb_input);

plt_pbs       = cell(1,length(ld1_pb_inputs));

for i=1:length(ld1_pb_inputs)

    switch ld1_pb_inputs{i}

        case 'Ch0Band0'; pb ='Ch0_powerBinInHz0';
        case 'Ch0Band1'; pb ='Ch0_powerBinInHz1';

        case 'Ch1Band0'; pb ='Ch1_powerBinInHz2';
        case 'Ch1Band1'; pb ='Ch1_powerBinInHz3';

        case 'Ch2Band0'; pb ='Ch2_powerBinInHz4';
        case 'Ch2Band1'; pb ='Ch2_powerBinInHz5';

        case 'Ch3Band0'; pb ='Ch3_powerBinInHz6';
        case 'Ch3Band1'; pb ='Ch3_powerBinInHz7';

    end
    plt_pbs{i}= pb;
end

if isempty(plt_pbs)

    meta.by_ld1_pb = 'No Power-bands inputted to LD1';

else
    form_pbs        = cellfun(@(x,y) [x,' | ', y, newline], ...
                    plt_pbs,...
                    table2cell(plt_ss_tbl_oi(1,  plt_pbs)), 'UniformOutput',false);

    meta.by_ld1_pb   = ['LD1 Inputs:', newline, form_pbs{:}];
end

%% 
ld0_vars   = all_vars(contains(all_vars, {'LD0'}) & ... 
                   strcmp(varfun(@class, plt_ss_tbl_oi(:,all_vars),...
                        'OutputFormat','cell'), 'double')...
                   );
    
if plt_ss_tbl_oi.LD0_biasTerm0 == plt_ss_tbl_oi.LD0_biasTerm1

    ld0_vars = ld0_vars(~contains(ld0_vars, 'LD0_biasTerm1'));

end

plt_ld_vars            = ld0_vars(plt_ss_tbl_oi{:, ld0_vars} ~=0  |...
                                  contains(ld0_vars, {'fractionalFixed', 'Duration', 'biasTerm'}))';

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

if isempty(by_ld_meta) || all([ld_meta_cell{2,:}] ==0)
    meta.ld0 = 'LD0 was not used';
else
    meta.ld0 = sprintf(['LD0\n',str_form], by_ld_meta{:});
end


%%
ld1_vars  = all_vars(contains(all_vars, {'LD1'}) & ... 
                   strcmp(varfun(@class, plt_ss_tbl_oi(:,all_vars),...
                        'OutputFormat','cell'), 'double')...
                   );

if plt_ss_tbl_oi.LD1_biasTerm0 == plt_ss_tbl_oi.LD1_biasTerm1

    ld1_vars = ld1_vars(~contains(ld1_vars, 'LD1_biasTerm1'));

end
    

plt_ld_vars            = ld1_vars(plt_ss_tbl_oi{:, ld1_vars} ~= 0 |...
                                  contains(ld1_vars, {'fractionalFixed', 'Duration','biasTerm'}))';

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

if isempty(by_ld_meta) || all([ld_meta_cell{2,:}] ==0)

    meta.ld1 = 'LD1 was not used';
else
    meta.ld1 = sprintf(['LD1\n',str_form], by_ld_meta{:});
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

meta.state   = sprintf(['State to Stim\n',str_form],  by_state_meta{:});


end