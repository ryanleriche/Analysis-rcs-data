function  [sense_meta, by_pb_meta, ld0_meta, ld1_meta]...
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

sense_meta = sprintf(['TD samp rate    | %.0f Hz',newline,...
                      'FFT interval    | %.0f ms', newline,...
                      'FFT size        | %.0f samples', newline,...
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

plt_ld_vars = ld0_vars(plt_ss_tbl_oi{:, ld0_vars} > 0);


ld_meta_cell    = table2cell(plt_ss_tbl_oi(1,  plt_ld_vars));

by_ld_meta      = cellfun(@(x,y) [x(5:end), ' | ', num2str(y)], ...
                          plt_ld_vars, ld_meta_cell,...
                          'UniformOutput',false);

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
    

plt_ld_vars = ld1_vars(plt_ss_tbl_oi{:, ld1_vars} > 0);


ld_meta_cell    = table2cell(plt_ss_tbl_oi(1,  plt_ld_vars));

by_ld_meta      = cellfun(@(x,y) [x(5:end), '    | ', num2str(y)], ...
                          plt_ld_vars, ld_meta_cell,...
                          'UniformOutput',false);

str_form = repmat(['%s', newline], 1, length(by_ld_meta));

if isempty(by_ld_meta)
    ld1_meta = 'LD1 was not used';
else
    ld1_meta = sprintf(['LD1\n',str_form], by_ld_meta{:});
end

end