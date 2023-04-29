function [N_days_back_tbl, sham_vs_stim_stats] = plot_RCS02_stage3(cfg, s3_tbl, stim_log)

% cfg= [];
% cfg.proc_dir          = [pia_dir, 'processed/pain_per_DBS_parameters/'];
% cfg.proc_subdir       = 'Stage 3 (ol-DBS_versus_sham)';
% cfg.plt_metrics       = {'painVAS', 'mayoNRS',  'MPQtotal', 'unpleasantVAS'};
% cfg.pt_id             = 'RCS02';
% 
% cfg.pt_meta           = pt_META;
% 
% cfg.N_days            = duration('3:00:00:00');
% 
% s3_tbl                = stimGroups.RCS02.s3{1};
%%
%% account for when RCS02 went on vacation
pt_meta      = cfg.pt_meta.(cfg.pt_id);


break_starts = pt_meta.dates(contains(pt_meta.desc, 'break_start'));
break_stops  = pt_meta.dates(contains(pt_meta.desc, 'break_stop'));


%% pull stimLog to get precise times of how long RCS02 was on a given stim group 
% before opting to try another

s3_start_date = pt_meta.dates(...
                    strcmp(...
                        pt_meta.desc, ...
                        's3_start_ol_versus_sham'));

s3_stop_date = pt_meta.dates(...
                    strcmp(...
                        pt_meta.desc, ...
                        's3_stop_ol_versus_sham'));

s3_stim_log = stim_log(...
                        ge(stim_log.time_stimLog, s3_start_date) &...
                        le(stim_log.time_stimLog, s3_stop_date),...
                        :);

s3_stim_log.break_or_s3   = repmat({'s3'}, height(s3_stim_log),1);

for j = 1:length(break_starts)

    i_break = find(ge(s3_stim_log.time_stimLog, break_starts(j)) &...
                   le(s3_stim_log.time_stimLog, break_stops(j)));

   
    s3_stim_log.break_or_s3(i_break) = {'break'};

end
%%
plt_vars    = {'stimContacts','ampInMilliamps','rateInHz','pulseWidthInMicroseconds',...
               'percentDutyCycle','therapyStatusDescription'};


s3_stim_log.u_stim_params = findgroups(s3_stim_log(:, plt_vars));

i_Off                    = find(strcmp(s3_stim_log.therapyStatusDescription, 'Off'));

s3_stim_log.u_stim_params(i_Off) = max(s3_stim_log.u_stim_params) + 1;


diff_params   = [1;...
                 diff(s3_stim_log.u_stim_params) ~=0 ...
                 ];

i_diffs  = find(diff_params);
i_starts = [];           i_ends   = [];

for h=1:length(i_diffs)-1
    i_starts(h) = i_diffs(h); %#ok<*AGROW> 
    i_ends(h)   = i_diffs(h+1)+1;
end

i_starts = [i_starts, i_ends(end)-1];
i_ends   = [i_ends,height(s3_stim_log)];

for h=1:length(i_starts)

    s3_stim_log.stim_params_by_days(i_starts(h):i_ends(h))  = h;
end








%%





sham_vs_stim_stats = struct;


tmp_s3_tbl = s3_tbl;


for j = 1:length(break_starts)

    i_break = find(ge(tmp_s3_tbl.time, break_starts(j)) &...
                le(tmp_s3_tbl.time, break_stops(j)));

    tmp_s3_tbl(i_break , :) = [];

end

%%
% plt_vars    = {'R_stimContacts','R_ampInMilliamps','R_rateInHz','R_pulseWidthInMicroseconds',...
%                'R_percentDutyCycle','R_therapyStatusDescription'};
% 
% 
% tmp_s3_tbl.u_stim_params = findgroups(tmp_s3_tbl(:, plt_vars));
% 
% i_Off                    = find(strcmp(tmp_s3_tbl.R_therapyStatusDescription, 'Off'));
% 
% tmp_s3_tbl.u_stim_params(i_Off) = max(tmp_s3_tbl.u_stim_params) + 1;
% 
% 
% diff_params   = [1;...
%                  diff(tmp_s3_tbl.u_stim_params) ~=0 ...
%                  ];
% 
% i_diffs  = find(diff_params);
% i_starts = [];           i_ends   = [];
% 
% for h=1:length(i_diffs)-1
%     i_starts(h) = i_diffs(h); %#ok<*AGROW> 
%     i_ends(h)   = i_diffs(h+1)+1;
% end
% 
% i_starts = [i_starts, i_ends(end)-1];
% i_ends   = [i_ends,height(tmp_s3_tbl)];
% 
% for h=1:length(i_starts)
% 
%     tmp_s3_tbl.stim_params_by_days(i_starts(h):i_ends(h))  = h;
% end
% 


for i = 1:length(cfg.plt_metrics)
    params_by_pain_cell = cell(length(i_starts),1);
    params_per_se_cell  = cell(length(i_starts),1);
    

    for h=1:length(i_starts)
    
        i_params                = find(tmp_s3_tbl.stim_params_by_days == h);
        params_by_pain_cell{h}  = tmp_s3_tbl{i_params, cfg.plt_metrics{i}};
         
        % return stim parameters over epoch
        tmp_vars = [(table2cell(...
                        unique(tmp_s3_tbl(i_params(1), plt_vars)), 'rows')),...
                     cellstr(...
                        datestr(...
                             dateshift(tmp_s3_tbl.time(i_params([1,end])), 'start', 'day'), 'mmm-dd')...
                        )'...
                     ];
    
        params_per_se_cell{h} = sprintf(['%s\\newline',...
                                         '%g mA\\newline',...
                                         '%g Hz\\newline',...
                                         '%g \\mus\\newline',...
                                         '%g%% DC\\newline',...
                                         '%s\\newline',...
                                         '%s\\newline%s'], tmp_vars{:});         
    end
    
    pain_by_params               = padcat(params_by_pain_cell{:});
    %%% plot everytime a stim parameter was changed
    save_dir = [cfg.proc_dir,cfg.pt_id,'/',cfg.proc_subdir,'/',cfg.plt_metrics{i},'/'];
    
    if ~isfolder(save_dir);     mkdir(save_dir); end
    
    figure('Units', 'Inches', 'Position', [0, 0, 30, 18]);
    sgtitle([cfg.pt_id, newline, 'day-by-day stage 3'], 'Fontsize',20);
    
    N_per_row = ceil(length(pain_by_params)/3);

    i_epochs = {1:N_per_row, N_per_row+1:N_per_row*2, N_per_row*2+1:length(pain_by_params)};
    
    for h=1:length(i_epochs)
        subplot(3,1,h)
        
        UnivarScatter(pain_by_params(:,i_epochs{h}),...
                        'Compression',50,... distance btwn boxes
                        'Width',1,...
                        'MarkerFaceColor','k',...
                        'MarkerEdgeColor', 'k',...
                        'PointSize', 8,...
                        'MeanColor', [.5 .5 .5],... color(s) of the mean line; def. black
                        'SEMColor', [.9 .9 .9],...  color(s) of the SEM box; def. dark gray
                        'StdColor',[.9 .9 .9]...
                        );
        
        ylabel(cfg.plt_metrics{i});     
        
        if contains(cfg.plt_metrics{i}, 'VAS');         ylim([0 100]);
        elseif contains(cfg.plt_metrics{i}, 'NRS');     ylim([0 10]);
        elseif contains(cfg.plt_metrics{i}, 'MPQ');     ylim([0 45]);  
        end
        
        xticklabels(params_per_se_cell(i_epochs{h}));      set(gca,'Fontsize', 14);
    end
    
    exportgraphics(gcf, [save_dir, cfg.plt_metrics{i}, '_day-by-day stage 3','.png'])    

    %% now group stim parameters REGARDLESS of when they occured (open-loop stim groups versus sham)
    u_param_days           = unique(tmp_s3_tbl.stim_params_by_days);
    s3_stim_group_cell     = {};
    tmp_s3_tbl.description = cell(height(tmp_s3_tbl), 1);
    for j = 1 : length(u_param_days)
    
        i_params = tmp_s3_tbl.stim_params_by_days == u_param_days(j);
        tmp_tbl  = tmp_s3_tbl(i_params,:);
    
        % ONLY return stim parameters tried for at least 2 days
        if ge(tmp_tbl.time(end) -  tmp_tbl.time(1), duration('3:00:00:00'))
    
            last_day      = dateshift(tmp_tbl.time(end), 'end', 'day');
    
            N_days_before = last_day - cfg.N_days;
    
            s3_stim_group_cell{end+1} = tmp_tbl(ge(tmp_tbl.time, N_days_before),:);
    
            tmp_s3_tbl.description(i_params) = {'s3 stim group'};
    
        else
    
            tmp_s3_tbl.description(i_params) = {'washout OR too short'};
    
        end
    end
    
    N_days_back_tbl = vertcat(s3_stim_group_cell{:});
    
    % edge case where before and after vacation was washout so washout appears
    % like many days w/ stim OFF , but few surveys
    i_Off                     = strcmp(N_days_back_tbl.R_therapyStatusDescription, 'Off');
    N_days_back_tbl(i_Off, :) = [];
    
    
    u_param_lbls = unique(N_days_back_tbl.u_stim_params);
    
    params_by_pain_cell = cell(length(u_param_lbls),1);
    params_per_se_cell  = cell(length(u_param_lbls),1);
    
    
    for h=1:length(u_param_lbls)
    
        i_params                = find(N_days_back_tbl.u_stim_params == u_param_lbls(h));
        params_by_pain_cell{h}  = N_days_back_tbl{i_params, cfg.plt_metrics{i}};
         
        % return stim parameters over epoch
        tmp_vars = table2cell(...
                        unique(N_days_back_tbl(i_params, plt_vars), 'rows'));
    
        params_per_se_cell{h} = sprintf(['%s\\newline',...
                                         '%g mA\\newline',...
                                         '%g Hz\\newline',...
                                         '%g \\mus\\newline',...
                                         '%g%% DC\\newline',...
                                         '%s\\newline',...
                                         '%s'], tmp_vars{1:end-1});         
    end
    
    pain_by_params               = padcat(params_by_pain_cell{:});
    %%% plot open-loop versus sham stim
    
    figure('Units', 'Inches', 'Position', [0, 0, 14, 10]);
    sgtitle(...
        sprintf('%s (stage 3 blinded testing)\nlast %g days for given week per stim group',...
            cfg.pt_id, days(cfg.N_days)),...
        'Fontsize',20);
    
    UnivarScatter(pain_by_params,...
            'Compression',50,... distance btwn boxes
            'Width',1, ...
                'MarkerFaceColor','k',...
                'MarkerEdgeColor', 'k',...
                'PointSize', 8,...
                'MeanColor', [.5 .5 .5],... color(s) of the mean line; def. black
                'SEMColor', [.9 .9 .9],...  color(s) of the SEM box; def. dark gray
                'StdColor',[.9 .9 .9]...
                );
        
    if        contains(cfg.plt_metrics{i}, 'VAS');        ylim([0 100]);
    elseif    contains(cfg.plt_metrics{i}, 'NRS');        ylim([0 10]);
    elseif    contains(cfg.plt_metrics{i}, 'MPQ');        ylim([0 45]);  
    end
    
    xticklabels(params_per_se_cell);      set(gca,'Fontsize', 14); grid on;
    ylabel(cfg.plt_metrics{i}, 'FontSize', 24);     

    exportgraphics(gcf, [save_dir, cfg.plt_metrics{i}, '_stage_3_ol_versus_sham','.png']);
%% NOW group into more parsimonious SHAM vs stim groups
    numeric_vars = {'R_rateInHz','R_ampInMilliamps','R_pulseWidthInMicroseconds',...
                    'R_percentDutyCycle'};

    %%% FIRST identify RACC 1+0- stim explicity, THEN write over w/ Sham
    %%% tag
     N_days_back_tbl.R_stim_groups(...
         strcmp(N_days_back_tbl.R_stimContacts, 'RACC 1+0-'))    = {'RACC 1+0-'};

     N_days_back_tbl.R_stim_groups(...
         strcmp(N_days_back_tbl.R_stimContacts, 'RACC 3+2-'))    = {'RACC 3+2-'};


     N_days_back_tbl.R_stim_groups(N_days_back_tbl.R_rateInHz == 2)    = {'Sham'};
    


    [stim_grps, u_params] = unique(N_days_back_tbl(:, cfg.seperate_by), 'rows');


    pain_params_by_cell     = cell(length(u_params), 1);
    params_per_se_cell      = pain_params_by_cell;

     for j = 1 : length(u_params)

         i_params = find(...
                            ismember(N_days_back_tbl(:, cfg.seperate_by),...
                                     N_days_back_tbl(u_params(j), cfg.seperate_by))...
                                     );

         pain_params_by_cell{j}     = N_days_back_tbl{i_params, cfg.plt_metrics{i}};

        % return stim parameters
        tmp_vars       = N_days_back_tbl{i_params, numeric_vars};
        stimContact_oi = unique(N_days_back_tbl{i_params, 'R_stimContacts'})';


        format_stim_params = [stimContact_oi, ...
                                 cellfun(@(x, y) sprintf('%.2f ± %.2f', x, y),...
                                    num2cell(mean(tmp_vars)),...
                                    num2cell(std(tmp_vars)),...
                                    'UniformOutput', false)...
                                    ];

        if length(stimContact_oi) > 1
            cont_str = repmat('%s | ',1, length(stimContact_oi));

            cont_str(end-2:end) = '';
        else
            cont_str = '%s';
        end
        tmp_lbl = sprintf([cont_str,'\\newline',...
                             '%s Hz\\newline',...
                             '%s mA\\newline',...
                             '%s \\mus\\newline',...
                             '%s%% DC\\newline'], format_stim_params{:});  

         params_per_se_cell{j} = replace(tmp_lbl, ...
                {' ± 0.00', '0.00m On/',     '0.00m Off', 'NaN ± NaN% duty cycle'},...
                {'',        'continous stim', '',          ''});

     end


    pain_by_params = padcat(pain_params_by_cell{:});

    %%%

    figure('Units', 'Inches', 'Position', [0, 0, 14, 8]);
    sgtitle(...
        sprintf('%s (stage 3 blinded testing)\nlast %g days for given week per stim group',...
            cfg.pt_id, days(cfg.N_days)),...
        'Fontsize',16);
    
    UnivarScatter(pain_by_params,...
            'Compression',50,... distance btwn boxes
            'Width',1, ...
                'MarkerFaceColor','k',...
                'MarkerEdgeColor', 'k',...
                'PointSize', 8,...
                'MeanColor', [.5 .5 .5],... color(s) of the mean line; def. black
                'SEMColor', [.9 .9 .9],...  color(s) of the SEM box; def. dark gray
                'StdColor',[.9 .9 .9]...
                );
        
    if        contains(cfg.plt_metrics{i}, 'VAS');        ylim([0 100]);
    elseif    contains(cfg.plt_metrics{i}, 'NRS');        ylim([0 10]);
    elseif    contains(cfg.plt_metrics{i}, 'MPQ');        ylim([0 45]);  
    end
    
    xticklabels(params_per_se_cell);      set(gca,'Fontsize', 10); grid on;
    ylabel(cfg.plt_metrics{i}, 'FontSize', 24); 
    exportgraphics(gcf, [save_dir, cfg.plt_metrics{i}, '_stage_3_ol_versus_sham (parsimonious grouping)','.png']);

    %% use Wilcoxin signed-rank test to determine sig. stim groups from Sham
    % from PFS, negative control, and "target pain
    % distribution" (i.e., PFS per pain metric - 50% mean of PFS)

    if contains(cfg.plt_metrics{i}, {'MPQ', 'NRS'})
        
      % from parsimonious stim group defintion above
      grp             = stim_grps{:,1};

      STIM_grps       = grp(~contains(grp, 'Sham'));

      i_sham          = strcmp(N_days_back_tbl.R_stim_groups, 'Sham');

      tmp_tbl              = table;
       for i_grp = 1   : length(STIM_grps)

            i_stim          = strcmp(N_days_back_tbl.R_stim_groups, STIM_grps(i_grp));

            [tmp_tbl.p(i_grp), ...
             tmp_tbl.h(i_grp),...
             tmp_tbl.stats{i_grp}] ...
             ...
             = ranksum(...
             ...
             ...
             N_days_back_tbl{i_sham, cfg.plt_metrics{i}}, ...
             N_days_back_tbl{i_stim, cfg.plt_metrics{i}},'method','exact');

       end

       tmp_tbl.Properties.RowNames = STIM_grps;

       sham_vs_stim_stats.(cfg.plt_metrics{i}) = tmp_tbl;
    end

%        tmp_tbl = table;
%        for i_grp = 1   : length(i_grp_oi)
% 
%             [tmp_tbl.p(i_grp), ...
%              tmp_tbl.h(i_grp),...
%              tmp_tbl.stats{i_grp}] ...
%              ...
%              = ranksum(...
%              ...
%              stim_grp{i_grp_oi(i_grp)}, neg_crl_tbl.(cfg.plt_metrics{i}), 'tail','left','method','exact');
% 
%        end
%         sham_vs_stim_stats.neg_crl                      = tmp_tbl;
%         sham_vs_stim_stats.neg_crl.Properties.RowNames  = roi_test;
% 
% 
%        analgesic_thres =  cfg.fluct_study{2}.(cfg.pts{d}){'half_improve', cfg.plt_metrics{i}};
%        fluct_shifted   = fluct_tbl.(cfg.plt_metrics{i})- analgesic_thres;
%        % if any are negative, bring up to floor of 0
%        fluct_shifted(fluct_shifted<0) = 0;
% 
%        tmp_tbl = table;
%        for i_grp = 1   : length(i_grp_oi)
% 
%             [tmp_tbl.p(i_grp), ...
%              tmp_tbl.h(i_grp),...
%              tmp_tbl.stats{i_grp}] ...
%              ...
%              = ranksum(...
%              ...
%              stim_grp{i_grp_oi(i_grp)}, fluct_shifted, 'tail','left','method','exact');
%        end
%        sham_vs_stim_stats.pre_trial_50_percent_mean_shifted = tmp_tbl;
%        sham_vs_stim_stats.pre_trial_50_percent_mean_shifted.Properties.RowNames =roi_test;
% 
% 











end


end