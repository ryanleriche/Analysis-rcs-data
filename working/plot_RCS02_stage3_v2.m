function [N_days_back_tbl, sham_vs_stim_stats] = plot_RCS02_stage3_v2(cfg, s3_tbl, stim_log)

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



t_diff = [diff(s3_stim_log.time_stimLog); duration('')];

s3_stim_log   = addvars(s3_stim_log, t_diff, 'After', 'time_stimLog');


i_le          = le(s3_stim_log.t_diff, duration('0:05:00'));

s3_stim_log(i_le, :)  = [];


break_starts = pt_meta.dates(contains(pt_meta.desc, 'break_start'));
break_stops  = pt_meta.dates(contains(pt_meta.desc, 'break_stop'));


s3_stim_log.break_or_s3   = repmat({'s3'}, height(s3_stim_log),1);

for j = 1:length(break_starts)

    i_break = find(ge(s3_stim_log.time_stimLog, break_starts(j)) &...
                   le(s3_stim_log.time_stimLog, break_stops(j)));

   
    s3_stim_log.break_or_s3(i_break) = {'break'};

end
%%%
plt_vars    = {'stimContacts','ampInMilliamps','rateInHz','pulseWidthInMicroseconds',...
               'percentDutyCycle','therapyStatusDescription'};


s3_stim_log.u_stim_params        = findgroups(s3_stim_log(:, plt_vars));
i_Off                            = find(strcmp(s3_stim_log.therapyStatusDescription, 'Off'));
s3_stim_log.u_stim_params(i_Off) = max(s3_stim_log.u_stim_params) + 1;


i_break                            = find(strcmp(s3_stim_log.break_or_s3, 'break'));
s3_stim_log.u_stim_params(i_break) = max(s3_stim_log.u_stim_params) + 1;

s3_stim_log.u_stim_params          = findgroups(s3_stim_log.u_stim_params);
%%%
diff_params   = [1;...
                 diff(s3_stim_log.u_stim_params) ~=0 ...
                 ];

i_diffs  = find(diff_params);
i_starts = [];           i_ends   = [];

for h=1:length(i_diffs)-1
    i_starts(h) = i_diffs(h); %#ok<*AGROW> 
    i_ends(h)   = i_diffs(h+1)-1;
end

i_starts = [i_starts, i_ends(end)+1];
i_ends   = [i_ends,height(s3_stim_log)];

for h=1:length(i_starts)
    s3_stim_log.stim_params_by_days(i_starts(h):i_ends(h))  = h;
end


for i = 1:length(cfg.plt_metrics)
    params_by_pain_cell = cell(length(i_starts),1);
    params_per_se_cell  = cell(length(i_starts),1);
    
    for h=1:length(i_starts)
    
        i_params                = find(s3_stim_log.stim_params_by_days == h);
        tmp_rcap                = vertcat(s3_stim_log{i_params, "redcap_reports"}{:});

        if isempty(tmp_rcap)
             params_by_pain_cell{h}  = NaN;
        else
            params_by_pain_cell{h}  = tmp_rcap.(cfg.plt_metrics{i});
        end
        

        % return stim parameters over epoch
        tmp_vars = [(table2cell(...
                        unique(s3_stim_log(i_params(1), plt_vars)), 'rows')),...
                     cellstr(...
                        datestr(...
                             dateshift(s3_stim_log.time_stimLog(i_params([1,end])), 'start', 'day'), 'mmm-dd')...
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
                    'Compression',100,... distance btwn boxes
                    'Width',1,...
                    'MarkerFaceColor','k',...
                    'MarkerEdgeColor', 'k',...
                    'PointSize', 8, ...
                    'MeanColor', [.5 .5 .5],... color(s) of the mean line; def. black
                    'SEMColor', [.9 .9 .9],...  color(s) of the SEM box; def. dark gray
                    'StdColor',[.9 .9 .9]...
                );
        
        ylabel(cfg.plt_metrics{i});     
        
        if contains(cfg.plt_metrics{i}, 'VAS');         ylim([0 100]);
        elseif contains(cfg.plt_metrics{i}, 'NRS');     ylim([0 10]);
        elseif contains(cfg.plt_metrics{i}, 'MPQ');     ylim([0 45]);  
        end
        
        xticklabels(params_per_se_cell(i_epochs{h}));      set(gca,'Fontsize', 10);
    end
    
    exportgraphics(gcf, [save_dir, cfg.plt_metrics{i}, '_day-by-day stage 3','.png'])    

    %% now group stim parameters REGARDLESS of when they occured (open-loop stim groups versus sham)
    u_param_days           = unique(s3_stim_log.stim_params_by_days);
    s3_stim_group_cell     = {};
    s3_stim_log.description = cell(height(s3_stim_log), 1);
    for j = 1 : length(u_param_days)
    
        i_params = s3_stim_log.stim_params_by_days == u_param_days(j);
        tmp_tbl  = s3_stim_log(i_params,:);
    
        % ONLY return stim parameters tried for at least 2 days
        if ge(tmp_tbl.time_stimLog(end) -  tmp_tbl.time_stimLog(1), duration('3:00:00:00'))
    
            last_day      = dateshift(tmp_tbl.time_stimLog(end), 'end', 'day');
    
            N_days_before = last_day - cfg.N_days;
    
            s3_stim_group_cell{end+1} = tmp_tbl(ge(tmp_tbl.time_stimLog, N_days_before),:);
    
            s3_stim_log.description(i_params) = {'s3 stim group'};
    
        else
    
            s3_stim_log.description(i_params) = {'washout OR too short'};
    
        end
    end
    
    N_days_back_tbl = vertcat(s3_stim_group_cell{:});
    
    % edge case where before and after vacation was washout so washout appears
    % like many days w/ stim OFF , but few surveys
    i_Off                     = strcmp(N_days_back_tbl.therapyStatusDescription, 'Off');
    N_days_back_tbl(i_Off, :) = [];
    
    
    u_param_lbls = unique(N_days_back_tbl.u_stim_params);
    
    params_by_pain_cell = cell(length(u_param_lbls),1);
    params_per_se_cell  = cell(length(u_param_lbls),1);
    
    
    for h=1:length(u_param_lbls)

        i_params                = find(N_days_back_tbl.u_stim_params == u_param_lbls(h));
        tmp_rcap                = vertcat(N_days_back_tbl{i_params, "redcap_reports"}{:});

        if isempty(tmp_rcap)
             params_by_pain_cell{h}  = NaN;
        else
            params_by_pain_cell{h}  = tmp_rcap.(cfg.plt_metrics{i});
        end

        if strcmp(...
                unique(N_days_back_tbl.break_or_s3(i_params)),...
                'break')

            params_per_se_cell{h} = 'Breaks';
        else
            % return stim parameters over epoch
            tmp_vars = table2cell(...
                            unique(N_days_back_tbl(i_params, plt_vars), 'rows'));
        
            params_per_se_cell{h} = sprintf(['%s\\newline',...
                                             '%g mA\\newline',...
                                             '%g Hz\\newline',...
                                             '%g \\mus\\newline',...
                                             '%g%% DC'], tmp_vars{1:end-1});
        end
    end
    
    pain_by_params               = padcat(params_by_pain_cell{:});
    %%% plot open-loop versus sham stim
    
    figure('Units', 'Inches', 'Position', [0, 0, 12, 10]);
    sgtitle(...
        sprintf('%s (stage 3 blinded testing)\nlast %g days for given week per stim group',...
            cfg.pt_id, days(cfg.N_days)),...
        'Fontsize',20);
    
    UnivarScatter(pain_by_params,...
            'Compression',25,... distance btwn boxes
            'Width',.75, ...
                'MarkerFaceColor','k',...
                'MarkerEdgeColor', 'k',...
                'PointSize', 2, 'PointStyle', 'o',...
                'MeanColor', [.5 .5 .5],... color(s) of the mean line; def. black
                'SEMColor', [.9 .9 .9],...  color(s) of the SEM box; def. dark gray
                'StdColor',[.9 .9 .9]...
                );
        
    if        contains(cfg.plt_metrics{i}, 'VAS');        ylim([0 100]);
    elseif    contains(cfg.plt_metrics{i}, 'NRS');        ylim([0 10]);
    elseif    contains(cfg.plt_metrics{i}, 'MPQ');        ylim([0 45]);  
    end
    
    xticklabels(params_per_se_cell);      set(gca,'Fontsize', 12); grid on;
    ylabel(cfg.plt_metrics{i}, 'FontSize', 24);     

    exportgraphics(gcf, [save_dir, cfg.plt_metrics{i}, '_stage_3_ol_versus_sham','.png']);
%% see N days per group
    if i == 1
        s3_stim_log.R_stim_groups = repmat({'Off'}, height(s3_stim_log), 1);
    
        i_on = strcmp(s3_stim_log.therapyStatusDescription, 'On');
         %%%
        s3_stim_log.R_stim_groups(i_on & ...
            strcmp(s3_stim_log.stimContacts, 'RACC 1+0-'))    = {'RACC 1+0-'};
        
        s3_stim_log.R_stim_groups(i_on & ...
            strcmp(s3_stim_log.stimContacts, 'RACC 3+2-'))    = {'RACC 3+2-'};
        
        s3_stim_log.R_stim_groups(i_on & ...
            s3_stim_log.rateInHz == 2)    = {'Sham'};
    
    
        i_break                     = strcmp(s3_stim_log.break_or_s3, 'break');
        s3_stim_log.R_stim_groups(i_break)     = {'break'};
    
        
        s3_stim_log.u_stim_grps          = findgroups(s3_stim_log.R_stim_groups);
        %%%
        diff_params   = [1;...
                         diff(s3_stim_log.u_stim_grps ) ~=0 ...
                         ];
        
        i_diffs  = find(diff_params);
        i_starts = [];           i_ends   = [];
        
        for h=1:length(i_diffs)-1
            i_starts(h) = i_diffs(h); %#ok<*AGROW> 
            i_ends(h)   = i_diffs(h+1)-1;
        end
        
        i_starts = [i_starts, i_ends(end)+1];
        i_ends   = [i_ends,height(s3_stim_log)];
    
    
        
        for h=1:length(i_starts)
            s3_stim_log.stim_grps_by_days(i_starts(h):i_ends(h))  = h;
        end
    
    
      %%%
        grp_dur_tbl = table;
    
        grp_dur_tbl.start     = s3_stim_log.time_stimLog(i_starts);
        grp_dur_tbl.end       = s3_stim_log.time_stimLog(i_ends);
        grp_dur_tbl.duration  = s3_stim_log.time_stimLog(i_ends) - ...
                                     s3_stim_log.time_stimLog(i_starts);
    
        grp_dur_tbl.grp_lbl   = s3_stim_log.R_stim_groups(i_starts);
    
        grp_lbls = unique(grp_dur_tbl.grp_lbl);
    
        grp_lbls = grp_lbls(~strcmp(grp_lbls, 'break'));
    
        grp_durs      = cell(length(grp_lbls), 1 );
        grp_prop_fin  = nan(length(grp_lbls), 1 );
    
        for j_grp = 1 : length(grp_lbls)
    
           j_grp_oi           = find(strcmp(grp_dur_tbl.grp_lbl, grp_lbls(j_grp)));
    
           N_days_vec         = days(grp_dur_tbl.duration(j_grp_oi));
           
           N_days_vec        = N_days_vec(N_days_vec >.2);
           
           grp_prop_fin(j_grp)       = 100 * sum(N_days_vec > 4.8) / length(N_days_vec);
           grp_durs{j_grp}           =  N_days_vec ;
    
        end
    
        grp_durs_mat = padcat(grp_durs{:});
    
    
    
        figure('Units', 'Inches', 'Position', [0, 0, 8, 10]);
        tiledlayout(2,1, 'TileSpacing', 'Compact');
        nexttile
    
        UnivarScatter(grp_durs_mat,...
                'Compression',25,... distance btwn boxes
                'Width',.75, ...
                    'MarkerFaceColor','k',...
                    'MarkerEdgeColor', 'k',...
                    'PointSize', 2, 'PointStyle', 'o',...
                    'MeanColor', [.5 .5 .5],... color(s) of the mean line; def. black
                    'SEMColor', [.9 .9 .9],...  color(s) of the SEM box; def. dark gray
                    'StdColor',[.9 .9 .9]...
                    );
        xticklabels([])
        ylabel('N days')
    
        set(gca, 'FontSize', 14, 'TickLength', [0,0])
    
        nexttile
            bar(grp_prop_fin, 'FaceColor', [0.5, 0.5, 0.5]);
            
            xlim([.5, length(grp_prop_fin)+.5])
            ylim([0,100]);    ylabel('Percent of weeks w/o early washout')
            
            
            xticklabels(grp_lbls); set(gca, 'FontSize', 14, 'TickLength', [0,0])
            
        sgtitle(...
        sprintf('%s (stage 3 blinded testing)\nN days before switching stim groups',...
        cfg.pt_id), 'Fontsize',14); hold on
    
         exportgraphics(gcf, [cfg.proc_dir,cfg.pt_id,'/',cfg.proc_subdir,'/stage_3_time_per_group.png']);
    end
%%  NOW group into more parsimonious SHAM vs stim groups

    %%% remove periods of break
    i_break                     = strcmp(N_days_back_tbl.break_or_s3, 'break');
    N_days_back_tbl(i_break, :) = [];


    numeric_vars = {'rateInHz','ampInMilliamps','pulseWidthInMicroseconds',...
                    'percentDutyCycle'};

    %%% identify RACC 1+0- stim explicity, THEN write over w/ Sham
    %%% tag
     N_days_back_tbl.R_stim_groups(...
         strcmp(N_days_back_tbl.stimContacts, 'RACC 1+0-'))    = {'RACC 1+0-'};

     N_days_back_tbl.R_stim_groups(...
         strcmp(N_days_back_tbl.stimContacts, 'RACC 3+2-'))    = {'RACC 3+2-'};


     N_days_back_tbl.R_stim_groups(N_days_back_tbl.rateInHz == 2)    = {'Sham'};
    

    [stim_grps, u_params]   = unique(N_days_back_tbl(:, cfg.seperate_by), 'rows');

    params_by_pain_cell     = cell(length(u_params), 1);
    params_per_se_cell      = params_by_pain_cell ;

     for h = 1 : length(u_params)

     
         i_params = find(...
                            ismember(N_days_back_tbl(:, cfg.seperate_by),...
                                     N_days_back_tbl(u_params(h), cfg.seperate_by))...
                                     );

        tmp_rcap                = vertcat(N_days_back_tbl{i_params, "redcap_reports"}{:});

        if isempty(tmp_rcap)
             params_by_pain_cell{h}  = NaN;
        else
            params_by_pain_cell{h}  = tmp_rcap.(cfg.plt_metrics{i});
        end

        % return stim parameters
        tmp_vars       = N_days_back_tbl{i_params, numeric_vars};
        stimContact_oi = unique(N_days_back_tbl{i_params, 'stimContacts'})';


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

         params_per_se_cell{h} = replace(tmp_lbl, ...
                {' ± 0.00', '0.00m On/',     '0.00m Off', 'NaN ± NaN% duty cycle'},...
                {'',        'continous stim', '',          ''});

     end


    pain_by_params = padcat(params_by_pain_cell{:});

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


            sham_tbl        = vertcat(N_days_back_tbl{i_sham, "redcap_reports"}{:});
            stim_tbl        = vertcat(N_days_back_tbl{i_stim, "redcap_reports"}{:});


            [tmp_tbl.p(i_grp), ...
             tmp_tbl.h(i_grp),...
             tmp_tbl.stats{i_grp}] ...
             ...
             = ranksum(...
             ...
             ...
            sham_tbl.(cfg.plt_metrics{i}), stim_tbl.(cfg.plt_metrics{i}), 'method','exact');

       end

       tmp_tbl.Properties.RowNames = STIM_grps;

       sham_vs_stim_stats.(cfg.plt_metrics{i}) = tmp_tbl;
    end
end
end