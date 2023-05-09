function stim_sum_cell= plot_freq_amp_pw_cyc(cfg, pt_id, cont)

% bilateral stim--has longer name since both sides are included
if length(cfg.stim_group_label) >20


else  % unilateral stim
    
    side            = cfg.stim_group_label(1);

    compare_vars = cellfun(@(x) [side, '_', x], cfg.seperate_by, 'UniformOutput', false);

    cont.([side, '_', 'cycleOnInSecs'])(cont.([side, '_', 'cycleEnabled']) == 0)  = 0;
    cont.([side, '_', 'cycleOffInSecs'])(cont.([side, '_', 'cycleEnabled']) == 0) = 0;

    cont.([side, '_', 'percentDutyCycle'])...
        (...
        cont.([side, '_', 'cycleEnabled']) == 0 &...
        ~strcmp(cont.([side, '_', 'activeGroup']), 'D')...
        )...
        ...
            = 100;

    cont.([side, '_', 'cycleOnInSecs'])    = cont.([side, '_', 'cycleOnInSecs'])/60;
    cont.([side, '_', 'cycleOffInSecs'])   = cont.([side, '_', 'cycleOffInSecs'])/60;
end
%%
[~, u_params] = unique(cont{:, compare_vars}, 'rows');

%%% for closed-loop stim only break-down by frequency and pulse width
%%% since amplitude and duty cycle are dynamic and better examined in the
%%% INS text logs
if cfg.include_cl && ~strcmp(cfg.cl_cont_group{1}, 'None')
    cl_cont = cfg.cl_cont_group{2};
    
    cl_vars     = cellfun(@(x) [side, '_', x], {'rateInHz', 'pulseWidthInMicroseconds'}, ...
                      'UniformOutput', false);
    
    [~, u_cl_params] = unique(cl_cont{:, cl_vars}, 'rows');

end

%% main plotting loop
for i = 1: length(cfg.plt_metrics)

    params_by_pain_cell = cell(1, length(u_params));

    i_keep                = [];
    ind_per_u_param       = {};
    stim_sum_cell         = {};

    for j = 1 : length(u_params)

         i_params = find(ismember(cont(:, compare_vars), ...
                                  cont(u_params(j), compare_vars), 'rows'));

         tmp_pain_by_params = cont{i_params, cfg.plt_metrics{i}};

        % ensure there are enough non-NaN parameter values to plot
        if sum(~isnan(tmp_pain_by_params)) >= cfg.min_n_reports_subspace
        
            params_by_pain_cell{1,j} = tmp_pain_by_params;

            i_keep                   = [i_keep, j]; %#ok<*AGROW> 
            ind_per_u_param{end+1}   = i_params;        
        end   
    end

%     if length(i_keep) < 2
%         continue % try next pain metric in loop if there's NOT two of them per group
%     end

    u_params        = u_params(i_keep);

    if length(i_keep) ==1
        pain_by_params  = params_by_pain_cell{i_keep};

    elseif isempty(i_keep)
        fprintf('%s | %s has less than %g surveys\n', pt_id, cfg.stim_group_label, cfg.min_n_reports_subspace);
        stim_sum_cell = {table};
        continue

    else
        pain_by_params  = padcat(params_by_pain_cell{i_keep});
    end
    

    if i == 1
        [~, i_sort]  = sort(mean(pain_by_params, 'omitnan'));

        ind_per_u_param = ind_per_u_param(i_sort);

    end

    %%% plot all unique parameters
    if count(cfg.stim_group_label ,'+') == 1    % only unilateral stim for simplicity
        
         stim_sum_tbl = table();

         tick_labels = cell(length(u_params),1);

        for j = 1:length(u_params)

            stim_vars    = cellfun((@(x) [side,'_', x]), ...
                              {'rateInHz','ampInMilliamps', 'pulseWidthInMicroseconds',...
                              'cycleOnInSecs','cycleOffInSecs', 'percentDutyCycle'}, ...
                              ...
                              ...
                              'UniformOutput', false);

            stim_param_oi = cont{ind_per_u_param{j}, stim_vars};

            format_stim_params = cellfun(@(x, y) sprintf('%.2f ± %.2f', x, y),...
                                    num2cell(mean(stim_param_oi)),...
                                    num2cell(std(stim_param_oi)),...
                                    'UniformOutput', false);
        
            stim_sum_tbl(j, stim_vars) = format_stim_params;

            stim_sum_tbl(j, 'ind_u_params') = ind_per_u_param(j);

            if any(strcmp(cfg.seperate_by, 'percentDutyCycle'))
                % string w/o explicit cycling settings
                tmp_str = '%s Hz\\newline%s mA\\newline%s \\mus\\newline%s%% duty cycle';
                i_stim = [1:3,6];
            else
                % string w/ explicit cycling settings
                tmp_str = '%s Hz\\newline%s mA\\newline%s \\mus\\newline%sm On/%sm Off\\newline%s%% duty cycle';
                i_stim  = 1:6;
            end

            tmp_lbl  = sprintf(tmp_str, format_stim_params{i_stim});

            tick_labels{j}  = replace(tmp_lbl, ...
                {' ± 0.00', '0.00m On/',     '0.00m Off', 'NaN ± NaN% duty cycle'},...
                {'',        'continous stim', '',          ''});   
        end
    end
    
    save_dir = [cfg.proc_dir,pt_id,'/',cfg.proc_subdir, '/',cfg.plt_metrics{i},'/'];

    %%% now iterate through closed-loop unique stim parameters
    if cfg.include_cl && ~strcmp(cfg.cl_cont_group{1}, 'None')

        cl_params_by_pain_cell = cell(1, length(u_cl_params));
        
        %%% loop through unique aDBS frequency and pulse-width and save out
        %%% pain surveys and formatted stim parameters labels
        for j_cl = 1 : length(u_cl_params)
           
            i_cl_params = find(...
                            ismember(cl_cont(:, cl_vars), ...
                                     cl_cont(u_cl_params(j_cl), cl_vars), 'rows'));
            
            tmp_pain_by_params               = cl_cont{i_cl_params, cfg.plt_metrics{i}};
            cl_params_by_pain_cell{1,j_cl}   = tmp_pain_by_params;
            
            cl_stim_param_oi             = cl_cont{i_cl_params, stim_vars};

            format_stim_params = cellfun(@(x, y) sprintf('%.2f ± %.2f', x, y),...
                                    num2cell(mean(cl_stim_param_oi)),...
                                    num2cell(std(cl_stim_param_oi)),...
                                    'UniformOutput', false);

            % across unique freq + PW combos include range of amplitudes
            format_stim_params{2}  = sprintf('%.2f–%.2f', ...
                                        min(cl_cont{:, [side, '_ampInMilliamps']}),...
                                        max(cl_cont{:, [side, '_ampInMilliamps']}));

            format_stim_params{6} = 'cl-DBS ';
        
            cl_label              = sprintf(tmp_str, format_stim_params{i_stim});


            tick_labels{j + j_cl} = replace(cl_label, {' ± 0.00'}, {''});  
        end
        
        pain_by_params  = padcat(params_by_pain_cell{i_keep}, cl_params_by_pain_cell{:}); 
        i_plt = [i_sort, (1:length(u_cl_params))+length(i_sort)];

    else
        i_plt = i_sort;

    end
        
    %%% returns specified stim parameter subspace
    stim_sum_cell = {stim_sum_tbl};

    if size(pain_by_params,2) > 10
        figure('Units', 'Inches', 'Position', [0, 0, 24, 6]);

        UnivarScatter(pain_by_params(:, i_plt(1:10)),...
                    'Compression',100,... distance btwn boxes
                    'Width',1,...
                    'MarkerFaceColor',[0 0 0],...
                    'PointSize', 12);

        ylabel(cfg.plt_metrics{i}); grid off
        
        if contains(cfg.plt_metrics{i}, 'VAS')
        ylim([0 100]);
        elseif contains(cfg.plt_metrics{i}, 'NRS')
        ylim([0 10]);
        elseif contains(cfg.plt_metrics{i}, 'MPQ')
        ylim([0 45]);  
        end
        
        title([pt_id, newline, cfg.stim_group_label], 'Fontsize',20);
        
        xticklabels(tick_labels);       set(gca,'Fontsize', 9)
        
        if ~isfolder(save_dir);    mkdir(save_dir);       end
        
        exportgraphics(gcf, [save_dir, cfg.stim_group_label,'_',cfg.plt_metrics{i},' (most_analgesic).png']);

        % repeat for remaining parameters
        figure('Units', 'Inches', 'Position', [0, 0, 24, 6]);
       
        UnivarScatter(pain_by_params(:, i_plt(11:end)),...
                    'Compression',100,... distance btwn boxes
                    'Width',1,...
                    'MarkerFaceColor',[0 0 0],...
                    'PointSize', 12);

        ylabel(cfg.plt_metrics{i}); grid off
        
        if contains(cfg.plt_metrics{i}, 'VAS')
        ylim([0 100]);
        elseif contains(cfg.plt_metrics{i}, 'NRS')
        ylim([0 10]);
        elseif contains(cfg.plt_metrics{i}, 'MPQ')
        ylim([0 45]);  
        end
        
        title([pt_id, newline, cfg.stim_group_label], 'Fontsize',20);

        xticklabels(tick_labels(11:end));       set(gca,'Fontsize', 9)
        
        if ~isfolder(save_dir);         mkdir(save_dir);       end
        
        exportgraphics(gcf, [save_dir, cfg.stim_group_label,'_',cfg.plt_metrics{i}, ' (least_analgesic).png']);

    else

        figure('Units', 'Inches', 'Position', [0, 0, 24, 6]);

        UnivarScatter(pain_by_params(:, i_plt),...
            'Compression',100,... distance btwn boxes
            'Width',1,...
            'MarkerFaceColor',[0 0 0],...
            'PointSize', 12);

        ylabel(cfg.plt_metrics{i}); grid off

        if       contains(cfg.plt_metrics{i}, 'VAS');       ylim([0 100]);
        elseif   contains(cfg.plt_metrics{i}, 'NRS');       ylim([0 10]);
        elseif   contains(cfg.plt_metrics{i}, 'MPQ');       ylim([0 45]);  
        end
    
        title([pt_id, newline, cfg.stim_group_label], 'Fontsize',20);
    
        xticklabels(tick_labels);       set(gca,'Fontsize', 11)
    
        if ~isfolder(save_dir);          mkdir(save_dir);       end
    
        exportgraphics(gcf, [save_dir, cfg.stim_group_label,'_',cfg.plt_metrics{i},'.png']);   
    end
    
    
end

return
%%
% 
% % only look at unique parameters that occur more than N times
% u_params   = u_params(accumarray(i_org_cont,1) >=3);
% 
% i_params              = cell(length(u_params),1);
% nrs_by_params         = [];
% 
% i_keep                = [];
% 
% for j = 1 : length(u_params)
% 
%    i_params{j,1} = find(ismember(cont(:, compare_vars), ...
%              cont(u_params(j),compare_vars), 'rows'));
% 
%    temp_nrs_by_params = {cont.mayoNRS(i_params{j,1})};
% 
%    % ensure there are enough non-NaN parameter values to plot
%    if sum(~isnan(temp_nrs_by_params{1})) > 3
% 
%         nrs_by_params        = [nrs_by_params, {cont.mayoNRS(i_params{j,1})}];
%         i_keep               = [i_keep, j];
% 
%    end      
% end
% 
% u_params   = u_params(i_keep);
% 
% 
% if  isempty(nrs_by_params) 
%     disp([pt_id,...
%         ' |  Not enough stim parameters within ' ...
%                 cont.R_stimContacts{1}]);
%     return
% 
% elseif length(nrs_by_params) > 1
%     nrs_by_params = padcat(nrs_by_params{:});
% 
% else
%     nrs_by_params = nrs_by_params{:};
% end
% 
% [~, i_nrs_sort]         = sort(mean(nrs_by_params,      'omitnan'));
% 
% nrs_by_params           = nrs_by_params(:, i_nrs_sort);
% u_params                = u_params(i_nrs_sort);

% 
% 
% if length(u_params) > 10
% 
%     % for contacts w/ more than 10 unique stim parameters, seperate into
%     % two plots
%     figure('Units', 'Inches', 'Position', [0, 0, 40, 9]);
%     
%     UnivarScatter(nrs_by_params(:, 1:10),...
%         'Compression', 100,... distance btwn boxes
%         'Width', 1,...
%         'MarkerFaceColor',[0 0 0],...
%         'PointSize', 8);
% 
%     title([pt_id, newline, cfg.stim_group_label, ...
%         newline, 'best 10 out of ', num2str(length(u_params)), ' unique parameters'], 'Fontsize',20);
% 
%     i_char_lbl_end = regexp(tick_labels, '\n');
% 
%     xticklabels(tick_labels(1 : i_char_lbl_end(10)));
% 
%     ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
%     set(gca,'Fontsize', 12)
%     
%     saveas(gcf, [by_parm_space_dir, cfg.stim_group_label, '_nrs_overall.png']);
% 
% 
%     % rest of unique stim parameters
%     figure('Units', 'Inches', 'Position', [0, 0, (length(u_params)-10)*4 , 9]);
%     
%     UnivarScatter(nrs_by_params(:, 11:end),...
%         'Compression', 100,... distance btwn boxes
%         'Width', 1,...
%         'MarkerFaceColor',[0 0 0],...
%         'PointSize', 8);
% 
%     title([pt_id, newline, cfg.stim_group_label, ...
%         newline, 'stim parameters NOT in top 10'], 'Fontsize',20);
% 
%     xticklabels(tick_labels(i_char_lbl_end(10)+1 :end));
% 
%     ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
%     set(gca,'Fontsize', 12)
%     
%     saveas(gcf, [by_parm_space_dir, cfg.stim_group_label, '_nrs_overall (least analgesic).png']);
% 
% 
% elseif length(u_params) > 1
% 
%     figure('Units', 'Inches', 'Position', [0, 0, length(u_params)*4 , 9]);
%     
%     UnivarScatter(nrs_by_params,...
%         'Compression', 100,... distance btwn boxes
%         'Width', 1,...
%         'MarkerFaceColor',[0 0 0],...
%         'PointSize', 8);
% 
%     title([pt_id, newline, cfg.stim_group_label], 'Fontsize',20);
% 
%     xticklabels(tick_labels);
% 
%     ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
%     set(gca,'Fontsize', 12)
%     
%     saveas(gcf, [by_parm_space_dir, cfg.stim_group_label, '_nrs_overall.png']);
% 
% 
% % formatting case of single unique parameter differently
% else
% 
%     figure('Units', 'Inches', 'Position', [0, 0, 6 , 9]);
%     UnivarScatter(nrs_by_params,...
%                 'Compression', 5,... distance btwn boxes
%                 'Width', 2,...
%                 'MarkerFaceColor',[0 0 0],...
%                 'PointSize', 8);
% 
%     title([pt_id, newline, cfg.stim_group_label], 'Fontsize',20);
% 
%     xticklabels(tick_labels);
% 
%     ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
%     set(gca,'Fontsize', 12)
%     
%     saveas(gcf, [by_parm_space_dir, cfg.stim_group_label, '_nrs_overall.png']);
% 
% end
end