function plot_freq_amp_pw_cyc(cfg, cont)

% bilateral stim--has longer name since both sides are included
if length(cfg.stim_group_label) >20

% unilateral L stim
elseif strcmp(cfg.stim_group_label(1), 'L')

    compared_vars = {'L_ampInMilliamps', 'L_rateInHz', 'L_pulseWidthInMicroseconds'...
                 'L_cycleOnInSecs', 'L_cycleOffInSecs'};


% unilateral R stim
elseif strcmp(cfg.stim_group_label(1), 'R')

    compared_vars = {'R_ampInMilliamps', 'R_rateInHz', 'R_pulseWidthInMicroseconds'...
                 'R_cycleOnInSecs', 'R_cycleOffInSecs'};

end


[~, u_params, i_org_cont] = unique(cont{:, compared_vars}, 'rows');


% only look at unique parameters that occur more than N times
u_params   = u_params(accumarray(i_org_cont,1) >=3);


i_params              = cell(length(u_params),1);
nrs_by_params         = [];

i_keep                 = [];

for j = 1 : length(u_params)

   i_params{j,1} = find(ismember(cont(:, compared_vars), ...
             cont(u_params(j),compared_vars), 'rows'));

   temp_nrs_by_params = {cont.mayoNRS(i_params{j,1})};

   % ensure there are enough non-NaN parameter values to plot
   if sum(~isnan(temp_nrs_by_params{1})) > 3

        nrs_by_params        = [nrs_by_params, {cont.mayoNRS(i_params{j,1})}];
        i_keep               = [i_keep, j];

   end      
end

u_params   = u_params(i_keep);


if  isempty(nrs_by_params) 
    disp([cfg.pt_id,...
        ': Not enough stim parameters wn ' ...
                cont.R_stimContacts{1}]);
    return

elseif length(nrs_by_params) > 1
    nrs_by_params = padcat(nrs_by_params{:});

else
    nrs_by_params = nrs_by_params{:};
end

[~, i_nrs_sort]         = sort(mean(nrs_by_params,      'omitnan'));

nrs_by_params           = nrs_by_params(:, i_nrs_sort);
u_params                = u_params(i_nrs_sort);

%% plot all unique parameters
% bilateral stim--has longer name since both sides are included
if length(cfg.stim_group_label) >20


% unilateral L stim
elseif strcmp(cfg.stim_group_label(1), 'L')

    amps         = cont.L_ampInMilliamps(u_params);
    freqs        = cont.L_rateInHz(u_params);
    pws          = cont.L_pulseWidthInMicroseconds(u_params);
    
    cycleOn      = cont.L_cycleOnInSecs(u_params)./60;
    cycleOff     = cont.L_cycleOffInSecs(u_params)./60;

% unilateral R stim
elseif strcmp(cfg.stim_group_label(1), 'R')

    amps         = cont.R_ampInMilliamps(u_params);
    freqs        = cont.R_rateInHz(u_params);
    pws          = cont.R_pulseWidthInMicroseconds(u_params);
    
    cycleOn      = cont.R_cycleOnInSecs(u_params)./60;
    cycleOff     = cont.R_cycleOffInSecs(u_params)./60;


end



label        = [compose('%.1f', freqs'); compose('%.1f', amps'); compose('%.0f', pws');...
               compose('%.1f',  cycleOn'); compose('%.1f',  cycleOff')];

tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s \\mus\\newline%s m On/%s m Off\n', label{:}));
  

by_parm_space_dir = [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id,'/by_freq_amp_pw_cyc/'];

if ~exist(by_parm_space_dir, 'dir')  
    
    mkdir(by_parm_space_dir)
end

if length(u_params) > 10

    % for contacts w/ more than 10 unique stim parameters, seperate into
    % two plots
    figure('Units', 'Inches', 'Position', [0, 0, 40, 9]);
    
    UnivarScatter(nrs_by_params(:, 1:10),...
        'Compression', 100,... distance btwn boxes
        'Width', 1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 8);

    title([cfg.pt_id, newline, cfg.stim_group_label, ...
        newline, 'best 10 out of ', num2str(length(u_params)), ' unique parameters'], 'Fontsize',20);

    i_char_lbl_end = regexp(tick_labels, '\n');

    xticklabels(tick_labels(1 : i_char_lbl_end(10)));

    ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
    set(gca,'Fontsize', 12)
    
    saveas(gcf, [by_parm_space_dir, cfg.stim_group_label, '_nrs_overall.png']);


    % rest of unique stim parameters
    figure('Units', 'Inches', 'Position', [0, 0, (length(u_params)-10)*4 , 9]);
    
    UnivarScatter(nrs_by_params(:, 11:end),...
        'Compression', 100,... distance btwn boxes
        'Width', 1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 8);

    title([cfg.pt_id, newline, cfg.stim_group_label, ...
        newline, 'stim parameters NOT in top 10'], 'Fontsize',20);

    xticklabels(tick_labels(i_char_lbl_end(10)+1 :end));

    ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
    set(gca,'Fontsize', 12)
    
    saveas(gcf, [by_parm_space_dir, cfg.stim_group_label, '_nrs_overall (least analgesic).png']);


elseif length(u_params) > 1

    figure('Units', 'Inches', 'Position', [0, 0, length(u_params)*4 , 9]);
    
    UnivarScatter(nrs_by_params,...
        'Compression', 100,... distance btwn boxes
        'Width', 1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 8);

    title([cfg.pt_id, newline, cfg.stim_group_label], 'Fontsize',20);

    xticklabels(tick_labels);

    ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
    set(gca,'Fontsize', 12)
    
    saveas(gcf, [by_parm_space_dir, cfg.stim_group_label, '_nrs_overall.png']);


% formatting case of single unique parameter differently
else

    figure('Units', 'Inches', 'Position', [0, 0, 6 , 9]);
    UnivarScatter(nrs_by_params,...
                'Compression', 5,... distance btwn boxes
                'Width', 2,...
                'MarkerFaceColor',[0 0 0],...
                'PointSize', 8);

    title([cfg.pt_id, newline, cfg.stim_group_label], 'Fontsize',20);

    xticklabels(tick_labels);

    ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
    set(gca,'Fontsize', 12)
    
    saveas(gcf, [by_parm_space_dir, cfg.stim_group_label, '_nrs_overall.png']);

end
end