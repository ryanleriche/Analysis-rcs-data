function plot_freq_amp_pw_cyc(cfg, cont)


compared_vars = {'R_stimAmp', 'R_stimfreq', 'R_stimPW'...
                 'R_cycleOnTime', 'R_cycleOffTime'};


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


%% plot all unique parameters
amps         = cont.R_stimAmp(u_params);
freqs        = cont.R_stimfreq(u_params);
pws          = cont.R_stimPW(u_params);

cycleOn      = cont.R_cycleOnTime(u_params)./60;
cycleOff     = cont.R_cycleOffTime(u_params)./60;

label        = [compose('%.1f', freqs'); compose('%.1f', amps'); compose('%.0f', pws');...
               compose('%.1f',  cycleOn'); compose('%.1f',  cycleOff')];

tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s \\mus\\newline%s m On/%s m Off\n', label{:}));
  

if length(u_params) > 1

    figure('Units', 'Inches', 'Position', [0, 0, length(u_params)*4 , 9]);
    
    UnivarScatter(nrs_by_params,...
        'Compression', 100,... distance btwn boxes
        'Width', 1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);

% formatting case of single unique parameter differently
else

    figure('Units', 'Inches', 'Position', [0, 0, 6 , 9]);
    UnivarScatter(nrs_by_params,...
                'Compression', 5,... distance btwn boxes
                'Width', 2,...
                'MarkerFaceColor',[0 0 0],...
                'PointSize', 12);
end

title([cfg.pt_id, newline, cfg.stim_group_label], 'Fontsize',20);


ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
xticklabels(tick_labels);


set(gca,'Fontsize', 12)

by_parm_space_dir = [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id,'/by_freq_amp_pw_cyc/'];

mkdir(by_parm_space_dir)

saveas(gcf, [by_parm_space_dir, cfg.stim_group_label, '_nrs_overall.png']);

end