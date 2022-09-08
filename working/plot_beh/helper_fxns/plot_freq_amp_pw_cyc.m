function plot_freq_amp_pw_cyc(cfg, cont)


compared_vars = {'R_stimAmp', 'R_stimfreq', 'R_stimPW'...
                 'R_cycleOnTime', 'R_cycleOffTime'};


[~, u_params, i_org_cont] = unique(cont{:, compared_vars}, 'rows');


% only look at unique parameters that occur more than N times
u_params   = u_params(accumarray(i_org_cont,1) >= 5);


i_params              = cell(length(u_params),1);
nrs_by_params         = [];


for j = 1 : length(u_params)

   i_params{j,1} = find(ismember(cont(:, compared_vars), ...
             cont(u_params(j),compared_vars), 'rows'));

   nrs_by_params        = [nrs_by_params, {cont.mayoNRS(i_params{j,1})}];
        
end



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
               compose('%.0f',  cycleOn'); compose('%.0f',  cycleOff')];

tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s \\mus\\newline%s m On/%s m Off\n', label{:}));


figure('Units', 'Inches', 'Position', [0, 0, length(u_params)*4 , 9]);

title([cfg.pt_id, newline, cfg.stim_group_label], 'Fontsize',20);


UnivarScatter(nrs_by_params,...
    'Compression',100,... distance btwn boxes
    'Width',1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 12);

ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
xticklabels(tick_labels);


set(gca,'Fontsize', 12)

saveas(gcf, [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id,'_',cfg.stim_group_label, '.png']);

end