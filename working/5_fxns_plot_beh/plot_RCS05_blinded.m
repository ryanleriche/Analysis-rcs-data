function plt_struct = plot_RCS05_blinded(cfg, pt_id, stimGroups, INS_logs_API_t_synced)

pt_id = 'RCS05';

blinded_tbl        = stimGroups.(pt_id).blinded_testing{1};

%%

pt_meta = cfg.pt_meta.(pt_id);

blind_start       = pt_meta.dates(strcmp(pt_meta.desc, 'blinded_testing_start'));
blind_stop        = pt_meta.dates(strcmp(pt_meta.desc, 'blinded_testing_stop'));

%% from processed stimLog.json and AppLog.txt, focus on blinded testing indices
L_stim_log        = stimLog.RCS05L;
i_stim_log        = find(ge(L_stim_log.time_stimLog, blind_start) &...
                         le(L_stim_log.time_stimLog, blind_stop));


L_stim_log        = L_stim_log (i_stim_log, :);
%
L_app_log         = INS_logs_API_t_synced.RCS05L.app;
i_app_log         = find(ge(L_app_log.time_INS, blind_start) &...
                         le(L_app_log.time_INS, blind_stop));

L_app_log         = L_app_log(i_app_log, :);
%% 


u_sess = unique(L_app_log.sess_w_same_settings, 'stable');


for i_sett = 1 : length(u_sess)


    tmp_app = L_app_log(...
        L_app_log.sess_w_same_settings == u_sess(i_sett), :);


    tmp_app.prog0mA



end

%{
* w/ left and right stim logs

    * merge into single stim log
        * per L stim log entry
            * add R stim log entry at that time
            * see if any R stim log entries occured btwn given L stim log
            entry and next
                * if so then save those rows to be added below given L stim
                log entry

    * return just subset of blinded testing
    * return stim


%}




%%

% cfg.plt_stim_vars = {'stimContacts', 'rateInHz', 'ampInMilliamps', ...
%                      'pulseWidthInMicroseconds','percentDutyCycle', 'cl_stim'};

cfg.plt_stim_vars = {'ampInMilliamps'};


R_vars           = cellfun(@(x) compose('R_%s',x), cfg.plt_stim_vars);
L_vars           = cellfun(@(x) compose('L_%s',x), cfg.plt_stim_vars);

L_vars = [L_vars, 'sham_stim_washout'];
R_vars = [R_vars, 'sham_stim_washout'];

%%% for plotting purposes, make cl stim an explict string and duty cycle as
%%% 0 
i_R_cl = blinded_tbl.R_cl_stim ==1 ;
i_L_cl = blinded_tbl.L_cl_stim ==1 ;

blinded_tbl.R_cl_stim = [];
blinded_tbl.L_cl_stim = [];

blinded_tbl.R_ol_cl_stim(i_R_cl) = {'closed-loop'};
blinded_tbl.L_ol_cl_stim(i_L_cl) = {'closed-loop'};

blinded_tbl.R_ol_cl_stim(~i_R_cl) = {'open-loop'};
blinded_tbl.L_ol_cl_stim(~i_L_cl) = {'open-loop'};



L_DC_str = cellstr(num2str(blinded_tbl.L_percentDutyCycle));
R_DC_str = cellstr(num2str(blinded_tbl.R_percentDutyCycle));


blinded_tbl.L_percentDutyCycle = [];
blinded_tbl.R_percentDutyCycle = [];

blinded_tbl.L_percentDutyCycle = L_DC_str;
blinded_tbl.R_percentDutyCycle = R_DC_str;



%%

R_off = strcmp(blinded_tbl.R_therapyStatusDescription, 'Off');
L_off = strcmp(blinded_tbl.L_therapyStatusDescription, 'Off');


blinded_tbl.R_off  = R_off;
blinded_tbl.L_off  = L_off;


i_L_0mA = blinded_tbl.L_ampInMilliamps ==0;
i_R_0mA = blinded_tbl.R_ampInMilliamps ==0;


blinded_tbl.sham_stim_washout = repmat({'Check entry'}, height(blinded_tbl), 1);

blinded_tbl.sham_stim_washout(~i_R_cl & ~i_R_0mA)   = {'Right open-loop'};
blinded_tbl.sham_stim_washout(~i_L_cl & ~i_L_0mA)   = {'Left open-loop'};

i_sham =    ((~i_R_cl & i_R_0mA) | R_off) &...
            ((~i_L_cl & i_L_0mA) | L_off) &...
            ~(R_off & L_off);

blinded_tbl.sham_stim_washout(i_sham) = {'Sham (0mA | Off)'};


blinded_tbl.sham_stim_washout(i_R_cl & ~R_off & L_off)     = {'Right closed-loop'};
blinded_tbl.sham_stim_washout(i_L_cl & ~L_off & R_off)     = {'Left closed-loop'};


blinded_tbl.sham_stim_washout(R_off & L_off)     = {'washout (both Off)'};

% blinded_tbl.L_0mA  = L_0mA;


blinded_tbl.stim_groups = findgroups(blinded_tbl.sham_stim_washout);


diff_params   = [1;...
                 diff(blinded_tbl.stim_groups) ~=0 ...
                 ];

i_diffs  = find(diff_params);
i_starts = [];           i_ends   = [];

for h=1:length(i_diffs)-1
    i_starts(h) = i_diffs(h); %#ok<*AGROW> 
    i_ends(h)   = i_diffs(h+1)-1;
end

i_starts = [i_starts, i_ends(end)+1];
i_ends   = [i_ends,height(blinded_tbl)];

for h=1:length(i_starts)
    blinded_tbl.stim_groups_by_changes(i_starts(h):i_ends(h))  = h;
end
%%
% stim_str = ['%s\\newline',...
%              '%g Hz\\newline',...
%              '%g mA\\newline',...
%              '%g \\mus\\newline',...
%              '%s%% DC\\newline',...
%              '%s\\newline',...
%              '%s–%s\\newline'];



stim_str = [...
             '%smA\\newline',...
             '%s\\newline',...
             '%s–%s\n'];

plt_vars ={};
for i = 1:length(cfg.plt_metrics)
    params_by_pain_cell = cell(length(i_starts),1);
    params_per_se_cell  = cell(length(i_starts),1);
    
    for h=1:length(i_starts)
    
        i_params                = find(blinded_tbl.stim_groups_by_changes == h);
        

        R_off_boolean = unique(blinded_tbl{i_params, {'R_off'}});
        L_off_boolean = unique(blinded_tbl{i_params, {'L_off'}});

        t_range = cellstr(...
                    datestr(dateshift(...
                    blinded_tbl.time(i_params([1,end])), 'start', 'day'),...
                    'mmm dd')...
                    );

        % SHAM stimulation
        if R_off_boolean && L_off_boolean

            params_per_se_cell{h} = sprintf('Both Off\\newline%s–%s\n', t_range{:});

        % left side ON
        elseif R_off_boolean && ~L_off_boolean

            plt_tbl  = blinded_tbl(i_params, L_vars);
            tmp_vars = table2cell(unique(plt_tbl, 'rows'));

            i_num = cellfun(@isnumeric, tmp_vars(1,:));
            plt_vars(~i_num) = unique(tmp_vars(:, ~i_num), 'rows');
            plt_vars{i_num}  = sprintf('%.1f–%.1f', min(tmp_vars{:, i_num}), max(tmp_vars{:, i_num}));

            params_per_se_cell{h} = sprintf(stim_str, plt_vars{:}, t_range{:});
        % right side ON
        elseif ~R_off_boolean && L_off_boolean

            plt_tbl  = blinded_tbl(i_params, R_vars);
            tmp_vars = table2cell(unique(plt_tbl, 'rows'));

         
            i_num = cellfun(@isnumeric, tmp_vars(1,:));
            plt_vars(~i_num) = unique(tmp_vars(:, ~i_num), 'rows');
            plt_vars{i_num}  = sprintf('%.1f–%.1f', min(tmp_vars{:, i_num}), max(tmp_vars{:, i_num}));

            params_per_se_cell{h} = sprintf(stim_str, plt_vars{:}, t_range{:});

        else

            params_per_se_cell{h} = sprintf('BOTH sides ON\\newlineunknown mA\n');

        end
      
        params_by_pain_cell{h} =  blinded_tbl.(cfg.plt_metrics{i})(i_params);
    end
    
    pain_by_params               = padcat(params_by_pain_cell{:});

    save_dir = [cfg.proc_dir,pt_id,'/',cfg.proc_subdir,'/',cfg.plt_metrics{i},'/'];
    
    if ~isfolder(save_dir);     mkdir(save_dir); end
    
    figure('Units', 'Inches', 'Position', [0, 0, 10, 8]);
    sgtitle([pt_id, newline, 'blinded testing'], 'Fontsize',16);
    
    UnivarScatter(pain_by_params,...
                'Compression',60,... distance btwn boxes
                'Width',1,...
                'MarkerFaceColor',[.5 .5 .5],...
                'MarkerEdgeColor', [.5 .5 .5],...
                'PointSize', 7,...
                ...
                'PlotStat', 'median with IQR',...
                ...
                'PlotStatColor', 'k',...                   color(s) of the mean or median lines line; def. black
                'DistributionColor_Outer', [.9 .9 .9],...  color(s) of the SEM or (IQR^-n) box
                'DistributionColor_Inner',[.9 .9 .9]...    color(s) of the STD or IQR box
                );
        
        ylabel(cfg.plt_metrics{i});     
        
        if contains(cfg.plt_metrics{i},     'VAS');         ylim([0 100]);
        elseif contains(cfg.plt_metrics{i}, 'NRS');     ylim([0 10]);
        elseif contains(cfg.plt_metrics{i}, 'MPQ');     ylim([0 45]);  
        end
        
        xticklabels([params_per_se_cell{:}]);      set(gca,'Fontsize', 10);
    
    exportgraphics(gcf, [save_dir, cfg.plt_metrics{i}, '_blinded_testing','.png'])    

    %%%
    plt_struct.data = pain_by_params;
    plt_struct.lbls = params_per_se_cell;
end

end