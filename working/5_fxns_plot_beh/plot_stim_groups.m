function varargout = plot_stim_groups(cfg, pt_id, stimGroups)


%%
set(0,'DefaultFigureVisible','off');            close all

stim_groups = stimGroups.(pt_id);

n_reports   = cellfun(@(x) height(x), table2cell(stim_groups));

i_kept      = n_reports >= cfg.min_n_reports;

n_reports   = n_reports(i_kept);
stim_groups = stim_groups(1, i_kept);

%% visualize pain metrics by contact

% grab pain metrics, and organize as matrix while accounting for different
% N of pain metrics per contact by padding with NaNs
by_cont_dir = [cfg.proc_dir,'/' ,pt_id,'/','by_contacts/'];

if ~isfolder(by_cont_dir);       mkdir(by_cont_dir);      end

%% if desired exlcude surveys w/ VAS 50s and NRS ~= 5
if cfg.exclude_VAS_50s.decision
    VAS_50s_per_group = table;
    VAS_50s_per_group.N       = nan(width(stim_groups),1);
    VAS_50s_per_group.total   = nan(width(stim_groups),1);
    VAS_50s_per_group.percent = nan(width(stim_groups),1);

    for j = 1 : width(stim_groups)
    
         tmp_grp = stim_groups{1, j}{1};

         i_ignore = find(((tmp_grp.painVAS == 50 & tmp_grp.unpleasantVAS == 50)...
                           ));

        stim_groups{1, j}{1}(i_ignore, :) = [];

        VAS_50s_per_group.N(j)     = length(i_ignore);
        VAS_50s_per_group.total(j) = height(tmp_grp); 
        VAS_50s_per_group.percent(j)     = 100*length(i_ignore) / height(tmp_grp);

    end
end
%% main plotting per pain metric (over stim groups) loop

for i = 1: length(cfg.plt_metrics)

    conts_by_pain_lbl = cell(1, width(stim_groups));

    for j = 1 : width(stim_groups)

        tmp_grp = stim_groups{1, j}{1};

        if contains(cfg.plt_metrics{i}, tmp_grp.Properties.VariableNames)
        
            conts_by_pain_lbl{1,j} = stim_groups{1, j}{1}.(cfg.plt_metrics{i});
        else
            conts_by_pain_lbl{1,j} = NaN;
        end

         
    end

    pain_by_cont   = padcat(conts_by_pain_lbl{:});

    if i == 1
        [~, i_sort]  = sort(mean(pain_by_cont, 'omitnan'));
    end

    figure('Units', 'Inches', 'Position', [0, 0, 18 , 8]);
    
    UnivarScatter(pain_by_cont(:, i_sort),...
        'Label', stim_groups.Properties.VariableNames(i_sort),...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[.5 .5 .5],...
        'MarkerEdgeColor', [.5 .5 .5],...
        'PointSize', 7,...
        'MeanColor', 'k',... color(s) of the mean line; def. black
        'SEMColor', [.9 .9 .9],...  color(s) of the SEM box; def. dark gray
        'StdColor',[.9 .9 .9]...
        );
        
    
    ylabel(cfg.plt_metrics{i});

    if       contains(cfg.plt_metrics{i}, 'VAS');         ylim([0 100]);
    elseif   contains(cfg.plt_metrics{i}, 'NRS');         ylim([0 10]);
    elseif   contains(cfg.plt_metrics{i}, 'MPQ');         ylim([0 45]); 
    end

    if cfg.exclude_VAS_50s.decision && cfg.exclude_VAS_50s.plot

        per_grp_txt = compose('%.1f%%\nunanswered\nsurveys', ...
            VAS_50s_per_group.percent);

        title(sprintf(...
            '%s\nUnanswered surveys defined as:\npainVAS = 50 & unplesantVAS = 50', cfg.pt_lbls.(pt_id)), ...
            'Fontsize',14);

        ax = gca;
        x_loc = .75:width(stim_groups)-.25;
        y_loc = repmat(ax.YLim(2)*.95, width(stim_groups), 1);

        text(x_loc,y_loc, per_grp_txt(i_sort), 'FontSize', 7)
    else
        title(cfg.pt_lbls.(pt_id), 'Fontsize',16)
    end
     exportgraphics(gcf, [by_cont_dir, cfg.plt_metrics{i},'.png'])       
end

N_per_contact = cellfun(@(x,y) [y, ' (N = ' , num2str(x), ')'],...
    num2cell(n_reports(i_sort)), stim_groups.Properties.VariableNames(i_sort),...
    'UniformOutput',false)';

disp(N_per_contact);

%% repeat w/ percent clinically relevant, analgesic, and remissive trials



%% for contacts w/ enough reports, break down freq-amp-PW space itself
%{

redcap.R_cycleOnTime(isnan(redcap.R_cycleOnTime))    = 0;

redcap.R_cycleOffTime(isnan(redcap.R_cycleOffTime))  = 0;

redcap.R_i_ol_DBS = contains(redcap.R_activeGroup, {'A', 'B', 'C'});




for i = 1 : size(stimGroups(1,contains(stimGroups.Properties.VariableNames, '+')),2)

    contact_pair = stimGroups.Properties.VariableNames{i};
    if length(contact_pair) > 20 % bilateral stim
      

    elseif strcmp('R', contact_pair(1))
       
        x = stimGroups.(i){1}.R_stimfreq;
        y = stimGroups.(i){1}.R_stimAmp;
        z = stimGroups.(i){1}.R_stimPW;

    
    elseif strcmp('L', contact_pair(1))

        x = stimGroups.(i){1}.L_stimfreq;
        y = stimGroups.(i){1}.L_stimAmp;
        z = stimGroups.(i){1}.L_stimPW;
    

     elseif strcmp('a', contact_pair(1))

        if strcmp('R', contact_pair(6))
            x = stimGroups.(i){1}.R_stimfreq;
            y = stimGroups.(i){1}.R_stimAmp;
            z = stimGroups.(i){1}.R_stimPW;

        elseif strcmp('L', contact_pair(6))
            
            x = stimGroups.(i){1}.L_stimfreq;
            y = stimGroups.(i){1}.L_stimAmp;
            z = stimGroups.(i){1}.L_stimPW;

        end
    end

    if strcmp('R', contact_pair(1)) || strcmp('L', contact_pair(1)) || strcmp('a', contact_pair(1))

        figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
        sgtitle([pt_id, newline, contact_pair, ' parameter space'], 'Fontsize',16)
        

        scatter3(x, y, z,'filled');

        xlabel('freq'); ylabel('amp'); zlabel('pw'); set(gca, 'Fontsize', 14)
        xlim([10, 200]); ylim([0, 8]); zlim([50, 600])

        A = [x,y,z];         [s1,s2,s3] = size(A);

        % Reshape the slices into rows (by working down column-wise, and then across columns)
        A_r = reshape(A,s1*s2, s3,1)';

        % Find the unique rows (which corresponds to unique slices of original array)
        A_ru = unique(A_r,'rows','stable');

        % Reshape the rows back to 3-d
        A_u = reshape(A_ru',s1,s2,[]);

        % find N reports per unique freq-amp-PW AND plot
        dx = 2; dy = 0.05; dz = 3;

        for j = 1 : length(A_u)
            
            A_i = (A_u(j,1) == A(:,1)) & (A_u(j,2) == A(:,2)) & (A_u(j,3) == A(:,3));

            text(A_u(j,1)+dx, A_u(j,2)+dy, A_u(j,3)+dz, ['N = ', num2str(sum(A_i))]);
            hold on

        end
    end
end
%}
%% break down stim_groups (with contacts) further by stim parameters

by_parm_space_dir = [cfg.proc_dir,pt_id,'/',cfg.proc_subdir, '/'];

if ~exist(by_parm_space_dir, 'dir');  mkdir(by_parm_space_dir);     end

stim_group_label       = stim_groups.Properties.VariableNames;
% clDBS_conts            = stim_group_label(contains(stim_group_label, 'clDBS'));
stim_subspace = table();
for i = 1 : size(stim_groups(1,contains(stim_group_label, '+')),2)

    cfg.stim_group_label = stim_group_label{1,i};
     % conceputally unimportant to seperate clDBS by its freq-amp-pw-cyc OR
    % bilateral stim
    if count(cfg.stim_group_label ,'+') == 1 && ~strcmp('cl', cfg.stim_group_label(1:2))

        if cfg.include_cl

            i_cl = find(contains(stim_group_label, 'clDBS ')...
                       & contains(stim_group_label, cfg.stim_group_label));

            if isempty(i_cl)
                 cfg.cl_cont_group = {'None'};
            else
            
                cfg.cl_cont_group = {stim_group_label{i_cl}, stim_groups.(i_cl){1}};
            end
        end

        stim_subspace.(cfg.stim_group_label) = ...
            plot_freq_amp_pw_cyc(cfg, pt_id, stim_groups.(i){1});
    end
end
%%


varargout{1} = stim_groups;
varargout{2} = stim_subspace;

end
