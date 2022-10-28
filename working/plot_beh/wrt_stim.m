% USE THIS TO PLOT PAIN METRICS by stim group
%  Prasad Shirvalkar Oct 25, 2022

% 
%% Set config variables and load databases, RCS logs, redcap etc 
clc
clear
cfg.pt_id                   = 'RCS02';

% Run the config script based on patient 
RCSconfig 

%% Load data 
%  This will load the Json databases, REDcap metrics and textlogs
RCSloaddata

%% RUN THIS IF you want to make the databases/ text logs again 
RCScompiledata 

% need to fix RCS logs to get all data including adaptive 
%  and load old log file and add to it

%% concatenate StimLog.json outputs from RCS database and then MAKE the stim groups
fprintf('Stim log concatenation:\n')
matchStimGroupToRedcap

%% sorted contacts by lowest to highest pain (left to right)

figure('Units', 'Inches', 'Position', [0, 0, 18 , 8]);
sgtitle([cfg.pt_id], 'Fontsize',16)

UnivarScatter(vas_by_contacts(:,sortedmetrics.vas),...
    'Label', stimGroups.Properties.VariableNames(sortedmetrics.vas),...
    'Compression',100,... distance btwn boxes
    'Width',1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 12);

ylabel({'VAS'});  ylim([0,100]);

saveas(gcf, [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id, '_all_contacts_vas', '.png']);

figure('Units', 'Inches', 'Position', [0, 0, 18 , 8]);
sgtitle([cfg.pt_id], 'Fontsize',16)

UnivarScatter(mpq_tot_by_contacts(:,sortedmetrics.mpq),...
    'Label', stimGroups.Properties.VariableNames(sortedmetrics.mpq),...
    'Compression',100,... distance btwn boxes
    'Width',1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 12);

ylabel({'MPQ Total'}, 'Fontsize', 14);  ylim([0,45]);

saveas(gcf, [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id, '_all_contacts_mpq_total', '.png']);


figure('Units', 'Inches', 'Position', [0, 0, 18 , 8]);
sgtitle([cfg.pt_id], 'Fontsize',16)

UnivarScatter(vas_unp_by_contacts(:,sortedmetrics.vasunp),...
    'Label', stimGroups.Properties.VariableNames(sortedmetrics.vasunp),...
    'Compression',100,... distance btwn boxes
    'Width',1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 12);

ylabel({'VAS Unpleasantness'}); ylim([0,100]);

saveas(gcf, [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id, '_all_contacts_vas_unpl', '.png']);

figure('Units', 'Inches', 'Position', [0, 0, 20 , 8]);
sgtitle([cfg.pt_id], 'Fontsize',16)

UnivarScatter(nrs_by_contacts(:,sortedmetrics.nrs),...
    'Label', stimGroups.Properties.VariableNames(sortedmetrics.nrs),...
    'Compression',100,... distance btwn boxes
    'Width',1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 12);

ylabel({'NRS'}, 'Fontsize', 14); ylim([0,10]);

N_per_contact = cellfun(@(x,y) [y, ' (N = ' , num2str(x), ')'],...
    num2cell(sortedmetrics.nreports(sortedmetrics.nrs)), stimGroups.Properties.VariableNames(sortedmetrics.nrs),...
    'UniformOutput',false)';

saveas(gcf, [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id, '_all_contacts_nrs', '.png']);

disp(N_per_contact);

%% for contacts w/ enough reports, break down freq-amp-PW space


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
        sgtitle([cfg.pt_id, newline, contact_pair, ' parameter space'], 'Fontsize',16)
        

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

%% based of N report summary, manually set groupings per pt

if strcmp(cfg.pt_id, 'RCS04')

for i = 1 : size(stimGroups(...
                 1,contains(stimGroups.Properties.VariableNames, '+'))...
                    ,2)

    cfg.stim_group_label = stimGroups.Properties.VariableNames{1,i};

    plot_freq_amp_pw_cyc(cfg, stimGroups.(i){1})

 end
  %{
%% RCS04 LCaud
    % c+9-
    cont = stimGroups.("LCaud: c+9-"){1};

    a_mA  = cont.mayoNRS(cont.L_stimAmp == 1   & cont.L_stimfreq == 125 & cont.L_stimPW== 300);
    b_mA  = cont.mayoNRS(cont.L_stimAmp == 1.5 & cont.L_stimfreq == 125 & cont.L_stimPW== 300);

    freqs        = [125, 125];
    amps         = [  1, 1.5];
    pws          = [300, 300];

    label        = [compose('%.0f', freqs); compose('%.1f', amps); compose('%.0f', pws)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s us\n', label{:}));


    figure('Units', 'Inches', 'Position', [0, 0, 20 , 12]);
    sgtitle(cfg.pt_id, 'Fontsize',16);

    subplot(1,4,1)
    title('LCaud: c+9-','Fontsize',16);
    
    UnivarScatter(padcat(a_mA, b_mA),...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylabel({'NRS'}, 'Fontsize', 14); ylim([0,10]);    grid on;     box off
    xticklabels(tick_labels);

    set(gca,'Fontsize', 14 )

    % c+11-
    cont = stimGroups.("LCaud: c+11-"){1};

    a_mA  = cont.mayoNRS(cont.L_stimAmp == 1   & cont.L_stimfreq == 100 & cont.L_stimPW== 200);
    b_mA  = cont.mayoNRS(cont.L_stimAmp == 1.2 & cont.L_stimfreq == 100 & cont.L_stimPW== 200);

    amps         = [  1, 1.2];
    freqs        = [100, 100];
    pws          = [200, 200];

    label        = [compose('%.0f', freqs); compose('%.1f', amps); compose('%.0f', pws)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s us\n', label{:}));

    subplot(1, 4, 2)
    title('LCaud: c+11-', 'Fontsize',16);
    
    UnivarScatter(padcat(a_mA, b_mA),...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylabel({'NRS'}, 'Fontsize', 18); ylim([0,10]);    grid on;     box off
    
    xticklabels(tick_labels)
    set(gca,'Fontsize', 14 )

    

    % 9+11-
    cont = stimGroups.("LCaud: 9+11-"){1};

    a_freq       = cont.mayoNRS(cont.L_stimfreq == 100 & cont.L_stimAmp == 1 & cont.L_stimPW == 200);

    b_freq_a_mA  = cont.mayoNRS(cont.L_stimfreq == 125 & cont.L_stimAmp == 1   & cont.L_stimPW == 300);
    b_freq_b_mA  = cont.mayoNRS(cont.L_stimfreq == 125 & cont.L_stimAmp == 1.2 & cont.L_stimPW == 300);
    b_freq_c_mA  = cont.mayoNRS(cont.L_stimfreq == 125 & cont.L_stimAmp == 1.5 & cont.L_stimPW == 300);
    b_freq_d_mA  = cont.mayoNRS(cont.L_stimfreq == 125 & cont.L_stimAmp == 2   & cont.L_stimPW == 300);

    freqs        = [100, 125, 125, 125, 125];
    amps         = [  1,   1, 1.2, 1.5,   2];
    pws          = [200, 300, 300, 300, 300];

    label        = [compose('%.0f', freqs); compose('%.1f', amps); compose('%.0f', pws)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s us\n', label{:}));

  
    subplot(1,4,[3,4])
    title('LCaud: 9+11-','Fontsize',16);

    UnivarScatter(padcat(a_freq, b_freq_a_mA, b_freq_b_mA, b_freq_c_mA, b_freq_d_mA),...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylabel({'NRS'}, 'Fontsize', 18); ylim([0,10]);    grid on;     box off

    xticklabels(tick_labels)

    set(gca,'Fontsize', 14 )
 
% LACC: 1+2-; has single stim settings--no need to plot twice (08/29/22)

%{
    cont = freq_amp_PW.("LACC: 1+2-"){1};

    x_mA  = cont.mayoNRS;

    figure('Units', 'Inches', 'Position', [0, 0, 10 , 8]);
    title([cfg.pt_id, newline, 'LACC: 1+2-'], 'Fontsize',16)
    
    UnivarScatter(x_mA,...
        'Compression',100,... distance btwn boxes
        'Width',.25,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylabel({'NRS'}, 'Fontsize', 14); ylim([0,10]);    grid on;     box off

    label        = [compose('%.0f', 100); compose('%.1f', 1.2); compose('%.0f', 200)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s us\n', label{:}));

    xticklabels(tick_labels);
%}

%% RCS04 Right side

    % RThal: 9+11-
    cont = stimGroups.("RThal: 9+11-"){1};

    a_mA  = cont.mayoNRS(cont.R_stimAmp == 0.5 & cont.R_stimfreq == 125 & cont.R_stimPW == 300);
    b_mA  = cont.mayoNRS(cont.R_stimAmp == 1   & cont.R_stimfreq == 125 & cont.R_stimPW == 300);
    c_mA  = cont.mayoNRS(cont.R_stimAmp == 1.2 & cont.R_stimfreq == 125 & cont.R_stimPW == 300);
    d_mA  = cont.mayoNRS(cont.R_stimAmp == 1.5 & cont.R_stimfreq == 125 & cont.R_stimPW == 300);


    amps         = [0.5,   1, 1.2, 1.5];
    freqs        = [125, 125, 125, 125];
    pws          = [300, 300, 300, 300];

    label        = [compose('%.0f', freqs); compose('%.1f', amps); compose('%.0f', pws)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s us\n', label{:}));


    figure('Units', 'Inches', 'Position', [0, 0, 24, 12]);
    sgtitle(cfg.pt_id, 'Fontsize', 16);

    subplot(1, 5, [1,2])
    title('RThal: 9+11-', 'Fontsize',16)
    
    UnivarScatter(padcat(a_mA, b_mA, c_mA, d_mA),...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylabel({'NRS'}, 'Fontsize', 14); ylim([0,10]);    grid on;     box off
 
    xticklabels(tick_labels)

    set(gca,'Fontsize', 14 )

    % RACC: 0+3-
    cont = stimGroups.("RACC: 0+3-"){1};

    % note spacing due to very close frequencies requiring range to capture
    a_freq       = cont.mayoNRS(cont.R_stimAmp == 2   & cont.R_stimfreq == 100                        & cont.R_stimPW == 300);

    b_freq       = cont.mayoNRS(cont.R_stimAmp == 2   & cont.R_stimfreq == 115.7                      & cont.R_stimPW == 250);


    c_freq_a_mA  = cont.mayoNRS(cont.R_stimAmp == 1.5 & cont.R_stimfreq > 149 & cont.R_stimfreq < 151 & cont.R_stimPW >= 200 & cont.R_stimPW <= 250);      

    c_freq_b_mA  = cont.mayoNRS(cont.R_stimAmp == 2   & cont.R_stimfreq > 149 & cont.R_stimfreq < 151 & cont.R_stimPW >= 200 & cont.R_stimPW <= 250);      

    c_freq_c_mA  = cont.mayoNRS(cont.R_stimAmp == 3   & cont.R_stimfreq == 150.6                      & cont.R_stimPW == 200);

    c_freq_d_mA  = cont.mayoNRS(cont.R_stimAmp == 3.5 & cont.R_stimfreq == 150.6                      & cont.R_stimPW == 200);


    amps         = [  2,     2,          1.5,           2,   3,    3.5];
    freqs        = [100, 115.7,          150,         150, 150,    150];
    pws          ={'300', '250', '200 | 250', '200 | 250', '200', '200'};

    label        = [compose('%.1f', freqs); compose('%.1f', amps); pws];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s us\n', label{:}));

    subplot(1,5,[3, 4, 5])
    title('RACC: 0+3-', 'Fontsize',16)
    
    UnivarScatter(padcat(a_freq,...
                         b_freq,...
                         c_freq_a_mA,...
                            c_freq_b_mA,...
                            c_freq_c_mA, ...
                            c_freq_d_mA),...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylabel({'NRS'}, 'Fontsize', 14); ylim([0,10]);    grid on;     box off
 
    xticklabels(tick_labels)

    set(gca,'Fontsize', 14 )


%%
  %}

elseif strcmp(cfg.pt_id, 'RCS05')

%% RCS05 LCaud
    % c+2-
    cont = stimGroups.("LCaud c+2-"){1};

    a_freq_a_mA_a_cyc  = cont.mayoNRS(cont.L_stimAmp == 1     & cont.L_stimfreq == 50  & cont.L_stimPW== 200  ...
                               & cont.L_cycleOnTime== 60 & cont.L_cycleOffTime== 120);

    a_freq_c_mA_a_cyc  = cont.mayoNRS(cont.L_stimAmp == 2     & cont.L_stimfreq == 50  & cont.L_stimPW== 200 ...
                               & cont.L_cycleOnTime== 60 & cont.L_cycleOffTime== 120);


    a_freq_d_mA_a_cyc  = cont.mayoNRS(cont.L_stimAmp == 3     & cont.L_stimfreq == 50  & cont.L_stimPW== 200 ...
                               & cont.L_cycleOnTime== 60 & cont.L_cycleOffTime== 120);


    b_freq_a_cyc        = cont.mayoNRS(cont.L_stimAmp == 3     & cont.L_stimfreq == 100 & cont.L_stimPW== 200 ...
                               & cont.L_cycleOnTime== 60 & cont.L_cycleOffTime == 120);


    amps         = [  1,   2,   3,    3];
    freqs        = [ 50,  50,  50,  100];
    pws          = [200, 200, 200,  200];

    cycleOn      = [  1,   1,   1,   1];
    cycleOff     = [  2,   2,   2,   2];

    label        = [compose('%.0f', freqs); compose('%.1f', amps); compose('%.0f', pws);...
                   compose('%.0f',  cycleOn); compose('%.0f',  cycleOff)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s \\mus\\newline%s m On/%s m Off\n', label{:}));

    figure('Units', 'Inches', 'Position', [0, 0, 28 , 16]);
    sgtitle(cfg.pt_id, 'Fontsize',24);

    subplot(1,10,1:4)
    title('LCaud: c+2- (1 min On, 2 min Off)')
    
    UnivarScatter(padcat(a_freq_a_mA_a_cyc, a_freq_c_mA_a_cyc, a_freq_d_mA_a_cyc, b_freq_a_cyc),...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
    xticklabels(tick_labels);

    set(gca,'Fontsize', 15)

    % repeat w/ longer cycling
    a_freq_a_mA_b_cyc  = cont.mayoNRS(cont.L_stimAmp == 1     & cont.L_stimfreq == 50  & cont.L_stimPW== 200 ...
                               & cont.L_cycleOnTime== 120 & cont.L_cycleOffTime== 1800);

    a_freq_b_mA_b_cyc  = cont.mayoNRS(cont.L_stimAmp == 1.5   & cont.L_stimfreq == 50  & cont.L_stimPW== 200 ...
                               & cont.L_cycleOnTime== 120 & cont.L_cycleOffTime== 1800);


    a_freq_c_mA_b_cyc  = cont.mayoNRS(cont.L_stimAmp == 2     & cont.L_stimfreq == 50  & cont.L_stimPW== 200 ...
                               & cont.L_cycleOnTime== 120 & cont.L_cycleOffTime== 1800);

    a_freq_d_mA_b_cyc  = cont.mayoNRS(cont.L_stimAmp == 3     & cont.L_stimfreq == 50  & cont.L_stimPW== 200 ...
                               & cont.L_cycleOnTime== 120 & cont.L_cycleOffTime== 1800);


    b_freq_b_cyc        = cont.mayoNRS(cont.L_stimAmp == 3     & cont.L_stimfreq == 100 & cont.L_stimPW== 200 ...
                               & cont.L_cycleOnTime== 120 & cont.L_cycleOffTime == 1800);



    amps         = [  1,  1.5,  2,   3,    3];
    freqs        = [ 50,  50,  50,  50,  100];
    pws          = [200, 200, 200, 200,  200];

    cycleOn      = [  2,   2,   2,   2,    2];
    cycleOff     = [ 30,  30,  30,  30,   30];

    label        = [compose('%.0f', freqs); compose('%.1f', amps); compose('%.0f', pws);...
                   compose('%.0f',  cycleOn); compose('%.0f',  cycleOff)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s \\mus\\newline%s m On/%s m Off\n', label{:}));


    subplot(1,10,5:9)
    title('LCaud: c+2- (2 min On, 30 min Off)');
    
    UnivarScatter(padcat(a_freq_a_mA_b_cyc,a_freq_b_mA_b_cyc, a_freq_c_mA_b_cyc, a_freq_d_mA_b_cyc, b_freq_b_cyc),...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylim([0,10]);
    xticklabels(tick_labels);

    set(gca,'Fontsize', 15)



    % aDBS c+2-
    cont = stimGroups.("aDBS LCaud c+2-"){1};

    a_mA  = cont.mayoNRS(cont.L_cycleOnTime== 120 & cont.L_cycleOffTime == 1800);


%     amps         = [];
    freqs        = 130.2;
    pws          = 200;

    cycleOn      = 2;
    cycleOff     = 30;

    label        = [compose('%.0f', freqs); "0-1.5"; compose('%.0f', pws);...
                   compose('%.0f',  cycleOn); compose('%.0f',  cycleOff)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s \\mus\\newline%s m On/%s m Off\n', label{:}));


    subplot(1,10,10)
    title('cl-DBS LCaud: c+2- (2 min On, 30 min)')
    
    UnivarScatter(a_mA,...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylim([0,10]);
    xticklabels(tick_labels);

    set(gca,'Fontsize', 15)

    saveas(gcf, [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id, '_LCaud_c2-_ol-DBS_cl-DBS', '.png']);
%% RCS05 RThal
    %% c+1-
    cont = stimGroups.("RThal c+1-"){1};

    a_mA  = cont.mayoNRS(cont.R_stimAmp == 1     & cont.R_stimfreq == 100  & cont.R_stimPW== 200 ...
                  & cont.R_cycleOnTime== 60 & cont.R_cycleOffTime== 120);

    b_mA  = cont.mayoNRS(cont.R_stimAmp == 2     & cont.R_stimfreq == 100  & cont.R_stimPW== 200 ...
                  & cont.R_cycleOnTime== 60 & cont.R_cycleOffTime== 120);

    c_mA  = cont.mayoNRS(cont.R_stimAmp == 3     & cont.R_stimfreq == 100  & cont.R_stimPW== 200 ...
                  & cont.R_cycleOnTime== 60 & cont.R_cycleOffTime== 120);



    amps         = [   1,  2,   3];
    freqs        = [100, 100, 100];
    pws          = [200, 200, 200];
    cycleOn      = [  1,   1,   1];
    cycleOff     = [  2,   2,   2];

    label        = [compose('%.0f', freqs); compose('%.1f', amps); compose('%.0f', pws);...
                   compose('%.0f',  cycleOn); compose('%.0f',  cycleOff)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s \\mus\\newline%s m On/%s m Off\n', label{:}));

    figure('Units', 'Inches', 'Position', [0, 0, 8 , 16]);
    sgtitle([cfg.pt_id, newline, 'RThal: c+1-'], 'Fontsize',24);


    
    UnivarScatter(padcat(a_mA, b_mA, c_mA),...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylabel({'Numerical Rating Scale (NRS)'}); ylim([0,10]);
    xticklabels(tick_labels);

    set(gca,'Fontsize', 15 )

    saveas(gcf, [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id, '_RThal_c1-ol-DBS', '.png']);
    %% c+3-
    cont = stimGroups.("RThal c+3-"){1};

    a_mA_a_cyc  = cont.mayoNRS(cont.R_stimAmp == 1     & cont.R_stimfreq == 100  & cont.R_stimPW== 200 ...
                  & cont.R_cycleOnTime== 60 & cont.R_cycleOffTime== 120);

    b_mA_a_cyc  = cont.mayoNRS(cont.R_stimAmp == 2     & cont.R_stimfreq == 100  & cont.R_stimPW== 200 ...
                  & cont.R_cycleOnTime== 60 & cont.R_cycleOffTime== 120);

    c_mA_a_cyc  = cont.mayoNRS(cont.R_stimAmp == 3     & cont.R_stimfreq == 100  & cont.R_stimPW== 200 ...
                  & cont.R_cycleOnTime== 60 & cont.R_cycleOffTime== 120);


    amps         = [   1,  2,   3];
    freqs        = [100, 100, 100];
    pws          = [200, 200, 200];

    cycleOn      = [  1,   1,   1];
    cycleOff     = [  2,   2,   2];

    label        = [compose('%.0f', freqs); compose('%.1f', amps); compose('%.0f', pws);...
                   compose('%.0f',  cycleOn); compose('%.0f',  cycleOff)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s \\mus\\newline%s m On/%s m Off\n', label{:}));


    figure('Units', 'Inches', 'Position', [0, 0, 24 , 16]);
    sgtitle(cfg.pt_id, 'Fontsize',24);

    subplot(1,3,1)
    title('RThal: c+3- (1m On/2m Off)');
    
    UnivarScatter(padcat(a_mA_a_cyc, b_mA_a_cyc, c_mA_a_cyc),...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylabel({'Numerical Rating Scale (NRS)'});      ylim([0,10]);    
    xticklabels(tick_labels);

    set(gca,'Fontsize', 15)


    a_mA_b_cyc  = cont.mayoNRS(cont.R_stimAmp == 1     & cont.R_stimfreq == 100  & cont.R_stimPW== 200 ...
                  & cont.R_cycleOnTime== 60 & cont.R_cycleOffTime== 1800);

    b_mA_b_cyc  = cont.mayoNRS(cont.R_stimAmp == 2     & cont.R_stimfreq == 100  & cont.R_stimPW== 200 ...
                  & cont.R_cycleOnTime== 60 & cont.R_cycleOffTime== 1800);

    c_mA_b_cyc  = cont.mayoNRS(cont.R_stimAmp == 3     & cont.R_stimfreq == 100  & cont.R_stimPW== 200 ...
                  & cont.R_cycleOnTime== 60 & cont.R_cycleOffTime== 1800);



    amps         = [   1,  2,   3];
    freqs        = [100, 100, 100];
    pws          = [200, 200, 200];

    cycleOn      = [  1,   1,   1];
    cycleOff     = [  30,  30,  30];

    label        = [compose('%.0f', freqs); compose('%.1f', amps); compose('%.0f', pws);...
                   compose('%.0f',  cycleOn); compose('%.0f',  cycleOff)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s \\mus\\newline%s m On/%s m Off\n', label{:}));


    subplot(1,3,2)
    title('RThal: c+3- (1m On/30m Off)')
    
    UnivarScatter(padcat(a_mA_b_cyc, b_mA_b_cyc, c_mA_b_cyc),...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylim([0,10]); 
    xticklabels(tick_labels);

    set(gca,'Fontsize', 15)

    % aDBS c+3-
    cont = stimGroups.("aDBS RThal: c+3-"){1};

    a_mA  = cont.mayoNRS(cont.R_stimPW== 200 ...
                  & cont.R_cycleOnTime== 60 & cont.R_cycleOffTime== 1800);


    % amps         = [    1,    1.5,  0];

    freqs        = 130.2;
    pws          = 200;

    cycleOn      = 1;
    cycleOff     = 30;

    label        = [compose('%.0f', freqs); "0-1.5"; compose('%.0f', pws);...
                   compose('%.0f',  cycleOn); compose('%.0f',  cycleOff)];
    tick_labels  = strcat(sprintf('%s Hz\\newline%s mA\\newline%s \\mus\\newline%s m On/%s m Off\n', label{:}));



    subplot(1,3,3)
    title('cl-RThal: c+3- (1m On/30m Off)');
    
    UnivarScatter(a_mA,...
        'Compression',100,... distance btwn boxes
        'Width',1,...
        'MarkerFaceColor',[0 0 0],...
        'PointSize', 12);
    
    ylim([0,10]);    
    xticklabels(tick_labels);

    set(gca,'Fontsize', 15)

     saveas(gcf, [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id, '_RThal_c3-_ol-DBS_cl-DBS', '.png']);
%%


elseif strcmp(cfg.pt_id, 'RCS06')

elseif strcmp(cfg.pt_id, 'RCS02')

    % only break down stimGroups w/ contacts further by stim parameters
    for i = 1 : size(stimGroups(...
                     1,contains(stimGroups.Properties.VariableNames, '+'))...
                        ,2)
    
        cfg.stim_group_label = stimGroups.Properties.VariableNames{1,i};
    
        plot_freq_amp_pw_cyc(cfg, stimGroups.(i){1})
    
    end
end

%end

%%
% redcap = groupSum.RCS05.REDcap;
% 
% 
% groupSum.RCS05.stimGroups   = stimGroups;
% groupSum.RCS05.freq_amp_pw  = stimGroups; 
% 
% groupSum.RCS05.REDcap  = redcap;
% 
% % groupSum.RCS04.sII  = sII;
% % groupSum.RCS04.freq_amp_pw = freq_amp_pw; 
%  
% groupSum.RCS04.REDcap  = redcap;



%end