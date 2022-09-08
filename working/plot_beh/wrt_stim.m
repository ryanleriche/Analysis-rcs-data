% function [redcap, stimGroups, stimGroups] ...
%     ...
%     =  wrt_stim(...
%     ...
%     cfg, stimLog_w_redcap_RCSXXL, stimLog_w_redcap_RCSXXR, redcap, visits_tbl)

% cfg             = [];
% cfg.pt_id       = 'RCS04';
% cfg.stimRegL    = [{'LACC: ', ["0","1","2","3"]}; {'LCaud: ', ["8","9","10","11"]}];
% cfg.stimRegR    = [{'RACC: ', ["0","1","2","3"]}; {'RThal: ', ["8","9","10","11"]}];
% 
% stimLog_w_redcap_RCSXXL = stimLog_w_redcap.RCS04L;
% 
% stimLog_w_redcap_RCSXXR = stimLog_w_redcap.RCS04R;
% 
% redcap                  = REDcap.RCS04;
% 
% visits_tbl                  = visits.RCS04;



%%
% cfg             = [];
% cfg.pt_id       = 'RCS05';
% cfg.stimRegL    = [{'LCaud: ', ["0","1","2","3"]}; {'LACC: ', ["8","9","10","11"]}];
% cfg.stimRegR    = [{'RThal: ', ["0","1","2","3"]}; {'RIFG: ', ["8","9","10","11"]}];
% 
% stimLog_w_redcap_RCSXXL = stimLog_w_redcap.RCS05L;
% 
% stimLog_w_redcap_RCSXXR = stimLog_w_redcap.RCS05R;
% 
% redcap                  = REDcap.RCS05;
% 
% visits_tbl              = visits.RCS05;

cfg                         = [];
cfg.pt_id                   = 'RCS02';
cfg.stimRegR                = [{'RACC ', ["0","1","2","3"]}; {'RThal ', ["8","9","10","11"]}];

stimLog_w_redcap_RCSXXL     = [];

stimLog_w_redcap_RCSXXR     = stimLog_w_redcap.RCS02R;

redcap                      = REDcap.RCS02;
visits_tbl                  = visits.RCS02;


%%
beh_stim_R      = timetable2table(stimLog_w_redcap_RCSXXR);

% adding in side + region for unambiguous contacts when comparing both sides
% repeat for right side
i_con        = cellfun(@(x) ~isempty(x) , beh_stim_R.stimContacts);

ind_contacts = cellfun(@(x) regexp(x,'\d*','Match'), beh_stim_R.stimContacts(i_con), 'UniformOutput', false);

ind_contacts = cellfun(@(x) x(1), ind_contacts);

i_small      = cellfun(@(x) any(strcmp(x, cfg.stimRegR{1,2})), ind_contacts);
i_large      = cellfun(@(x) any(strcmp(x, cfg.stimRegR{2,2})), ind_contacts);

i_con        = find(i_con);


beh_stim_R.stimContacts(i_con(i_small)) =...
    ...
    cellfun(@(x) [cfg.stimRegR{1,1}, x], ...
    beh_stim_R.stimContacts(i_con(i_small)), 'UniformOutput', false);

beh_stim_R.stimContacts(i_con(i_large)) =...
    ...
    cellfun(@(x) [cfg.stimRegR{2,1}, x], ...
    beh_stim_R.stimContacts(i_con(i_large)), 'UniformOutput', false);


if ~strcmp(cfg.pt_id, 'RCS02')

    beh_stim_L   = timetable2table(stimLog_w_redcap_RCSXXL);
    
    i_con        = cellfun(@(x) ~isempty(x), beh_stim_L.stimContacts);
    
    ind_contacts = cellfun(@(x) regexp(x,'\d*','Match'), beh_stim_L.stimContacts(i_con), 'UniformOutput', false);
    ind_contacts = cellfun(@(x) x(1), ind_contacts);
    
    i_small      = cellfun(@(x) any(strcmp(x, cfg.stimRegL{1,2})), ind_contacts);
    i_large      = cellfun(@(x) any(strcmp(x, cfg.stimRegL{2,2})), ind_contacts);
    
    i_con        = find(i_con);
    
    beh_stim_L.stimContacts(i_con(i_small)) =...
        ...
        cellfun(@(x) [cfg.stimRegL{1,1}, x], ...
        beh_stim_L.stimContacts(i_con(i_small)), 'UniformOutput', false);
    
    beh_stim_L.stimContacts(i_con(i_large)) =...
        ...
        cellfun(@(x) [cfg.stimRegL{2,1}, x], ...
        beh_stim_L.stimContacts(i_con(i_large)), 'UniformOutput', false);
end

%% see where StimLogs overlap (i.e., bilateral stim)

% replaces empty cells w/ empty character arrays for ease of handling
i_cell                          = cellfun(@iscell, beh_stim_R.stimContacts);
beh_stim_R.stimContacts(i_cell) = repmat({' '}, sum(i_cell),1);


% redcap.R_time_stimLog = repmat(NaT, height(redcap), 1);

redcap.R_activeGroup  = repmat({' '}, height(redcap), 1);
redcap.R_therapyStatusDescription  = repmat({' '}, height(redcap), 1);
redcap.R_stimContacts  = repmat({' '}, height(redcap), 1);

redcap.R_stimAmp        = zeros(height(redcap), 1);
redcap.R_stimPW         = zeros(height(redcap), 1);
redcap.R_stimfreq       = zeros(height(redcap), 1);

redcap.R_cycleOnTime    = zeros(height(redcap), 1);
redcap.R_cycleOffTime   = zeros(height(redcap), 1);

% repeat for left
% redcap.L_time_stimLog = repmat(NaT, height(redcap), 1);

if ~strcmp(cfg.pt_id, 'RCS02')
    i_cell                          = cellfun(@iscell, beh_stim_L.stimContacts);
    beh_stim_L.stimContacts(i_cell) = repmat({' '}, sum(i_cell),1);
    
    redcap.L_activeGroup  = repmat({' '}, height(redcap), 1);
    redcap.L_therapyStatusDescription  = repmat({' '}, height(redcap), 1);
    redcap.L_stimContacts  = repmat({' '}, height(redcap), 1);
    
    redcap.L_stimAmp  = zeros(height(redcap), 1);
    redcap.L_stimPW   = zeros(height(redcap), 1);
    redcap.L_stimfreq  = zeros(height(redcap), 1);
    
    redcap.L_cycleOnTime    = zeros(height(redcap), 1);
    redcap.L_cycleOffTime   = zeros(height(redcap), 1);

end

% using already found REDcap reports wrt StimLog.jsons, incrp stim
% parameters from both sides w/ REDcap reports
tic
for i = 1 : height(redcap)

    for j = 1 : height(beh_stim_R.i_redcap)
    
        if any(beh_stim_R.i_redcap{j} == i)
        
            redcap.R_activeGroup(i)              = beh_stim_R.activeGroup(j);
            redcap.R_therapyStatusDescription(i) = beh_stim_R.therapyStatusDescription(j);
            redcap.R_stimContacts(i)             = beh_stim_R.stimContacts(j);
            redcap.R_stimfreq(i)                 = beh_stim_R.stimfreq(j);
            
            redcap.R_stimAmp(i)                  = beh_stim_R.stimAmp(j);
            redcap.R_stimPW(i)                   = beh_stim_R.stimPW(j);

            redcap.R_cycleOnTime(i)              = beh_stim_R.cycleOnTime(j);
            redcap.R_cycleOffTime(i)             = beh_stim_R.cycleOffTime(j);



        end
    end

if ~strcmp(cfg.pt_id, 'RCS02')

    for j = 1 : height(beh_stim_L.i_redcap)
    
        if any(beh_stim_L.i_redcap{j} == i)
             
            redcap.L_activeGroup(i)              = beh_stim_L.activeGroup(j);
            redcap.L_therapyStatusDescription(i) = beh_stim_L.therapyStatusDescription(j);
            redcap.L_stimContacts(i)             = beh_stim_L.stimContacts(j);
            redcap.L_stimfreq(i)                 = beh_stim_L.stimfreq(j);
            
            redcap.L_stimAmp(i)                  = beh_stim_L.stimAmp(j);
            redcap.L_stimPW(i)                   = beh_stim_L.stimPW(j);

            redcap.L_cycleOnTime(i)              = beh_stim_L.cycleOnTime(j);
            redcap.L_cycleOffTime(i)             = beh_stim_L.cycleOffTime(j);


        end
    end
end
    
    % see REDcap report occured during inpatient, inclinic, or home testing
    visit_sess_diff      =  visits_tbl.dates - redcap.time(i);
    
    i_visit_day    = find(le(visit_sess_diff, duration('0:00:00')) &...
                  ge(visit_sess_diff, '-24:00:00'));
    
    if ~isempty(i_visit_day)
    
        redcap.visits(i)  = visits_tbl.desc(i_visit_day);
    
    else
    
        redcap.visits(i)  = {'   '};
    
    end


end
toc


%%
% all the "types" of stimulation NOT on

i_R_off           = ~strcmp(redcap.R_therapyStatusDescription, 'On');

i_R_eq0           = redcap.R_stimAmp == 0;
  
i_R_stim_on_ge_0  = ~i_R_off & ~i_R_eq0;


if ~strcmp(cfg.pt_id, 'RCS02')

    i_L_eq0          = redcap.L_stimAmp == 0;
    i_L_off          = ~strcmp(redcap.L_therapyStatusDescription, 'On'); 
     
    i_L_stim_on_ge_0 = ~i_L_off & ~i_L_eq0;
    
    i_bilat_on       = i_R_stim_on_ge_0 & i_L_stim_on_ge_0;
    
    i_stim_eq0       = i_R_eq0 & i_L_eq0;
    i_stim_off       = i_R_off & i_L_off;

end

i_s1             = contains(redcap.visits, 's1');
i_s2             = contains(redcap.visits, 's2');
i_s3             = contains(redcap.visits, 's3');

% "fill in" s1 as all reports btwn the last s1 inpatient stay to right
% before first s2 clinic visit

s1               = find(i_s1);
s1_s2_diff       = find(i_s1 - i_s2 == -1);

if ~isempty(s1_s2_diff)
    i_s1(s1(1): s1_s2_diff(1) - 1)             = 1;

else
    i_s1(s1(1):end) = 1;
end


s2               = find(i_s2);
s2_s3_diff       = find(i_s2 - i_s3 == -1); 

if ~isempty(s2_s3_diff)
    i_s2(s2(1): s2_s3_diff(1) -1)                            = 1;

    i_s3(s2_s3_diff(1):end)                                  = 1;
else
    i_s2(s2(1):end) = 1;
end


% remove if during washout testing or during a clinic/home visit
i_washout         = contains(redcap.visits, 'washout_testing');

i_s2(i_washout | ~strcmp(redcap.visits, {'   '})) = 0;
i_s1(i_washout | ~strcmp(redcap.visits, {'   '})) = 0;


if ~strcmp(cfg.pt_id, 'RCS02')


    % aDBS unilateral R w/ mA > 0 | mA == 0
    i_aDBS_ROn_Loff    =  strcmp(redcap.R_activeGroup, 'D') & ...
                            ~i_R_off & ...
                           (i_L_off | (i_L_eq0 & ~i_L_off));
    
    % aDBS unilateral L w/ mA > 0 | mA == 0
    i_aDBS_LOn_Roff    = strcmp(redcap.L_activeGroup, 'D') & ...
                            ~i_L_off &...
                           (i_R_off | (i_R_eq0 & ~i_R_off));
    
    
    i_B_aDBS            =  ~i_R_off & ~i_L_off &...
                           (    strcmp(redcap.L_activeGroup, 'D') ...
                                |...
                                strcmp(redcap.R_activeGroup, 'D') ...
                           );

    if sum(i_B_aDBS) > 0
        disp([cfg.pt_id, ': ', num2str(sum(i_B_aDBS)),' reports w/ bilateral stim w/ aDBS'])
    end

    % open-loop unilateral stim
    i_R_on           = ~i_R_off & ~i_R_eq0 & ~i_bilat_on & ~i_aDBS_ROn_Loff...
                           & (i_L_off | i_L_eq0) & i_s2; 

    i_L_on           = ~i_L_off & ~i_L_eq0 & ~i_bilat_on & ~i_aDBS_LOn_Roff ...
                           & (i_R_off | i_R_eq0) & i_s2; 
    
    
    both_contacts    = cellfun(@(x,y) [x, ' ', y], ...
                       redcap.R_stimContacts, redcap.L_stimContacts, 'UniformOutput', false);
    
    
    
    aDBS_stimCont    = cellfun(@(x) ['aDBS ' x], [unique(redcap.R_stimContacts(i_aDBS_ROn_Loff)),...
                               unique(redcap.L_stimContacts(i_aDBS_LOn_Roff))], 'UniformOutput', false)';
    
%     aDBS_stimCont    = [aDBS_stimCont ; cellfun(@(x) ['aDBS 0 mA ' x], [unique(redcap.R_stimContacts(i_aDBS_Req0_Loff)),...
%                                unique(redcap.L_stimContacts(i_aDBS_Leq0_Roff))], 'UniformOutput', false)'];

    stimContacts     = [unique(redcap.R_stimContacts(i_R_on)); ...
                        unique(redcap.L_stimContacts(i_L_on));...
                        unique(both_contacts(i_bilat_on));...
                        aDBS_stimCont];
end




%%
stimGroups = table;
if ~strcmp(cfg.pt_id, 'RCS02')

    for i = 1 : length(stimContacts)
    
        if strcmp('R', stimContacts{i}(1)) && length(stimContacts{i}) < 20
    
            stimGroups.(stimContacts{i}) = {redcap(strcmp(redcap.R_stimContacts, stimContacts{i}) &... 
                                                                  i_R_on , :)};
        
        elseif strcmp('L', stimContacts{i}(1))  && length(stimContacts{i}) < 20
    
            stimGroups.(stimContacts{i}) = {redcap(strcmp(redcap.L_stimContacts, stimContacts{i}) &... 
                                                                  i_L_on , :)};
    
        elseif strcmp('a', stimContacts{i}(1))
    
            if strcmp('R', stimContacts{i}(6))
                stimGroups.(stimContacts{i}) = {redcap(i_aDBS_ROn_Loff , :)};
    
            elseif strcmp('L', stimContacts{i}(6))
                stimGroups.(stimContacts{i}) = {redcap(i_aDBS_LOn_Roff , :)};
    
%             % for stimoff | 0 mA stimulation
%             elseif contains(stimContacts{i}, '0 mA')
%     
%                 if strcmp('R', stimContacts{i}(11))
%     
%                     stimGroups.(stimContacts{i}) = {redcap(i_aDBS_Req0_Loff , :)};
%     
%                 elseif strcmp('L', stimContacts{i}(11))
%                     stimGroups.(stimContacts{i}) = {redcap(i_aDBS_Leq0_Roff , :)};
%     
%                 end
             end
    
    
        else   % bilateral stim
    
            stimGroups.(stimContacts{i}) = {redcap(strcmp(both_contacts, stimContacts{i}) &... 
                                                                  i_bilat_on , :)};
        end
    end
    
    % stimGroups.('s1: both Off')         = {redcap(i_s1 & i_R_off & i_L_off,:)};

    i_s1_both_0mA                       = i_R_eq0 & i_L_eq0 & i_s1 & ~i_R_off & ~i_L_off & ...
                                         ~(i_aDBS_ROn_Loff | i_aDBS_LOn_Roff);

    stimGroups.('s1: both 0mA')         = {redcap(i_s1_both_0mA ,:)};
    
    i_s2_both_off                       = i_R_off & i_L_off  & i_s2 &...
                                        ~(i_aDBS_ROn_Loff | i_aDBS_LOn_Roff);

    stimGroups.('s2: both Off')         = {redcap(i_s2_both_off ,:)};

    i_s2_both_0mA                       = i_R_eq0 & i_L_eq0  & i_s2 & ~i_R_off & ~i_L_off & ...
                                        ~(i_aDBS_ROn_Loff | i_aDBS_LOn_Roff);

    stimGroups.('s2: both 0mA')         = {redcap(i_s2_both_0mA,:)};


elseif strcmp(cfg.pt_id, 'RCS02')

    
    i_s2_aDBS_ROn      = ~i_R_off & strcmp(redcap.R_activeGroup, 'D') & i_s2;

    aDBS_stimCont   = cellfun(@(x) ['aDBS ' x], ...
                       unique(redcap.R_stimContacts(i_s2_aDBS_ROn)), 'UniformOutput', false)';
    
    stimContacts     = [unique(redcap.R_stimContacts); ...
                        aDBS_stimCont'];

    stimContacts     = stimContacts(~contains(stimContacts, '3+0-2-'));

    for i = 1 : length(stimContacts)

    if ~isempty(stimContacts{i})

        if strcmp('R', stimContacts{i}(1))
    
            stimGroups.(stimContacts{i}) = {redcap(strcmp(redcap.R_stimContacts, stimContacts{i}) &... 
                                                                  ~i_R_off & ~i_R_eq0 &...
                                                                  i_s2 & ~strcmp(redcap.R_activeGroup, 'D') , :)};
    
        elseif contains(stimContacts{i}, 'aDBS')
    
            if strcmp('R', stimContacts{i}(6))
    
                stimGroups.(stimContacts{i}) = {redcap(contains(redcap.R_stimContacts, stimContacts{i}(end-8:end)) ...
                                                & i_s2_aDBS_ROn, :)};
            end
        end
    end
    end


    stimGroups.('s1 0mA | Off')     = {redcap(i_s1 & (i_R_off | i_R_eq0) ,:)};


    stimGroups.('s2 Off')     = {redcap(i_R_off  & i_s2, :)};

    stimGroups.('s2 0mA')     = {redcap(i_R_eq0 & i_s2 & ...
                                  ~i_R_off & ~strcmp(redcap.R_activeGroup, 'D') ,:)};

    stimGroups.('washout testing')     = {redcap(i_washout ,:)};

end
n_reports   = cellfun(@(x) height(x), table2cell(stimGroups));

times       = cellfun(@(x) x.time, table2cell(stimGroups), 'UniformOutput', false);
times       = vertcat(times{:});


disp([cfg.pt_id, ' N reports: ', num2str(sum(n_reports)),...
    '; N unique reports: ', num2str(length(...
                                         unique(times)))])

times       = cellfun(@(x) x.time, table2cell(stimGroups), 'UniformOutput', false);
times       = vertcat(times{:});


% [~,ia] = setdiff(redcap.time, times);
% stimGroups.('0mA other Off')   = {redcap(ia,:)};

% a = {redcap(ia,:)}

stimGroups = stimGroups(1, (n_reports >=5));
%%
nrs_by_contacts         = [];
vas_by_contacts         = [];
vas_unp_by_contacts     = [];
mpq_tot_by_contacts     = [];

for i = 1 : size(stimGroups,2) 
        nrs_by_contacts         = [nrs_by_contacts, {stimGroups.(i){1}.mayoNRS}];
        
        vas_by_contacts         = [vas_by_contacts; {stimGroups.(i){1}.painVAS}];
        vas_unp_by_contacts     = [vas_unp_by_contacts; {stimGroups.(i){1}.unpleasantVAS}];
        
        mpq_tot_by_contacts     = [mpq_tot_by_contacts; {stimGroups.(i){1}.MPQtotal}];
end

nrs_by_contacts         = padcat(nrs_by_contacts{:});
vas_by_contacts         = padcat(vas_by_contacts{:});
vas_unp_by_contacts     = padcat(vas_unp_by_contacts{:});
mpq_tot_by_contacts     = padcat(mpq_tot_by_contacts{:});


[~, i_vas_sort]         = sort(mean(vas_by_contacts,      'omitnan'));
[~, i_mpq_tot_sort]     = sort(mean(mpq_tot_by_contacts,  'omitnan'));
[~, i_vas_unp_sort]     = sort(mean(vas_unp_by_contacts,  'omitnan'));
[~, i_nrs_sort]         = sort(mean(nrs_by_contacts,      'omitnan'));


%% sorted contacts by lowest to highest pain (left to right)

figure('Units', 'Inches', 'Position', [0, 0, 18 , 8]);
sgtitle([cfg.pt_id], 'Fontsize',16)

UnivarScatter(vas_by_contacts(:,i_vas_sort),...
    'Label', stimGroups.Properties.VariableNames(i_vas_sort),...
    'Compression',100,... distance btwn boxes
    'Width',1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 12);

ylabel({'VAS'});  ylim([0,100]);

saveas(gcf, [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id, '_all_contacts_vas', '.png']);

figure('Units', 'Inches', 'Position', [0, 0, 18 , 8]);
sgtitle([cfg.pt_id], 'Fontsize',16)

UnivarScatter(mpq_tot_by_contacts(:,i_mpq_tot_sort),...
    'Label', stimGroups.Properties.VariableNames(i_mpq_tot_sort),...
    'Compression',100,... distance btwn boxes
    'Width',1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 12);

ylabel({'MPQ Total'}, 'Fontsize', 14);  ylim([0,45]);

saveas(gcf, [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id, '_all_contacts_mpq_total', '.png']);


figure('Units', 'Inches', 'Position', [0, 0, 18 , 8]);
sgtitle([cfg.pt_id], 'Fontsize',16)

UnivarScatter(vas_unp_by_contacts(:,i_vas_unp_sort),...
    'Label', stimGroups.Properties.VariableNames(i_vas_unp_sort),...
    'Compression',100,... distance btwn boxes
    'Width',1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 12);

ylabel({'VAS Unpleasantness'}); ylim([0,100]);

saveas(gcf, [cd,'/plot_beh/figs/stim_groups/', cfg.pt_id, '_all_contacts_vas_unpl', '.png']);

figure('Units', 'Inches', 'Position', [0, 0, 20 , 8]);
sgtitle([cfg.pt_id], 'Fontsize',16)

UnivarScatter(nrs_by_contacts(:,i_nrs_sort),...
    'Label', stimGroups.Properties.VariableNames(i_nrs_sort),...
    'Compression',100,... distance btwn boxes
    'Width',1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 12);

ylabel({'NRS'}, 'Fontsize', 14); ylim([0,10]);

N_per_contact = cellfun(@(x,y) [y, ' (N = ' , num2str(x), ')'],...
    num2cell(n_reports(i_nrs_sort)), stimGroups.Properties.VariableNames(i_nrs_sort),...
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
    cont = stimGroups.("LCaud: c+2-"){1};

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
    cont = stimGroups.("aDBS LCaud: c+2-"){1};

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
    cont = stimGroups.("RThal: c+1-"){1};

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
    cont = stimGroups.("RThal: c+3-"){1};

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