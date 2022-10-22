 function [redcap, stimGroups] ...
     ...
     =  wrt_stim(...
     ...
     cfg, stimLog_w_redcap_RCSXXL, stimLog_w_redcap_RCSXXR, redcap, visits_tbl)
    %{
    cfg             = [];
    cfg.pt_id       = 'RCS04';
    cfg.stimRegL    = [{'LACC ', ["0","1","2","3"]}; {'LCaud ', ["8","9","10","11"]}];
    cfg.stimRegR    = [{'RACC ', ["0","1","2","3"]}; {'RThal ', ["8","9","10","11"]}];
    
    stimLog_w_redcap_RCSXXL = stimLog_w_redcap.RCS04L;
    
    stimLog_w_redcap_RCSXXR = stimLog_w_redcap.RCS04R;
    
    redcap                  = REDcap.RCS04;
    
    visits_tbl                  = visits.RCS04;
    
    n_min_reports           = cfg.min_n_reports_per_contact5;
    %}
    %%
    % cfg             = [];
    % cfg.pt_id       = 'RCS05';
    % cfg.stimRegL    = [{'LCaud ', ["0","1","2","3"]}; {'LACC ', ["8","9","10","11"]}];
    % cfg.stimRegR    = [{'RThal ', ["0","1","2","3"]}; {'RIFG ', ["8","9","10","11"]}];
    % 
    % stimLog_w_redcap_RCSXXL = stimLog_w_redcap.RCS05L;
    % 
    % stimLog_w_redcap_RCSXXR = stimLog_w_redcap.RCS05R;
    % 
    % redcap                  = REDcap.RCS05;
    % 
    % visits_tbl              = visits.RCS05;
    
    % cfg                         = [];
    % cfg.pt_id                   = 'RCS02';
    % cfg.stimRegR                = [{'RACC ', ["0","1","2","3"]}; {'RThal ', ["8","9","10","11"]}];
    % 
    % stimLog_w_redcap_RCSXXL     = [];
    % 
    % stimLog_w_redcap_RCSXXR     = stimLog_w_redcap.RCS02R;
    % 
    % redcap                      = REDcap.RCS02;
    % visits_tbl                  = visits.RCS02;
    
    
    %%
    beh_stim_R      = timetable2table(stimLog_w_redcap_RCSXXR);
    
    % adding in side + region for unambiguous contacts when comparing both sides
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
    
    % all other pts have bilateral implants--repeat for left side
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
    
    redcap.R_activeGroup                = repmat({' '}, height(redcap), 1);
    redcap.R_therapyStatusDescription   = repmat({' '}, height(redcap), 1);
    redcap.R_stimContacts               = repmat({' '}, height(redcap), 1);
    
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
        
        redcap.L_activeGroup               = repmat({' '}, height(redcap), 1);
        redcap.L_therapyStatusDescription  = repmat({' '}, height(redcap), 1);
        redcap.L_stimContacts              = repmat({' '}, height(redcap), 1);
        
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
    
    % repeat for left side
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
        
    
    % see if REDcap report occured during inpatient, inclinic, or home testing
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
    
    
    %% broadly define stim varianst: stim ON, stim OFF, stim @ 0 mA, bilateral stim, etc
    % right side
    
    i_R_off            = ~strcmp(redcap.R_therapyStatusDescription, 'On');
    
    i_R_eq0            = redcap.R_stimAmp == 0;
      
    i_R_stim_on_ge_0   = ~i_R_off & ~i_R_eq0;
    
    % left side--no pts have unilateral left implant bilateral stim is a possibility
    if ~strcmp(cfg.pt_id, 'RCS02')
    
        i_L_eq0           = redcap.L_stimAmp == 0;
        i_L_off           = ~strcmp(redcap.L_therapyStatusDescription, 'On'); 
         
        i_L_stim_on_ge_0  = ~i_L_off & ~i_L_eq0;
        
        i_bilat_on        = i_R_stim_on_ge_0 & i_L_stim_on_ge_0;
        
        i_stim_eq0        = i_R_eq0 & i_L_eq0;
        i_stim_off        = i_R_off & i_L_off;
    
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
    
    % seperate stage 1 from first week after Stage 1 (i.e., post-surgery analgesia)
    
    implant_date            = visits_tbl.dates(strcmp(visits_tbl.desc, 's1_implant'));
    s1_first_week           = implant_date + caldays(7);
    
    
    i_s1_first_week         = ge(redcap.time, implant_date) & le(redcap.time, s1_first_week);
    
    % Stage 1 unless otherwise specified ignores the first week post-surgery
    i_s1(i_s1_first_week)   = 0;
    
    s2                      = find(i_s2);
    s2_s3_diff              = find(i_s2 - i_s3 == -1); 
    
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
        % unilateral cl-DBS includes if the other side
    
        % clDBS unilateral R w/ mA > 0 | mA == 0
        i_clDBS_ROn_Loff    =  strcmp(redcap.R_activeGroup, 'D') & ...
                                ~i_R_off & ...
                               (i_L_off | (i_L_eq0 & ~i_L_off));
        
        % clDBS unilateral L w/ mA > 0 | mA == 0
        i_clDBS_LOn_Roff    = strcmp(redcap.L_activeGroup, 'D') & ...
                                ~i_L_off &...
                               (i_R_off | (i_R_eq0 & ~i_R_off));
        
        
        i_B_clDBS            =  ~i_R_off & ~i_L_off &...
                               (    strcmp(redcap.L_activeGroup, 'D') ...
                                    |...
                                    strcmp(redcap.R_activeGroup, 'D') ...
                               );
    
        if sum(i_B_clDBS) > 0
            disp([cfg.pt_id, ': ', num2str(sum(i_B_clDBS)),' reports w/ bilateral stim w/ clDBS'])
        end
    
        % open-loop unilateral stim
        i_R_on           = ~i_R_off & ~i_R_eq0 & ~i_bilat_on & ~strcmp(redcap.R_activeGroup, 'D')...
                               & (i_L_off | i_L_eq0); 
    
        i_L_on           = ~i_L_off & ~i_L_eq0 & ~i_bilat_on & ~strcmp(redcap.L_activeGroup, 'D')...
                               & (i_R_off | i_R_eq0); 
        
        
        both_contacts    = cellfun(@(x,y) [x, ' ', y], ...
                           redcap.R_stimContacts, redcap.L_stimContacts, 'UniformOutput', false);
        
        
        clDBS_stimCont    = cellfun(@(x) ['clDBS ' x], [unique(redcap.R_stimContacts(i_clDBS_ROn_Loff)),...
                                   unique(redcap.L_stimContacts(i_clDBS_LOn_Roff))], 'UniformOutput', false)';
        
    
        stimContacts     = [unique(redcap.R_stimContacts(i_R_on)); ...
                            unique(redcap.L_stimContacts(i_L_on));...
                            unique(both_contacts(i_bilat_on));...
                            clDBS_stimCont];
    end
    
    
    
    
    %% seperate stim Groups by contact and On/Off status
    stimGroups = table;
    if ~strcmp(cfg.pt_id, 'RCS02')
    
        for i = 1 : length(stimContacts)
        
            if strcmp('R', stimContacts{i}(1)) && length(stimContacts{i}) < 20
        
                stimGroups.(stimContacts{i}) = {redcap(strcmp(redcap.R_stimContacts, stimContacts{i}) &... 
                                                                      i_R_on & i_s2, :)};
            
            elseif strcmp('L', stimContacts{i}(1))  && length(stimContacts{i}) < 20
        
                stimGroups.(stimContacts{i}) = {redcap(strcmp(redcap.L_stimContacts, stimContacts{i}) &... 
                                                                      i_L_on & i_s2 , :)};
        
            elseif strcmp('cl', stimContacts{i}(1:2))
        
                if strcmp('R', stimContacts{i}(7))
                    stimGroups.(stimContacts{i}) = {redcap(i_clDBS_ROn_Loff & i_s2 , :)};
        
                elseif strcmp('L', stimContacts{i}(7))
                    stimGroups.(stimContacts{i}) = {redcap(i_clDBS_LOn_Roff & i_s2, :)};
        
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
                                                                      i_bilat_on & i_s2, :)};
            end
        end
        
        % stimGroups.('s1: both Off')         = {redcap(i_s1 & i_R_off & i_L_off,:)};
    
        i_s1_both_0mA                       = i_R_eq0 & i_L_eq0 & i_s1 & ~i_R_off & ~i_L_off & ...
                                             ~(i_clDBS_ROn_Loff | i_clDBS_LOn_Roff);
    
        stimGroups.('s1 both 0mA')         = {redcap(i_s1_both_0mA ,:)};
    
    
        i_s1_first_week_both_0mA            = i_R_eq0 & i_L_eq0 & i_s1_first_week & ~i_R_off & ~i_L_off & ...
                                             ~(i_clDBS_ROn_Loff | i_clDBS_LOn_Roff);
    
        stimGroups.('s1, week 1, both 0 mA')         = {redcap(i_s1_first_week_both_0mA ,:)};
        
        i_s2_both_off                       = i_R_off & i_L_off  & i_s2 &...
                                            ~(i_clDBS_ROn_Loff | i_clDBS_LOn_Roff);
    
        stimGroups.('s2 both Off')         = {redcap(i_s2_both_off ,:)};
    
    
        i_s2_both_0mA                       = i_R_eq0 & i_L_eq0  & i_s2 & ~i_R_off & ~i_L_off & ...
                                            ~(i_clDBS_ROn_Loff | i_clDBS_LOn_Roff);
    
        stimGroups.('s2 both 0mA')         = {redcap(i_s2_both_0mA,:)};
    
    
    elseif strcmp(cfg.pt_id, 'RCS02')
    
        
        i_s2_aDBS_ROn      = ~i_R_off & strcmp(redcap.R_activeGroup, 'D') & i_s2;
    
        clDBS_stimCont   = cellfun(@(x) ['clDBS ' x], ...
                           unique(redcap.R_stimContacts(i_s2_aDBS_ROn)), 'UniformOutput', false)';
        
        stimContacts     = [unique(redcap.R_stimContacts); ...
                            clDBS_stimCont'];
    
        stimContacts     = stimContacts(~contains(stimContacts, '3+0-2-'));
    
        for i = 1 : length(stimContacts)
    
        if ~isempty(stimContacts{i})
    
            if strcmp('R', stimContacts{i}(1))
        
                stimGroups.(stimContacts{i}) = {redcap(strcmp(redcap.R_stimContacts, stimContacts{i}) &... 
                                                                      ~i_R_off & ~i_R_eq0 &...
                                                                      i_s2 & ~strcmp(redcap.R_activeGroup, 'D') , :)};
        
            elseif contains(stimContacts{i}, 'clDBS')
        
                if strcmp('R', stimContacts{i}(6))
        
                    stimGroups.(stimContacts{i}) = {redcap(contains(redcap.R_stimContacts, stimContacts{i}(end-8:end)) ...
                                                    & i_s2_aDBS_ROn, :)};
                end
            end
        end
        end
    
    
        stimGroups.('s1, week 1, both 0 mA')       = {redcap(i_s1_first_week & (i_R_off | i_R_eq0) ,:)};
        
    
        stimGroups.('s1 0mA | Off')     = {redcap(i_s1 & (i_R_off | i_R_eq0) ,:)};
    
    
        stimGroups.('s2 Off')     = {redcap(i_R_off  & i_s2, :)};
    
        stimGroups.('s2 0mA')     = {redcap(i_R_eq0 & i_s2 & ...
                                      ~i_R_off & ~strcmp(redcap.R_activeGroup, 'D') ,:)};
    
        stimGroups.('washout testing')     = {redcap(i_washout ,:)};
    
    end
    n_reports   = cellfun(@(x) height(x), table2cell(stimGroups));
    
    times       = cellfun(@(x) x.time, table2cell(stimGroups), 'UniformOutput', false);
    times       = vertcat(times{:});
    
    
    [~,ia] = setdiff(redcap.time, times);
    stimGroups.('0mA other Off')   = {redcap(ia,:)};
    
    n_reports   = cellfun(@(x) height(x), table2cell(stimGroups));
    
    times       = cellfun(@(x) x.time, table2cell(stimGroups), 'UniformOutput', false);
    times       = vertcat(times{:});
    
    
    disp([cfg.pt_id, ' N reports: ', num2str(sum(n_reports)),...
        '; N unique reports assigned to stimGroups: ', num2str(length(...
                                             unique(times)))])
    
end