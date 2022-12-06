function [redcap, stim_groups] ...
         ...
         =  make_stim_groups(...
         ...
         pt_id, redcap_RCSXXL, redcap_RCSXXR, visits_tbl)
% for troubleshooting fxn as script
%{
i_epoch        = ge(REDcap.RCS04.time, visits.RCS04.dates(11));
pt_id          = 'RCS04';
redcap_RCSXXL  = REDcap.RCS04L(i_epoch,:);

redcap_RCSXXR  = REDcap.RCS04R(i_epoch,:);

visits_tbl     = visits.RCS04;
%}

%%
switch pt_id
    
    case 'RCS02'
        stimRegR      = [{'RACC ', ["0","1","2","3"]}; {'RThal ', ["8","9","10","11"]}];
       

        
    case 'RCS04'
         stimRegL    = [{'LACC ', ["0","1","2","3"]}; {'LCaud ', ["8","9","10","11"]}];
         stimRegR    = [{'RACC ', ["0","1","2","3"]}; {'RThal ', ["8","9","10","11"]}];
        
        
    case 'RCS05'
         stimRegL    = [{'LCaud ', ["0","1","2","3"]}; {'LACC ', ["8","9","10","11"]}];
         stimRegR    = [{'RThal ', ["0","1","2","3"]}; {'RIFG ', ["8","9","10","11"]}];

        
    case 'RCS06'
         stimRegL     = [{'LACC ',  ["0","1","2","3"]}; {'LCaud ', ["8","9","10","11"]}];
         stimRegR     = [{'RThal ', ["0","1","2","3"]}; {'RSFG ', ["8","9","10","11"]}];
        
    case 'RCS07'
        % SGC (not ACC) was bilaterally implanted (RCS data is wrong)
         stimRegL     = [{'LGPi ',  ["0","1","2","3"]}; {'LSGC ', ["8","9","10","11"]}];
         stimRegR     = [{'RThal ', ["0","1","2","3"]}; {'RSGC ', ["8","9","10","11"]}];
        

end

% give explicit right-sided names
i             = find(contains(redcap_RCSXXR.Properties.VariableNames, 'time'));
stim_params   = redcap_RCSXXR.Properties.VariableNames(i(2):end);

R_stim_params = cellfun(@(x) ['R_',x], stim_params, 'UniformOutput', false);

redcap_RCSXXR = renamevars(redcap_RCSXXR, stim_params, R_stim_params);


% left-sided names
if ~strcmp(pt_id, 'RCS02')

    i_stim        = find(strcmp(redcap_RCSXXL.Properties.VariableNames, 'time_stimLog'));
    stim_params   = redcap_RCSXXL.Properties.VariableNames(i_stim:end);
    
    L_stim_params = cellfun(@(x) ['L_',x], stim_params, 'UniformOutput', false);
    
    redcap_RCSXXL = renamevars(redcap_RCSXXL, stim_params, L_stim_params);

end



    % all other pts have bilateral implants--repeat for left side
if ~strcmp(pt_id, 'RCS02')


    % adding in side + region for unambiguous contacts when comparing both sides
    temp         = redcap_RCSXXR;
    
    i_con        = cellfun(@(x) ~isempty(x) , temp.R_stimContacts);
    
    ind_contacts = cellfun(@(x) regexp(x,'\d*','Match'), temp.R_stimContacts(i_con), 'UniformOutput', false);
    ind_contacts = cellfun(@(x) x(1), ind_contacts);
    
    i_small      = cellfun(@(x) any(strcmp(x, stimRegR{1,2})), ind_contacts);
    i_large      = cellfun(@(x) any(strcmp(x, stimRegR{2,2})), ind_contacts);
    
    i_con        = find(i_con);
    
    
    temp.R_stimContacts(i_con(i_small)) =...
        ...
        cellfun(@(x) [stimRegR{1,1}, x], ...
        temp.R_stimContacts(i_con(i_small)), 'UniformOutput', false);
    
    temp.R_stimContacts(i_con(i_large)) =...
        ...
        cellfun(@(x) [stimRegR{2,1}, x], ...
        temp.R_stimContacts(i_con(i_large)), 'UniformOutput', false);
    
     redcap_RCSXXR  = temp;


    temp   = redcap_RCSXXL;
    
    i_con        = cellfun(@(x) ~isempty(x), temp.L_stimContacts);
    
    ind_contacts = cellfun(@(x) regexp(x,'\d*','Match'), temp.L_stimContacts(i_con), 'UniformOutput', false);
    ind_contacts = cellfun(@(x) x(1), ind_contacts);
    
    i_small      = cellfun(@(x) any(strcmp(x, stimRegL{1,2})), ind_contacts);
    i_large      = cellfun(@(x) any(strcmp(x, stimRegL{2,2})), ind_contacts);
    
    i_con        = find(i_con);
    
    temp.L_stimContacts(i_con(i_small)) =...
        ...
        cellfun(@(x) [stimRegL{1,1}, x], ...
        temp.L_stimContacts(i_con(i_small)), 'UniformOutput', false);
    
    temp.L_stimContacts(i_con(i_large)) =...
        ...
        cellfun(@(x) [stimRegL{2,1}, x], ...
        temp.L_stimContacts(i_con(i_large)), 'UniformOutput', false);

    redcap_RCSXXL  = temp;

    % merge left and right sides
    redcap = [redcap_RCSXXR, redcap_RCSXXL(:, i_stim:end)];

    
    % add in visits
    for i = 1 : height(redcap)
    
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
else

    redcap = redcap_RCSXXR;
    redcap = renamevars(redcap, 'R_visits', 'visits');
end
       


%% broadly define stim variants: stim ON, stim OFF, stim @ 0 mA, bilateral stim, etc

% right side

i_R_off            = ~strcmp(redcap.R_therapyStatusDescription, 'On');

i_R_eq0            = redcap.R_ampInMilliamps == 0;
  
i_R_stim_on_ge_0   = ~i_R_off & ~i_R_eq0;

% left side--no pts have unilateral left implant bilateral stim is a possibility
if ~strcmp(pt_id, 'RCS02')

    i_L_eq0           = redcap.L_ampInMilliamps == 0;
    i_L_off           = ~strcmp(redcap.L_therapyStatusDescription, 'On'); 
     
    i_L_stim_on_ge_0  = ~i_L_off & ~i_L_eq0;
    
    i_bilat_on        = i_R_stim_on_ge_0 & i_L_stim_on_ge_0;
    

end

i_s1             = contains(redcap.visits, 's1');
i_s2             = contains(redcap.visits, 's2');
i_s3             = contains(redcap.visits, 's3');

% "fill in" s1 as all reports btwn the last s1 inpatient stay to right
% before first s2 clinic visit

s1               = find(i_s1);
s1_s2_diff       = find(i_s1 - i_s2 == -1);

if ~isempty(s1)
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


else 
    % no reports from S1, means no reports from the first week of S1
    i_s1_first_week = i_s1;
end


s2                      = find(i_s2);
s2_s3_diff              = find(i_s2 - i_s3 == -1); 

if ~isempty(s2)
    if ~isempty(s2_s3_diff)
        i_s2(s2(1): s2_s3_diff(1) -1)                            = 1;
    
        i_s3(s2_s3_diff(1):end)                                  = 1;
    else
        i_s2(s2(1):end) = 1;
    end

else
    i_s2 = 1;
end

% remove if during washout testing or during a clinic/home visit
i_washout         = contains(redcap.visits, 'washout_testing');



i_s1(i_washout | ~strcmp(redcap.visits, {'   '})) = 0;
i_s2(i_washout | ~strcmp(redcap.visits, {'   '})) = 0;
i_s3(i_washout | ~strcmp(redcap.visits, {'   '})) = 0;




if ~strcmp(pt_id, 'RCS02')
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
        disp([pt_id, ': ', num2str(sum(i_B_clDBS)),' reports w/ bilateral stim w/ clDBS'])
    end

    % open-loop unilateral stim
    i_R_ol_on           = ~i_R_off & ~i_R_eq0 & ~i_bilat_on & ~strcmp(redcap.R_activeGroup, 'D')...
                           & (i_L_off | i_L_eq0); 

    i_L_ol_on           = ~i_L_off & ~i_L_eq0 & ~i_bilat_on & ~strcmp(redcap.L_activeGroup, 'D')...
                           & (i_R_off | i_R_eq0); 
    
    
    both_contacts    = cellfun(@(x,y) [x, ' ', y], ...
                       redcap.R_stimContacts, redcap.L_stimContacts, 'UniformOutput', false);
    
    
    clDBS_stimCont    = cellfun(@(x) ['clDBS ' x], [unique(redcap.R_stimContacts(i_clDBS_ROn_Loff));...
                               unique(redcap.L_stimContacts(i_clDBS_LOn_Roff))], 'UniformOutput', false);
    

    stimContacts     = [unique(redcap.R_stimContacts(i_R_ol_on)); ...
                        unique(redcap.L_stimContacts(i_L_ol_on));...
                        unique(both_contacts(i_bilat_on));...
                        clDBS_stimCont];
end



%% seperate stim Groups by contact and On/Off status
stim_groups = table;

if ~strcmp(pt_id, 'RCS02')

    for i = 1 : length(stimContacts)
    
        if strcmp('R', stimContacts{i}(1)) && length(stimContacts{i}) < 20
    
            stim_groups.(stimContacts{i}) = {redcap(strcmp(redcap.R_stimContacts, stimContacts{i}) &... 
                                                                  i_R_ol_on & i_s2, :)};
        
        elseif strcmp('L', stimContacts{i}(1))  && length(stimContacts{i}) < 20
    
            stim_groups.(stimContacts{i}) = {redcap(strcmp(redcap.L_stimContacts, stimContacts{i}) &... 
                                                                  i_L_ol_on & i_s2 , :)};
    
        elseif strcmp('cl', stimContacts{i}(1:2))
    
            if strcmp('R', stimContacts{i}(7))
                stim_groups.(stimContacts{i}) = {redcap(i_clDBS_ROn_Loff & i_s2 , :)};
    
            elseif strcmp('L', stimContacts{i}(7))
                stim_groups.(stimContacts{i}) = {redcap(i_clDBS_LOn_Roff & i_s2, :)};

             end
    
    
        else   % bilateral stim
    
            stim_groups.(stimContacts{i}) = {redcap(strcmp(both_contacts, stimContacts{i}) &... 
                                                                  i_bilat_on & i_s2, :)};
        end
    end
    
    % stimGroups.('s1: both Off')         = {redcap(i_s1 & i_R_off & i_L_off,:)};

    i_s1_both_0mA                       = i_R_eq0 & i_L_eq0 & i_s1 & ~i_R_off & ~i_L_off & ...
                                         ~(i_clDBS_ROn_Loff | i_clDBS_LOn_Roff);

    stim_groups.('s1 both 0mA')         = {redcap(i_s1_both_0mA ,:)};


    i_s1_first_week_both_0mA            = i_R_eq0 & i_L_eq0 & i_s1_first_week & ~i_R_off & ~i_L_off & ...
                                         ~(i_clDBS_ROn_Loff | i_clDBS_LOn_Roff);

    stim_groups.('s1, week 1, both 0 mA')         = {redcap(i_s1_first_week_both_0mA ,:)};
    
    i_s2_both_off                       = i_R_off & i_L_off  & i_s2 &...
                                        ~(i_clDBS_ROn_Loff | i_clDBS_LOn_Roff);

    stim_groups.('s2 both Off')         = {redcap(i_s2_both_off ,:)};


    i_s2_both_0mA                       = i_R_eq0 & i_L_eq0  & i_s2 & ~i_R_off & ~i_L_off & ...
                                        ~(i_clDBS_ROn_Loff | i_clDBS_LOn_Roff);

    stim_groups.('s2 both 0mA')         = {redcap(i_s2_both_0mA,:)};
    stim_groups.('s3')                  = {redcap(i_s3 ,:)};


elseif strcmp(pt_id, 'RCS02')

    
    i_s2_clDBS_ROn      = ~i_R_off & strcmp(redcap.R_activeGroup, 'D') & i_s2;

    clDBS_stimCont   = cellfun(@(x) ['clDBS ' x], ...
                       unique(redcap.R_stimContacts(i_s2_clDBS_ROn)), 'UniformOutput', false)';
    
    stimContacts     = [unique(redcap.R_stimContacts); ...
                        clDBS_stimCont'];

    stimContacts     = stimContacts(~contains(stimContacts, '3+0-2-'));

    for i = 1 : length(stimContacts)

        if ~isempty(stimContacts{i})
    
            if strcmp('R', stimContacts{i}(1))
        
                stim_groups.(stimContacts{i}) = {redcap(strcmp(redcap.R_stimContacts, stimContacts{i}) &... 
                                                                      ~i_R_off & ~i_R_eq0 &...
                                                                      i_s2 & ~strcmp(redcap.R_activeGroup, 'D') , :)};
        
            elseif contains(stimContacts{i}, 'clDBS')
        
                if strcmp('R', stimContacts{i}(7))
        
                    stim_groups.(stimContacts{i}) = {redcap(contains(redcap.R_stimContacts, stimContacts{i}(end-8:end)) ...
                                                    & i_s2_clDBS_ROn, :)};
                end
            end
        end
    end


    stim_groups.('s1, week 1, both 0 mA')       = {redcap(i_s1_first_week & (i_R_off | i_R_eq0) ,:)};
    

    stim_groups.('s1 0mA | Off')                = {redcap(i_s1 & (i_R_off | i_R_eq0) ,:)};


    stim_groups.('s2 Off')              = {redcap(i_R_off  & i_s2, :)};

    stim_groups.('s2 0mA')              = {redcap(i_R_eq0 & i_s2 & ...
                                          ~i_R_off & ~strcmp(redcap.R_activeGroup, 'D') ,:)};

    stim_groups.('washout testing')     = {redcap(i_washout ,:)};

    stim_groups.('s3')                  = {redcap(i_s3 ,:)};


end

    
times       = cellfun(@(x) x.time, table2cell(stim_groups), 'UniformOutput', false);
times       = vertcat(times{:});

[~,ia]      = setdiff(redcap.time, times);

stim_groups.('inclinic, home visits, etc')   = {redcap(ia,:)};

n_reports   = cellfun(@(x) height(x), table2cell(stim_groups));

times       = cellfun(@(x) x.time, table2cell(stim_groups), 'UniformOutput', false);
times       = vertcat(times{:});


disp([pt_id, ' | ', num2str(sum(n_reports)), ' reports assigned to stim_groups']);

disp([pt_id, ' | ', num2str(length(unique(times))), ' unique reports possible']);

end