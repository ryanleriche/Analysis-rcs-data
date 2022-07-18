function    [painscores_out]  = RCS_redcap_painscores(varargin)
%{
  [painscores_out]  = RCS_painScores_REDcap(varargin)

 This will import redcap data for RCS patients using Prasad's
 RedCap API Token for daily surveys

INPUT
  1. (OPTIONAL) PATIENTID
      such as 'RCS01', 'FLUCT', etc.

  If omitted, will get ALL patients RCS01-05 (but not the FLUCTUATION data)

  2. (OPTIONAL) PLOT pain scores ?
        if second value = 1, will plot pain scores, otherwise not



 OUTPUT

  1. PAINSCORES_OUT - a structure, with one field per patient per pain
  scores or streaming notes. Each field contains a table of values


If you want to get the pain fluctuation data, use 'FLUCT' as the
PATIENTID.

EXAMPLE USAGE:

 painscores = RCS_redcap_painscores()

      OR

 painscores = RCS_redcap_painscores('RCS01',1)

      etc...



 ***NOT yet working for FLUCT****

Prasad Shirvalkar MD, PhD
Sept 16, 2021
UCSF
%}
tic


% NOTE THAT For RCS02,04,05 patients, there is a NEW Pain reporting survey
% is concatenated at end

% IMPLANT DATES
% RCS01: 11/19/19
% RCS02L: 9/8/2020 - Explanted 10/13/21
% RCS02R: 9/8/2020 
% RCS04L: 5/13/21
% RCS04R: 5/13/21
% RCS05L: 7/21/21
% RCS05R: 7/21/21

SERVICE            = 'https://redcap.ucsf.edu/api/';


if nargin == 0
    error('Token is required input')

elseif nargin == 1
    TOKEN = varargin{1};
    pt_id_list = ...
        {'RCS01','RCS02','RCS04','RCS05','RCS02new','RCS04new','RCS05new',...
        'RCS01_STREAMING','RCS02_STREAMING','RCS04_STREAMING',...
        'RCS04_STREAMING_v2','RCS05_STREAMING', 'RCS_Weekly', 'RCS_Monthly'};
    

elseif nargin == 2
    TOKEN      = varargin{1};
    pt_id_list = varargin{2};
      

end



for p = 1:numel(pt_id_list)
    
    
    pt_id = pt_id_list{p};
    
    clear redcap*
    
    % Uses REDCap API to fetch redcap pain data based off reportid and
    % patientID
    disp(['Pulling REDcap data for ' pt_id '....'])
    
    

    % Report ID determines which report set to load from. (Daily, Weekly, or Monthly)
    
    switch pt_id
        % old arms
        case 'RCS01'
            PATIENT_ARM = 'rcs01_daily_arm_1';
            reportid = '73191';
        case 'RCS02'
            PATIENT_ARM = 'rcs02_daily_arm_7';
            reportid = '87050';
        case 'RCS04'
            PATIENT_ARM = 'rcs04_daily_arm_10';
            reportid = '104806';
        case 'RCS05'
            PATIENT_ARM = 'rcs05_daily_arm_13';
            reportid = '109667';
            
            % NEW arms
        case  'RCS02new'
            PATIENT_ARM = 'rcs02_new_pain_rep_arm_17';
            reportid = '112131';
        case  'RCS04new'
            PATIENT_ARM = 'rcs04_new_pain_rep_arm_18';
            reportid = '112132';
        case  'RCS05new'
            PATIENT_ARM = 'rcs05_new_pain_rep_arm_19';
            reportid = '112133';
            
            % Streaming arms
        case 'RCS01_STREAMING'
            PATIENT_ARM = 'streaming_arm_9';
            reportid = '95083';
        case 'RCS02_STREAMING'
            PATIENT_ARM = 'streaming_arm_5';
            reportid = '80139';
        case 'RCS04_STREAMING'
            PATIENT_ARM = 'streaming_arm_11';
            reportid = '104807';
        case 'RCS04_STREAMING_v2'          %NOTE That RCS04 has a second streaming arm
            PATIENT_ARM = 'rcs04_streaming_ac_arm_16';
            reportid = '110756';
        case 'RCS05_STREAMING'
            PATIENT_ARM = 'rcs05_recording_se_arm_14';  %NOTE this field is named 'recording' unusually
            reportid = '109668';
            
        case 'FLUCT'
            PATIENT_ARM = 'dbs_and_nondbs_pat_arm_23';
            reportid = '84060';

        case 'RCS_Weekly'
            
            reportid = '135968';

        case 'RCS_Monthly'

            reportid = '135969';
            
            
        otherwise
            
            fprintf('\n Data not found for %s !     ...      Continuing ...  \n\n',pt_id)
            continue
            %
    end
    
    
    disp('************************');
    
    data = webwrite(...
        SERVICE,...
        'token', TOKEN, ...
        'content', 'report',...
        'report_id',reportid, ...
        'format', 'csv',...
        'type','flat',...
        'rawOrLabelHeaders','raw',...
        'exportCheckboxLabel','false',...
        'exportSurveyFields','true',...
        'returnformat','csv');
    
    
    alltable = data;
    
   
    
    if ~contains(pt_id, {'STREAMING', 'Weekly', 'Monthly'}) 
        
        timestampvars = alltable.Properties.VariableNames( find(contains(alltable.Properties.VariableNames,'timestamp')) );
        
        %FOR PAIN SCORE ARMS
         % remove all the extraneous rows ('events' that are different)
        keeprows = strcmp(alltable.redcap_event_name, PATIENT_ARM) & ...
            (arrayfun(@(x) ~isnat(x),alltable.(timestampvars{1})) | ...
            arrayfun(@(x) ~isnat(x),alltable.(timestampvars{2})));
        %for some reason, all patients are using the field rcs01_mpq_...
        
        clntable = alltable(keeprows,:);
        
        % Moves redcap timestamps into designated VAS or MPQ segments. Then combines
        % into consolidated timestamp structure.
        
        redcap_timestamp.vasnrs =  datetime(clntable.(timestampvars{1}));
        redcap_timestamp.mpq = datetime(clntable.(timestampvars{2}));
        time_transfer = find(isnat(redcap_timestamp.vasnrs) & ~isnat(redcap_timestamp.mpq));
        redcap_timestamp.alltimes = redcap_timestamp.vasnrs;
        redcap_timestamp.alltimes(time_transfer) = redcap_timestamp.mpq(time_transfer);
        
        %%%%%%%%% Define Time and Pain Scores %%%%%%%%%%%
        redcap_painscores.time = redcap_timestamp.alltimes;
        % specify timezone for explicit datetime variable
        redcap_painscores.time.TimeZone = '-07:00'; % corresponds to 'America/Los_Angeles' timezone;
        
%         dynamic field names for each subject
        varnames = clntable.Properties.VariableNames;
        nrs_field = varnames{contains(varnames,'intensity_nrs')};
        vas_field = varnames{contains(varnames,'intensity_vas')};
        unp_field = varnames{contains(varnames,'unpleasantness_vas')};
        
        
        % Populate the new flavors of painscores downloaded from redcap
        redcap_painscores.mayoNRS = (clntable.(nrs_field));
        redcap_painscores.painVAS = (clntable.(vas_field));
        redcap_painscores.unpleasantVAS = (clntable.(unp_field));
        
        
        if contains(pt_id,'new')
            worstnrs_field = varnames{contains(varnames,'worst_please_rate')};
            worstvas_field = varnames{contains(varnames,'please_rate_your_pain_inte') & ~contains(varnames,'worst')};
                        
            redcap_painscores.worstNRS = (clntable.(worstnrs_field));
            redcap_painscores.worstVAS = (clntable.(worstvas_field));
            
        else
            redcap_painscores.worstNRS = nan(numel(redcap_painscores.mayoNRS),1);
            redcap_painscores.worstVAS = nan(numel(redcap_painscores.mayoNRS),1);
        end
        
        
        % take the mcgill pain questionnaires
        table_header = clntable.Properties.VariableNames;
        mpq_start = find(contains(table_header,'rcs01_mpq_timestamp'))+1;
        mpq_end = find(contains(table_header,'cruel_punishing'));
        
        mpqhold = table2array(clntable(:,mpq_start:mpq_end));
        
        redcap_painscores.MPQsum = tsnansum((mpqhold),2);
        redcap_painscores.MPQthrobbing = (clntable.throbbing);
        redcap_painscores.MPQshooting = (clntable.shooting);
        redcap_painscores.MPQstabbing = (clntable.stabbing);
        redcap_painscores.MPQsharp = (clntable.sharp);
        redcap_painscores.MPQcramping = (clntable.cramping);
        redcap_painscores.MPQgnawing = (clntable.gnawing);
        redcap_painscores.MPQhot_burning = (clntable.hot_burning);
        redcap_painscores.MPQaching = (clntable.aching);
        redcap_painscores.MPQheavy = (clntable.heavy);
        redcap_painscores.MPQtender = (clntable.tender);
        redcap_painscores.MPQsplitting = (clntable.splitting);
        redcap_painscores.MPQtiring = (clntable.tiring_exhausting);
        redcap_painscores.MPQsickening = (clntable.sickening);
        redcap_painscores.MPQfearful = (clntable.fearful);
        redcap_painscores.MPQcruel = (clntable.cruel_punishing);

        
        % make 0 values NaN for RCS01 because of unreliable reporting
        if strcmp(pt_id,'RCS01')
            redcap_painscores.MPQsum(redcap_painscores.MPQsum == 0) = NaN;
        else
            redcap_painscores.MPQsum = redcap_painscores.MPQsum;
        end
        
        painscores_out.(pt_id) = redcap_painscores;
        

        
    elseif ~contains(pt_id, {'Weekly', 'Monthly'}) 
        %FOR STREAMING NOTES ARMS
        
        timevarNAME = alltable.Properties.VariableNames(    find(contains(alltable.Properties.VariableNames,'timestamp'))     );
        
        keeprows = strcmp(alltable.redcap_event_name, PATIENT_ARM) & ...
            (arrayfun(@(x) ~isnat(x),alltable.(timevarNAME{1})));
              
        clntable = alltable(keeprows,:);
        
        varnames = clntable.Properties.VariableNames;
        medchange_field = varnames{contains(varnames,'medication_changes')};
        activity_field = varnames{contains(varnames,'activity')};
        explain_field = varnames{contains(varnames,'explain')};
        stimon_field = varnames{contains(varnames,'is_stimulation_on')};
        stimprog_field = varnames{contains(varnames,'which_stimulation_program')};

%         **  IN future, will need to add all the individual med fields for
%         each patient **
        
        redcap_streaming.time = clntable.(timevarNAME{1});
        redcap_streaming.medchange = clntable.(medchange_field);
        redcap_streaming.activity = clntable.(activity_field);
        redcap_streaming.explain = clntable.(explain_field);
        redcap_streaming.stimON = clntable.(stimon_field);
        redcap_streaming.stimprog = clntable.(stimprog_field);
        
        painscores_out.(pt_id) = redcap_streaming;
       
        
    end
end

% make into structure of tables
holdfieldnames =  fields(painscores_out);
for n = 1:numel(holdfieldnames)
    painscores_out.(holdfieldnames{n}) = struct2table(painscores_out.(holdfieldnames{n}));
end


oldscores = painscores_out;
clear painscores_out
 
newscores.RCS01 =  oldscores.RCS01;
newscores.RCS02 = [oldscores.RCS02; oldscores.RCS02new];
newscores.RCS04 = [oldscores.RCS04; oldscores.RCS04new];
newscores.RCS05 = [oldscores.RCS05; oldscores.RCS05new];
    
newscores.RCS01_STREAMING = oldscores.RCS01_STREAMING;
newscores.RCS02_STREAMING = oldscores.RCS02_STREAMING;
newscores.RCS04_STREAMING = [oldscores.RCS04_STREAMING; oldscores.RCS04_STREAMING_v2];
newscores.RCS05_STREAMING = oldscores.RCS05_STREAMING;


% OUTPUT
painscores_out = newscores;



toc