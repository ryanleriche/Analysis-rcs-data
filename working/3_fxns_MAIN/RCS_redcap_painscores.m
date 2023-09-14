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
    rcs_TOKEN = varargin{1};
    pt_id_list = ...
        {'RCS02','RCS04','RCS05','RCS02new','RCS04new','RCS05new','RCS06', 'RCS07',...
        ...
        'RCS02_STREAMING','RCS04_STREAMING',...
        'RCS04_STREAMING_v2','RCS05_STREAMING', 'RCS06_STREAMING', 'RCS07_STREAMING'...
        ...
        'RCS_Weekly', 'RCS_Monthly',...
        'FlUCT'};
    

elseif nargin == 2
    rcs_TOKEN      = varargin{1};
    pcs_TOKEN      = varargin{2};
    pt_id_list = ...
        {'RCS01','RCS02','RCS04','RCS05','RCS02new','RCS04new','RCS05new','RCS06',...
        ...
        'RCS01_STREAMING','RCS02_STREAMING','RCS04_STREAMING',...
        'RCS04_STREAMING_v2','RCS05_STREAMING', 'RCS06_STREAMING',...
        ...
        'RCS_Weekly', 'RCS_Monthly',...
        'FlUCT'};

elseif nargin == 3

    rcs_TOKEN      = varargin{1};
    pcs_TOKEN      = varargin{2};
    pt_id_list     = varargin{3};


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
            PATIENT_ARM     = 'rcs01_daily_arm_1';
            reportid        = '73191';
        case 'RCS02'
            PATIENT_ARM     = 'rcs02_daily_arm_7';
            reportid        = '87050';
        case 'RCS04'
            PATIENT_ARM     = 'rcs04_daily_arm_10';
            reportid        = '104806';
        case 'RCS05'
            PATIENT_ARM      = 'rcs05_daily_arm_13';
            reportid         = '109667';

        case 'RCS06'
            PATIENT_ARM      = 'rcs06_pain_report_arm_27';
            reportid         = '135787';

        case 'RCS07'
            PATIENT_ARM      = 'rcs07_pain_report_arm_33';
            reportid         = '144629';
            
            % NEW arms
        case  'RCS02new'
            PATIENT_ARM = 'rcs02_new_pain_rep_arm_17';
            reportid = '112131';
        case  'RCS04new'
            PATIENT_ARM = 'rcs04_new_pain_rep_arm_18';
            reportid = '112132';

        case  'RCS05new'
            PATIENT_ARM      = 'rcs05_new_pain_rep_arm_19';
            reportid         = '112133';
            
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

        case 'RCS06_STREAMING'

            PATIENT_ARM      = 'rcs06_streaming_arm_28';
            reportid         = '135788';
            

        case 'RCS07_STREAMING'
      
            PATIENT_ARM      = 'rcs07_streaming_no_arm_36';
            reportid         = '146251';


        case 'FLUCT'
            PATIENT_ARM = {'dbs_and_nondbs_pat_arm_23', 'pain_fluctuation_s_arm_11'};
            reportid    = '128478';

        case 'RCS_Weekly'
            
            reportid = '135968';

        case 'RCS_Monthly'

            reportid = '135969';
            
            
        otherwise
            
            fprintf('\n Data not found for %s !     ...      Continuing ...  \n\n',pt_id)
            continue
            %
    end
    
    options = weboptions('TimeOut', 60);
    disp('************************');
    if strcmp(reportid , '128478') % FLUCT study requires different API token
        
        data = webwrite(...
            SERVICE,...
            'token', pcs_TOKEN, ...
            'content', 'report',...
            'report_id',reportid, ...
            'format', 'csv',...
            'type','flat',...
            'rawOrLabelHeaders','label',...
            'exportCheckboxLabel','false',...
            'exportSurveyFields','true',...
            'returnformat','csv', ...
            options);
        
        
        alltable = data;

    else
        
        data = webwrite(...
        SERVICE,...
        'token', rcs_TOKEN, ...
        'content', 'report',...
        'report_id',reportid, ...
        'format', 'csv',...
        'type','flat',...
        'rawOrLabelHeaders','raw',...
        'exportCheckboxLabel','false',...
        'exportSurveyFields','true',...
        'returnformat','csv',...
        options);
    
    
        alltable = data;
    end
   

    varnames = alltable.Properties.VariableNames;
    
    for i = 1 : length(varnames)
        col = varnames{i};
        
        % searchers for columns where the rows are cells, but w/n those
        % cells there is NOT character arrays

        % handles cases where datetimes, and doubles have frivolous braces
        if iscell(alltable.(col))
            if ~ischar(alltable.(col){1})
        
            alltable.(col) = cellfun(@(x) x(1:end), alltable.(col), 'UniformOutput', false);
            end
        end
    end


    if strcmp(pt_id, 'FLUCT')

       varnames = alltable.Properties.VariableNames;

        for i = 1 : length(varnames)
            if isnumeric(alltable.(varnames{i}))

                if all(isnan(alltable.(varnames{i})))
                    alltable = removevars(alltable, varnames{i});
                end
  
            elseif isdatetime(alltable.(varnames{i}))

                if all(isnat(alltable.(varnames{i})))
                     alltable = removevars(alltable, varnames{i});
                end
            end
        end

       alltable = removevars(...
                        alltable(contains(alltable.EventName, PATIENT_ARM),:),...
                        'EventName');
       %%
       raw_initals = alltable.PleaseEnterYourInitials_;

       tmp_initals = upper(replace(raw_initals, {' ','.'},{'',''}));

       rcs_initals.MB  = {'MB'};
       rcs_initals.AMW = {'AMW','AMW11PM 09-20-20', 'AME', 'ANW', 'AMP'};
       rcs_initals.EM  = {'EM'};    
       rcs_initals.SDZ = {'SDZ'};
       rcs_initals.SW  = {'SW'};

       initals = fieldnames(rcs_initals);
       pts                     = {'RCS02','RCS04', 'RCS05', 'RCS06', 'RCS07'};


       mpq_questions = {'Throbbing', 'Shooting', 'Stabbing',...
                            'Sharp', 'Cramping', 'Gnawing','Hot_burning',...
                            'Aching', 'Heavy', 'Tender', 'Splitting',...
                            'Tiring_Exhausting', 'Sickening', 'Fearful',...
                            'Cruel'};

       for i=1: length(initals)

            clntable = alltable(...
                       contains(tmp_initals, rcs_initals.(initals{i})), :);


            varnames = clntable.Properties.VariableNames;

            for j = 1 : length(varnames)
                if isnumeric(clntable.(varnames{j})) && ~contains(varnames{j}, mpq_questions)
    
                    if all(isnan(clntable.(varnames{j})))
                        clntable = removevars(clntable, varnames{j});
                    end
      
                elseif isdatetime(clntable.(varnames{j}))
                    
                    if all(isnat(clntable.(varnames{j})))
                         clntable = removevars(clntable, varnames{j});
                    end
                end
            end

            varnames     = clntable.Properties.VariableNames;
            complete_var = varnames(contains(varnames, 'Complete'));
            clntable     = clntable(clntable{:, complete_var} == 2, :);


            redcap = table;
            
            redcap.time = clntable.SurveyTimestamp;
            redcap.time.TimeZone = 'America/Los_Angeles';

            raw_initals    = clntable.PleaseEnterYourInitials_;
            redcap.initals = upper(replace(raw_initals, {' ','.'},{'',''}));

           switch initals{i}
               
               case {'MB', 'EM'}
                    %field names per subject
                    nrs_field         = 'PainIntensityRating_1';
                    nrs_unp_field     = 'PainUnpleasantness';
                    nrs_relief_field  = 'PainReliefRating_1';

                    nrs_relat_field   = 'RelativePainIntensityRating_pleaseRateYourPainRelativeToLastS_1';
                    
                    vas_field         = 'PainIntensityRating';
                    vas_relief_field  = 'PainReliefRating';
                    vas_relat_field   = 'RelativePainIntensityRating_pleaseRateYourPainRelativeToLastSun';

                    vas_unp_field     = 'PainUnpleasantnessRating';

                    redcap.reliefNRS     = (clntable.(nrs_relief_field));
                    redcap.relatNRS      = (clntable.(nrs_relat_field));
                    redcap.reliefVAS     = (clntable.(vas_relief_field));
                    redcap.relatVAS      = (clntable.(vas_relat_field));

               case {'AMW','SDZ', 'SW'}
                    %field names per subject
                    nrs_field         = 'PleaseRateYourCurrentPainIntensity_';
                    nrs_unp_field     = 'PleaseRateYourCurrentPainUnpleasantness__1';

                    vas_field         = 'PleaseRateYourCurrentPainIntensityBySlidingTheBoxAlongTheLine_';

                    vas_unp_field     = 'PleaseRateYourCurrentPainUnpleasantness_';

                    vas_mood_field    = 'PleaseRateYourCurrentMoodBySlidingTheBoxAlongTheLine_';
                
                    % No mood VAS was collected for some pts
                    redcap.moodVAS       = (clntable.(vas_mood_field));

                    redcap.reliefVAS     = nan(height(redcap),1);
               otherwise


           end

            %%%
            redcap.mayoNRS       = (clntable.(nrs_field));
            redcap.unpleasantNRS = (clntable.(nrs_unp_field));
            
            redcap.painVAS       = (clntable.(vas_field));
            redcap.unpleasantVAS = (clntable.(vas_unp_field));


            redcap.MPQthrobbing    = (clntable.Throbbing);
            redcap.MPQshooting     = (clntable.Shooting);
            redcap.MPQstabbing     = (clntable.Stabbing);
            redcap.MPQsharp        = (clntable.Sharp);
            redcap.MPQcramping     = (clntable.Cramping);
            redcap.MPQgnawing      = (clntable.Gnawing);
            redcap.MPQhot_burning  = (clntable.Hot_burning);
            redcap.MPQaching       = (clntable.Aching);
            redcap.MPQheavy        = (clntable.Heavy);
            redcap.MPQtender       = (clntable.Tender);
            redcap.MPQsplitting    = (clntable.Splitting);
            redcap.MPQtiring       = (clntable.Tiring_Exhausting);
            redcap.MPQsickening    = (clntable.Sickening);
            redcap.MPQfearful      = (clntable.Fearful);
            redcap.MPQcruel        = (clntable.Cruel_Punishing);

            i_mpq_fields           = contains(redcap.Properties.VariableNames, 'MPQ');
            MPQtotal               = sum(redcap{:, i_mpq_fields}, 2, 'omitnan');

            redcap = addvars(redcap, MPQtotal, 'Before', 'MPQthrobbing');

            % addresses cases when MPQtotal was likely unanswered
            redcap = redcap(~(redcap.MPQtotal == 0 & redcap.painVAS > 30), :);

%             save_dir = ...
%            ['/Users/Leriche/Dropbox (UCSF Department of Neurological Surgery)/',...
%                 'SUBNETS Dropbox/Chronic Pain - Activa and Summit 2.0/DATA ANALYSIS/',...
%                 'Ryan_Stage0_group_level_behavioral_analysis/processed_data/pain_fluctuation_study_RCS_pts/'];
% 
%             writetable(redcap, [save_dir,pts{i},'_pre_trial.xls'])

            painscores_out.(initals{i}) = redcap;
       end
        % rename as RCS pt code for clarity
       
        painscores_out          = cell2struct(struct2cell(painscores_out), pts);

       toc
       return
    end

    if ~contains(pt_id, {'STREAMING', 'Weekly', 'Monthly'}) 
        
        vars = alltable.Properties.VariableNames;        

        % if column was imported as cell array--> parse doubles and
        % datetimes as their own arrays


        for i=1:length(vars)

            if iscell(alltable.(vars{i}))
                if contains(vars{i}, 'timestamp')
                    
                    alltable.(vars{i}) = datetime(alltable.(vars{i}));    
                
                elseif ~ischar(alltable.(vars{i}){1}) % character arrays are fine in cells

                    alltable.(vars{i}) = [alltable.(vars{i}){:}]';

                end
            end
        end

        timestampvars = ...
            vars(find(contains(alltable.Properties.VariableNames,'timestamp')));
        

        %FOR PAIN SCORE ARMS
         % remove all the extraneous rows ('events' that are different)
        keeprows = strcmp(alltable.redcap_event_name, PATIENT_ARM) & ...
            (arrayfun(@(x) ~isnat(x),alltable.(timestampvars{1})) | ...
            arrayfun(@(x) ~isnat(x),alltable.(timestampvars{2})));


        %for some reason, all patients are using the field rcs01_mpq_...
        
        clntable = alltable(keeprows,:);
        
        % Moves redcap timestamps into designated VAS or MPQ segments. Then combines
        % into consolidated timestamp structure.
        
        redcap_timestamp.vasnrs          =  datetime(clntable.(timestampvars{1}));
        redcap_timestamp.mpq             = datetime(clntable.(timestampvars{2}));
        time_transfer                    = find(isnat(redcap_timestamp.vasnrs) & ~isnat(redcap_timestamp.mpq));
        redcap_timestamp.alltimes        = redcap_timestamp.vasnrs;
        redcap_timestamp.alltimes(time_transfer) = redcap_timestamp.mpq(time_transfer);
        
        %%%%%%%%% Define Time and Pain Scores %%%%%%%%%%%
        redcap.time = redcap_timestamp.alltimes;
        % specify timezone for explicit datetime variable
        redcap.time.TimeZone = 'America/Los_Angeles'; % corresponds to 'America/Los_Angeles' timezone;
   
              
        varnames = clntable.Properties.VariableNames;
    
         switch pt_id
            case 'RCS06' 
            
                redcap.mayoNRS = (clntable.('intensity_nrs_v2_v2_51a1a3'));
                
                redcap.nocNRS = (clntable.('worst_please_rate_your_pain_inte_v2_37cdf10'));
                
                redcap.npNRS = (clntable.('worst_please_rate_your_pain_inte_v2_37cdf9'));
                
                redcap.painVAS = (clntable.('intensity_vas_v2_v2_8a3273'));
                
                redcap.nocVAS = (clntable.('please_rate_your_pain_inte_v2_9d8b34'));
                redcap.npVAS = (clntable.('please_rate_your_pain_inte_v2_9d8b33'));
                
                
                unp_field = varnames{contains(varnames,'unpleasantness')};
                
                redcap.unpleasantVAS = (clntable.(unp_field));
            
            case 'RCS07'
            
                % Populate the new flavors of painscores downloaded from redcap
                redcap.mayoNRS          = clntable.('overall_intensity_nrs');
                redcap.unpNRS           = clntable.('unp_intensity_nrs');
                redcap.unp_abbrev_NRS   = clntable.('unpleasantness_abbrev');
                
                
                redcap.leftarmNRS       = clntable.('left_arm_intensity_nrs');
                redcap.leftlegNRS       = clntable.('left_leg_intensity_nrs');
                redcap.leftfaceNRS      = clntable.('left_face_intensity_nrs');

 
                redcap.painVAS          = clntable.('pain_intensity_vas');
                redcap.unpleasantVAS    = clntable.('unp_intensity_vas');
                
                redcap.moodVAS          = clntable.('mood_intensity_vas');
                

             case 'RCS05new'
            
                nrs_field = varnames{contains(varnames,'intensity_nrs')};
                vas_field = varnames{contains(varnames,'intensity_vas')};
                
                unp_field = varnames{contains(varnames,'unpleasantness')};
                
                % Populate the new flavors of painscores downloaded from redcap
                redcap.mayoNRS       = clntable.(nrs_field);
                redcap.painVAS       = clntable.(vas_field);
                redcap.unpleasantVAS = clntable.(unp_field);


                redcap.painSpikeBoolean = clntable.pain_spike_y_n;
                redcap.painSpikeNRS     = clntable.pain_spike_nrs;
 
            otherwise
                % field names per subject
                nrs_field = varnames{contains(varnames,'intensity_nrs')};
                vas_field = varnames{contains(varnames,'intensity_vas')};
                
                unp_field = varnames{contains(varnames,'unpleasantness')};
                
                % Populate the new flavors of painscores downloaded from redcap
                redcap.mayoNRS = (clntable.(nrs_field));
                redcap.painVAS = (clntable.(vas_field));
                
                redcap.unpleasantVAS = (clntable.(unp_field));

                redcap.painSpikeBoolean  = nan(numel(redcap.mayoNRS),1);
                redcap.painSpikeNRS     = nan(numel(redcap.mayoNRS),1);
      
         end

        
        if contains(pt_id,'new')
            worstnrs_field = varnames{contains(varnames,'worst_please_rate')};
            worstvas_field = varnames{contains(varnames,'please_rate_your_pain_inte') & ~contains(varnames,'worst')};
                        
            redcap.worstNRS = (clntable.(worstnrs_field));
            redcap.worstVAS = (clntable.(worstvas_field));
            
        else
            redcap.worstNRS = nan(numel(redcap.mayoNRS),1);
            redcap.worstVAS = nan(numel(redcap.mayoNRS),1);
        end
        
        
        % take the mcgill pain questionnaires
        table_header = clntable.Properties.VariableNames;
        mpq_start = find(contains(table_header,'rcs01_mpq_timestamp'))+1;
        mpq_end = find(contains(table_header,'cruel_punishing'));
        
        mpqhold = table2array(clntable(:,mpq_start:mpq_end));
        
        redcap.MPQtotal       = tsnansum((mpqhold),2);
        redcap.MPQthrobbing   = (clntable.throbbing);
        redcap.MPQshooting    = (clntable.shooting);
        redcap.MPQstabbing    = (clntable.stabbing);
        redcap.MPQsharp       = (clntable.sharp);
        redcap.MPQcramping    = (clntable.cramping);
        redcap.MPQgnawing     = (clntable.gnawing);
        redcap.MPQhot_burning = (clntable.hot_burning);
        redcap.MPQaching      = (clntable.aching);
        redcap.MPQheavy       = (clntable.heavy);
        redcap.MPQtender      = (clntable.tender);
        redcap.MPQsplitting   = (clntable.splitting);
        redcap.MPQtiring      = (clntable.tiring_exhausting);
        redcap.MPQsickening   = (clntable.sickening);
        redcap.MPQfearful     = (clntable.fearful);
        redcap.MPQcruel       = (clntable.cruel_punishing);

        
        % make 0 values NaN for RCS01 because of unreliable reporting
        if strcmp(pt_id,'RCS01')
            redcap.MPQtotal(redcap.MPQtotal == 0) = NaN;
        else
            redcap.MPQtotal = redcap.MPQtotal;
        end
        
        % if NRS is Nan--rather than zero--likely means survey was quickly opened/closed
        % therefore changed MPQtotal to NaN rather than 0
        redcap.MPQtotal(isnan(redcap.mayoNRS)) = NaN;

        painscores_out.(pt_id) = redcap;
        

        
    elseif ~contains(pt_id, {'Weekly', 'Monthly'}) 
        %FOR STREAMING NOTES ARMS
        vars = alltable.Properties.VariableNames;

         for i=1:length(vars)

            if iscell(alltable.(vars{i}))
                if contains(vars{i}, 'timestamp')
                    
                    alltable.(vars{i}) = datetime(alltable.(vars{i}));    
                
                elseif ~ischar(alltable.(vars{i}){1}) % character arrays are fine in cells

                    alltable.(vars{i}) = [alltable.(vars{i}){:}]';

                end
            end
         end
        
         timevarNAME = alltable.Properties.VariableNames( ...
            contains(alltable.Properties.VariableNames,'timestamp'));
        
        
        
        
        keeprows = strcmp(alltable.redcap_event_name, PATIENT_ARM) & ...
            (arrayfun(@(x) ~isnat(x),alltable.(timevarNAME{1})));
              
        clntable = alltable(keeprows,:);
        
        varnames = clntable.Properties.VariableNames;

        if contains(pt_id, 'RCS07')

            medchange_field = varnames{contains(varnames,'changes_to')};
            activity_field = varnames{contains(varnames,'activit')};

            explain_field = varnames{contains(varnames,'please_note_changes')};
            stimprog_field = varnames{contains(varnames,'stimulation_pro')};


        else

            medchange_field = varnames{contains(varnames,'medication_changes')};
            activity_field = varnames{contains(varnames,'activity')};

            explain_field = varnames{contains(varnames,'explain')};
            stimprog_field = varnames{contains(varnames,'which_stimulation_program')};

        end

            stimon_field = varnames{contains(varnames,'is_stimulation_on')};


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
 
%newscores.RCS01 =  oldscores.RCS01;

newscores.RCS02 = [oldscores.RCS02; oldscores.RCS02new];
newscores.RCS04 = [oldscores.RCS04; oldscores.RCS04new];
newscores.RCS05 = [oldscores.RCS05; oldscores.RCS05new];

newscores.RCS06 = oldscores.RCS06;
newscores.RCS07 = oldscores.RCS07;
    
    
%newscores.RCS01_STREAMING = oldscores.RCS01_STREAMING;
newscores.RCS02_STREAMING = oldscores.RCS02_STREAMING;
newscores.RCS04_STREAMING = [oldscores.RCS04_STREAMING; oldscores.RCS04_STREAMING_v2];
newscores.RCS05_STREAMING = oldscores.RCS05_STREAMING;

newscores.RCS06_STREAMING = oldscores.RCS06_STREAMING;
newscores.RCS07_STREAMING = oldscores.RCS07_STREAMING;

% OUTPUT
painscores_out = newscores;

toc