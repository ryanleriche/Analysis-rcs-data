function plot_timeline_no_PFS(cfg, REDcap, varargin)
%% Basic Information
% Index behavioral data of patient of interest 
redcap                  = REDcap.(cfg.pt_id);

% List out stage dates
cfg.stage_dates         = cfg.stage_dates{str2double(cfg.pt_id(end))}; % starts at Stage 1


% Define figure dimensions and dates
figure('Units', 'Inches', 'Position', [0, 0, 13, 9])

[~, date_range] = date_parser(cfg, redcap);

ds = datestr(date_range,'dd-mmm-yyyy');


%% NRS

    % Define plot dimensions     
    if cfg.subplot == true
        subplot(311);
        sgtitle([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);

    else
        title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);
        hold on
    end

hold on
    
    % Plotting behavioral data 
    switch cfg.dates
        case 'AllTime'
            
                % Plot moving mean based on 5 surveys 
                plot(redcap.time, movmean([redcap.mayoNRS,redcap.worstNRS], 5),...
                    'LineWidth', 2, 'HandleVisibility','off');
                
                % Set color of moving mean to match behavioral plot
                hold on; set(gca,'ColorOrderIndex',1);
                
                % Plot behavioral plot
                plot(redcap.time, [redcap.mayoNRS,redcap.worstNRS], '.', 'MarkerSize',5);

        case {'DateRange', 'PreviousDays'}
            
            % Number of surveys used in the moving mean
            n_back  = 21;
    
            % Plot the moving mean
            plot(redcap.time, movmean(redcap.mayoNRS, [n_back 0], 'omitnan'),...
                'LineWidth', 2,'HandleVisibility','off'); hold on;
            
            % Set moving mean and RECAP NRS values to the same color
            set(gca,'ColorOrderIndex',1)
    
            % Plot the RECAP NRS values
            scatter(redcap.time, redcap.mayoNRS, 125, 'filled','MarkerFaceAlpha', 0.6);    
        
        % Specific cases for RCS06 and RCS07 for their unique pain
        % survey questions 
        switch cfg.pt_id

            case 'RCS06'
                set(gca,'ColorOrderIndex',3)
                
                % Nociceptive NRS (nocNRS)
                plot(redcap.time, movmean(redcap.nocNRS, [n_back 0], 'omitnan'),...
                    'LineWidth', 2, 'HandleVisibility','off'); hold on;
    
                    set(gca,'ColorOrderIndex',3)
                scatter(redcap.time, redcap.nocNRS, 75, 'filled','MarkerFaceAlpha', 0.6); 
    
                % Neuropathic NRS (npNRS)
                plot(redcap.time, movmean(redcap.npNRS, [n_back 0], 'omitnan'),...
                    'LineWidth', 2, 'HandleVisibility','off'); hold on;
    
                    set(gca,'ColorOrderIndex',4)
                scatter(redcap.time, redcap.npNRS, 50, 'filled', 'MarkerFaceAlpha', 0.6);

           
            case 'RCS07'
                
                % Plot scatter plots for unpleasantness, left arm, face, and
                % leg NRS
                scatter(redcap.time, redcap.unpNRS, 50, 'filled','MarkerFaceAlpha', 0.6); 
                scatter(redcap.time, redcap.leftarmNRS, 65, 'filled', 'MarkerFaceAlpha', 0.6);
    
                scatter(redcap.time, redcap.leftlegNRS, 50 , 'filled', 'MarkerFaceAlpha', 0.6);  
                scatter(redcap.time, redcap.leftfaceNRS, 35, 'filled', 'MarkerFaceAlpha', 0.6);


            case {'RCS02', 'RCS04', 'RCS05'}
                % Plot scatter plots for worst NRS 
                scatter(redcap.time, redcap.worstNRS, 75, 'filled', 'MarkerFaceAlpha', 0.6);

        end
    end
  
     % Define y-axis
     ylabel('Numeric Rating Scale');     ylim([0,10]); yticks(1:2:10);

     % Define legends for specific patients
     switch cfg.pt_id
         case 'RCS06'

            legend({'NRS Intensity','NRS Nociceptive','NRS Neuropathic',}, ...
                    'Location','northoutside', 'NumColumns',3);

         case 'RCS07'

            legend({'NRS Intensity','NRS Unpleasantness',...
                    'NRS Left Arm', 'NRS Left Leg', 'NRS LeftFace'},...
                    'Location','northoutside', 'NumColumns',3);
         
         case {'RCS02', 'RCS04', 'RCS05'}

            legend({'NRS Intensity' , 'NRS Worst Intensity'},...
                    'Location','northoutside', 'NumColumns',3);
     end
  
     format_plot();
  
     % Define stimulation parameter (previously embedded stim on graph)
     if ~strcmp(cfg.stim_parameter, '')
         overlay_stim(gca, cfg.stim_parameter, redcap);   
     end
     

%% VAS

    % Define plot dimensions
    if cfg.subplot == true
        subplot(312);
    
    else
        figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
        title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16); 
        hold on

    end

hold on

    % Plot behavioral data
    switch cfg.dates
        case 'AllTime'
            % Plot moving mean of 5 surveys
            plot(redcap.time, movmean([redcap.painVAS,redcap.unpleasantVAS,redcap.worstVAS], 5),...
                'LineWidth', 2, 'HandleVisibility','off');
    
            % Make moving mean same color as behavioral data points
            set(gca,'ColorOrderIndex',1);
            
            % Plot behavioral data
            plot(redcap.time, [redcap.painVAS,redcap.unpleasantVAS,redcap.worstVAS], '.', 'MarkerSize',5);
        
        case {'DateRange', 'PreviousDays'}

            scatter(redcap.time, redcap.painVAS, 100, 'filled', 'MarkerFaceAlpha', 0.6);

        % Specific cases for RCS06 and RCS07 for their unique pain
        % survey questions 
        switch cfg.pt_id

            case 'RCS06'
                
                % Plot moving mean for nociceptive pain VAS
                set(gca,'ColorOrderIndex',3);
                    plot(redcap.time, movmean(redcap.nocVAS, [n_back 0], 'omitnan'),...
                'LineWidth', 2.5,'HandleVisibility','off');
    
                % Plot behavioral data for nociceptive pain VAS
                set(gca,'ColorOrderIndex',3);
                    scatter(redcap.time, redcap.nocVAS, 75, 'filled', 'MarkerFaceAlpha', 0.6);
                
                % Plot moving mean for neuropathic pain VAS
                set(gca,'ColorOrderIndex',4);
                    plot(redcap.time, movmean(redcap.npVAS, [n_back 0], 'omitnan'),...
                    'LineWidth', 1.5, 'HandleVisibility','off');
    
                % Plot behavioral data for neuropathic pain VAS
                set(gca,'ColorOrderIndex',4);
                    scatter(redcap.time, redcap.npVAS, 50, 'filled', 'MarkerFaceAlpha', 0.6);


            case {'RCS02', 'RCS04', 'RCS05'}
                set(gca,'ColorOrderIndex',1);
            
                % Plot moving mean for pain VAS
                plot(redcap.time, movmean(redcap.painVAS, [n_back 0], 'omitnan'),...
                'LineWidth', 2.5, 'HandleVisibility','off');
                
                % Plot behavioral data for pain VAS
                scatter(redcap.time, redcap.unpleasantVAS, 75, 'filled','MarkerFaceAlpha', 0.6);
                scatter(redcap.time, redcap.worstVAS, 50, 'filled','MarkerFaceAlpha', 0.6);

            case  'RCS07'
                    
                set(gca,'ColorOrderIndex',1);
                
                % Plot moving mean for pain VAS
                plot(redcap.time, movmean(redcap.painVAS, [n_back 0], 'omitnan'),...
                    'LineWidth', 2.5, 'HandleVisibility','off');
                
                % Plot behavioral data for unpleasantness pain VAS
                scatter(redcap.time, redcap.unpleasantVAS, 75, 'filled','MarkerFaceAlpha', 0.6);
            
                set(gca,'ColorOrderIndex',4);

                % Plot moving mean for mood VAS
                plot(redcap.time, movmean(redcap.moodVAS, [n_back 0], 'omitnan'),...
                'LineWidth', 2,'HandleVisibility','off'); hold on;
                
                % Plot behavioral data for mood VAS 
                set(gca,'ColorOrderIndex',4);
                scatter(redcap.time, redcap.moodVAS, 75, 'filled','MarkerFaceAlpha', 0.6);

        end
    end

    % Define y-axis 
     ylabel('Visual Analog Scale');         ylim([0,100]);  yticks(0:20:100);

    % Define legends for specific patients
    switch cfg.pt_id
         case 'RCS06'
            
            legend({'VAS Intensity', 'VAS Nociceptive','VAS Neuropathic'},...
                    'Location','northoutside', 'NumColumns',3);

         case 'RCS07'

            legend({'VAS Intensity', 'VAS Unpleasantness','VAS Mood'},...
                    'Location','northoutside', 'NumColumns',3);
         
         case {'RCS02', 'RCS04', 'RCS05'}

            legend({'VAS Intensity', 'VAS Unpleasantness','VAS Worst'},...
                    'Location','northoutside', 'NumColumns',3);
     end

     % Define stimulation parameter (previously embedded stim on graph)
     if ~strcmp(cfg.stim_parameter, '')
         overlay_stim(gca, cfg.stim_parameter, redcap);   
     end

     format_plot();

            
%% MPQ

    % Define plot dimensions
    if cfg.subplot == true

        subplot(313);
    else
        figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
        title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16); 
        hold on

    end
    
    % Calculate affective and somatic MPQ scores
    MPQaff       = sum([redcap.MPQsickening, redcap.MPQfearful, redcap.MPQcruel, redcap.MPQtiring],2,'omitnan');
    MPQsom       = redcap.MPQtotal - MPQaff;

    % Plotting behavioral data
    switch cfg.dates
        case 'AllTime'
           
            % Plot moving mean for MPQ total, somatic, and affective
            plot(redcap.time, movmean(redcap.MPQtotal , 5),'LineWidth', 4.5,'HandleVisibility','off');
            hold on;
            plot(redcap.time, movmean([MPQsom, MPQaff], 5),'LineWidth', 2,'HandleVisibility','off');
    
             set(gca,'ColorOrderIndex',1);
            
            % Plot behavioral data for MPQ total, somatic, and affective
            plot(redcap.time, redcap.MPQtotal , '.', 'MarkerSize',9);
            plot(redcap.time, [MPQsom, MPQaff], '.', 'MarkerSize',5);

        case {'DateRange', 'PreviousDays'}

            % Plot moving mean for MPQ total
            plot(redcap.time, movmean(redcap.MPQtotal, [n_back 0], 'omitnan'),...
                'LineWidth', 2,'HandleVisibility','off');
    
            % Set moving mean color same to MPQ total
            set(gca,'ColorOrderIndex',1);
            
            % Plot behavioral data for MPQ scores
            scatter(redcap.time, redcap.MPQtotal , 75, 'filled','MarkerFaceAlpha', 0.6);   
            hold on;
            scatter(redcap.time, MPQsom, 50, 'filled', 'MarkerFaceAlpha', 0.6);
            scatter(redcap.time, MPQaff, 50, 'filled', 'MarkerFaceAlpha', 0.6);

    end
    
    % Define y-axis
    ylabel('McGill Pain Questionaire');     ylim([0,45]); %yticks(0:5:45)
    
    legend({'MPQ Total',  'MPQ Somatic','MPQ Affective',}, 'Location','northoutside', 'NumColumns',3);
    
    format_plot();

    if ~strcmp(cfg.stim_parameter, '')
        overlay_stim(gca, cfg.stim_parameter, redcap);   
    end

%% Indicates stage date lines 

function format_plot()  

    % Define grids 
    grid on;    grid MINOR;   
    
    legend boxoff;    box off;

    set(gca,'FontSize',12, 'xlim', date_range, 'TickLength', [0 0],...
        'GridAlpha',0.4,'MinorGridAlpha',0.7, 'GridColor', 'k', 'MinorGridColor', 'k'); 

    t = datetime(cfg.stage_dates, 'TimeZone', '-07:00');

    if length(cfg.stage_dates) == 1
         xline(t, '-', {'Stage 1'},...
            'HandleVisibility','off');

    elseif length(cfg.stage_dates) == 2
        xline(t(1), '-', {'Stage 1'},...
            'HandleVisibility','off');
        
        xline(t(2), '-', {'Stage 2'},...
            'HandleVisibility','off');

    else
         xline(t(1), '-', {'Stage 1'},...
            'HandleVisibility','off');
        
        xline(t(2), '-', {'Stage 2'},...
            'HandleVisibility','off');
       
         xline(t(3), '-', {'Stage 3'},...
            'HandleVisibility','off');
    
    end
end
end
