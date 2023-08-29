function plt_percent_ON_overtime(t_vec, on_off_vec,...
                                step_dur, wash_out_tbl, redcap, colors,...
                                t_plt_start, t_plt_end, x_tick_dur, pt_side_id)


    plot(t_vec,   movmean(on_off_vec, [duration('01:00:00')/step_dur,0]), 'k', 'LineWidth',1.5); hold on
            
    i_toi      = isbetween(t_vec, t_plt_start, t_plt_end);
    avg_per_on = mean(on_off_vec(i_toi), 'omitnan');

    plot(t_vec, repmat(avg_per_on, length(t_vec), 1),'--', 'Color', 'k', 'LineWidth', 1.5);


    % attempt to show washout from consensus from 
    % EventLog.txt and DeviceSettings.json files
    
    for k = 1:height(wash_out_tbl)
        % only include handle of inital patch for ease of plotting
        % legend (legend accounts for EVERY object handle)
        if k == 1
            patch([wash_out_tbl.start(k), wash_out_tbl.start(k), wash_out_tbl.stop(k), wash_out_tbl.stop(k)],...
                   [0,100,100,0],[0.7, 0.7,0.7], ...
                    'FaceAlpha',0.5,'EdgeColor', 'none');  hold on
    
        else
            patch([wash_out_tbl.start(k), wash_out_tbl.start(k), wash_out_tbl.stop(k), wash_out_tbl.stop(k)],...
                   [0,100,100,0],[0.7, 0.7,0.7], ...
                'FaceAlpha',0.5,'EdgeColor', 'none',...
                'HandleVisibility','off'); 
        end
    end
    
    ylabel(['Percent time ON', newline,'(stim amplitude > 0 mA)'], 'FontSize',18)
    
    ylim([-5,105]);        
    
    yyaxis right; 
    
    plt_rcap = redcap;
    plt_rcap.painVAS  = redcap.painVAS/10;
    plt_rcap.MPQtotal = redcap.MPQtotal/4.5;
    
    switch pt_side_id(1:end-1)

        case 'RCS05'
            pain_metrics =  {'painSpikeNRS', 'mayoNRS', 'painVAS', 'MPQtotal'};
            symbols       = {'*','o', '+', '^'};

        case  'RCS06'
            pain_metrics =  {'npNRS','nocNRS', 'painVAS', 'MPQtotal'};
            symbols       = {'*','o', '+', '^'};
        otherwise
            pain_metrics =  {'mayoNRS', 'painVAS', 'MPQtotal'};
            symbols       = {'o', '+', '^'};


    end
    
    for i_pain =1:size(pain_metrics,2)
 

        switch pain_metrics{:, i_pain}

            case {'painSpikeNRS', 'npNRS'}
                 plt_color = 'r';
                 plt_size  = 200;

            otherwise
                 plt_color = 'b';
                 plt_size  = 75;
    
        end
        scatter(plt_rcap.time, plt_rcap.(pain_metrics{:,i_pain}), ...
                plt_size, symbols{i_pain}, ...
                'MarkerEdgeColor', plt_color,...
                'MarkerFaceColor', 'none',...
                'LineWidth',1.7...
            );
        hold on
    end
    
    
    y_lbls_txt = cellfun(@(x) sprintf('%g   %g   %g',x), ...
    num2cell([0:10; 0:10:100; 0:4.5:45]', 2), 'UniformOutput', false);
    
    yticks(0:10); ylim([0,10]);  yticklabels(y_lbls_txt);
    
    fig_h = gcf;

    fig_h.CurrentAxes.YAxis(1).FontSize = 14;     
    fig_h.CurrentAxes.YAxis(2).FontSize = 14; 
    
    fig_h.CurrentAxes.XAxis.FontSize    = 14;
    

    leg_txt = [{'Lagging Mean of 1 hour'},...
            {sprintf('Mean Percent On (%.1f%%)', avg_per_on)},...
            {'aDBS Off'}, ...
           pain_metrics(:)'];


    legend(leg_txt, 'FontSize',14, 'Location', 'northoutside', 'NumColumns', 8);
    
    grid on;

    xticks(t_plt_start:x_tick_dur:t_plt_end);


end