function plt_percent_ON_overtime(t_vec, on_off_vec,...
                                step_dur, wash_out_tbl, redcap, colors,...
                                t_plt_start, t_plt_end)


    plot(t_vec,   movmean(on_off_vec, [duration('01:00:00')/step_dur,0]), 'k', 'LineWidth',2); hold on
            
    avg_per_on = mean(on_off_vec, 'omitnan');

    plot(t_vec, repmat(avg_per_on, length(t_vec), 1),'--', 'Color', 'k', 'LineWidth', 1.5);


    if le(t_plt_end-t_plt_start, duration('60:00:00'))

        leg_txt = {'Lagging Mean of 1 hour', ...
                    sprintf('Mean Percent On (%.1f%%)', avg_per_on),...
                   'Lagging Mean of 1 min', 'aDBS Off', ...
                   'NRS intensity', 'VAS intensity', 'MPQ total'};

        stairs(t_vec, movmean(on_off_vec, [duration('00:01:00')/step_dur,0]),'Color', colors(3,:));

    else


        leg_txt = {'Lagging Mean of 1 hour',...
                    sprintf('Mean Percent On (%.1f%%)', avg_per_on),...
                    'aDBS Off', ...
                    'NRS intensity', 'VAS intensity', 'MPQ total'};

        
    end

    
    
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
    
    pain_metrics =  {'mayoNRS', 'painVAS', 'MPQtotal'};
    
    for i_pain =1:size(pain_metrics,2)
        plot(plt_rcap.time, movmean(plt_rcap.(pain_metrics{:,i_pain}), 3),...
            'LineWidth', 3, 'Color', colors(i_pain+3,:),...
            ...
            'Marker', 'none','LineStyle', '-',...
            'HandleVisibility','off'); hold on
    
        scatter(plt_rcap.time, plt_rcap.(pain_metrics{:,i_pain}), 130, ...
            'filled','MarkerFaceAlpha', 0.6, ...
            'MarkerEdgeColor', colors(i_pain+3,:), 'MarkerFaceColor', colors(i_pain+3,:)...
            );   
    end
    
    
    y_lbls_txt = cellfun(@(x) sprintf('%g   %g   %g',x), ...
    num2cell([0:10; 0:10:100; 0:4.5:45]', 2), 'UniformOutput', false);
    
    yticks(0:10); ylim([0,10]);  yticklabels(y_lbls_txt);
    
    fig_h = gcf;

    fig_h.CurrentAxes.YAxis(1).FontSize = 14;     
    fig_h.CurrentAxes.YAxis(2).FontSize = 14; 
    
    fig_h.CurrentAxes.XAxis.FontSize    = 14;
    
    legend(leg_txt, 'FontSize',14, 'Location', 'northoutside', 'NumColumns', 6);
    
    grid on;

end