function plt_state_mA_overtime(fig_h, plt_app_oi, t_plt_start, t_plt_end, x_tick_dur, h)

% ignore state "15" as this is therapyStatus Off
plt_app_oi.oldstate(plt_app_oi.oldstate == 15) = NaN;


% specific sizing
fig_h.CurrentAxes.Position([2, 4]) = [0.1, 0.22];

% create matrix w/ any state change time as columns and state as rows
state_chan     = plt_app_oi.oldstate(2:end);
state_chan_mat = nan(9, length(state_chan));

% keeping state itself to visualize w/ different colors per state
for j = 0:8
    state_chan_mat(j+1, state_chan == j) = state_chan(state_chan == j);
end

if ~isempty(state_chan) && height(plt_app_oi.time_INS) >2

   % state_col = brewermap(max(state_chan)+1, 'Spectral');
    
    t_chan        = plt_app_oi.time_INS(1:end-1);
    
     % now put current on second y-axis to intuitively see mA:state definition

     
    if h ~= 0

        yyaxis right
        
        set(gca, 'YColor', 'k');

        amp_chan = plt_app_oi.prog0mA(2:end);
        
        h_stair =  stairs(t_chan, amp_chan ,'-','LineWidth',.1, 'Color',  'k');    hold on
        
        alpha(0.1); 

        yticks(0:0.5:max(amp_chan)+.5);            ylim([0, max(amp_chan)+.5]); 
        
        ylabel('Current (mA)', 'Color', 'k');   

        yyaxis left
    end



    h_state = pcolor(t_chan,0:8, state_chan_mat);    colormap("lines")
    
    % remove pcolor gridlines and shift y-labels to be centered on "state
    % raster plot"
    set(h_state, 'EdgeColor', 'none');
    
    yticks(0.5:8.5); yticklabels(compose('%g', 0:8)); ylim([0,9]);
    
    ylabel('States (0–8)');
    
    set(gca, 'FontSize', 12);  grid on;
    

    set(gca,'TickLength', [0, 0],'FontSize', 12, 'YColor', 'k', ...
      'XLim', [t_plt_start, t_plt_end],...
      'GridAlpha',0.3);
    
% include xticks *exactly* like lagging average current plot
    xticks(t_plt_start:x_tick_dur:t_plt_end) 
end

end