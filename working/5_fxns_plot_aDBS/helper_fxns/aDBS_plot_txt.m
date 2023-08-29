function aDBS_plot_txt(fig_h, t_plt_start, t_plt_end, pt_side_id, meta, plt_app_oi, h)
        

set(gca,'xlim', [t_plt_start, t_plt_end],...
    'GridAlpha',0.3,...
    'YColor', 'k', 'TickLength', [0,0], 'FontSize', 12);


% reduce axis to 65% of its size so sense, and LD data can fit
% above
fig_h.CurrentAxes.Position(2)        = fig_h.CurrentAxes.Position(2) *.65;

% format xtick labels to fit on multiple lines
if h == 0 && ge(t_plt_end-t_plt_start, duration('72:00:00'))
    fig_h.CurrentAxes.XTickLabelRotation = 0;

    tmp_tick = cellfun(@(x) split(x,','), fig_h.CurrentAxes.XTickLabel, 'UniformOutput', false);

    if numel(tmp_tick(1)) >1
        new_tick = cellfun(@(x) sprintf('%s\\newline%s\n', x{1},x{2}), tmp_tick, 'UniformOutput', false);
        xticklabels([new_tick{:}]);
    end
end
%%%



dur_range =  t_plt_end -  t_plt_start;
    
text(t_plt_start - dur_range/15, fig_h.CurrentAxes.YLim(2)*1.9, pt_side_id , 'FontSize', 32);

text(t_plt_start - dur_range/15, fig_h.CurrentAxes.YLim(2)*1.4, meta.sense, 'FontSize', 10);

text(t_plt_start + dur_range/5.5, fig_h.CurrentAxes.YLim(2)*1.8, meta.by_ld0_pb, 'FontSize',10, 'Interpreter','none');
text(t_plt_start + dur_range/5.5, fig_h.CurrentAxes.YLim(2)*1.45, meta.by_ld1_pb, 'FontSize',10, 'Interpreter','none');



text(t_plt_end - dur_range/1.6, fig_h.CurrentAxes.YLim(2)*1.7, meta.ld0, 'FontSize',9, 'Interpreter','none');
text(t_plt_end - dur_range/2.5, fig_h.CurrentAxes.YLim(2)*1.7, meta.ld1, 'FontSize',9, 'Interpreter','none');

text(t_plt_end - dur_range/8,...
          fig_h.CurrentAxes.YLim(2) *1.6, ...
          ...
          [meta.state, sprintf('    GroupDrateInHz | %g', mode(plt_app_oi.rateHz(plt_app_oi.prog0mA > 0)))],...
          ...
          'FontSize',10, 'Interpreter','none');

end