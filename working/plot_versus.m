function plot_versus(cfg, RCSXX)

 figure('Units', 'Inches', 'Position', [0, 0, 15, 10])

    [RCSXX, date_range] = date_parser(cfg, RCSXX);

    ds  =    datestr(date_range,'dd-mmm-yyyy');
    title([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',16);
    hold on

    RCSXX.MPQtotal    = RCSXX.MPQsum;
    RCSXX.MPQaff       = sum([RCSXX.MPQsickening, RCSXX.MPQfearful, RCSXX.MPQcruel],2,'omitnan');
    RCSXX.MPQsom       = RCSXX.MPQtotal - RCSXX.MPQaff;



    RCSXX_trim_VAS_all = RCSXX(RCSXX.unpleasantVAS ~= 50 &...
                            RCSXX.painVAS ~= 50 & RCSXX.worstVAS ~= 50,:);



    RCSXX_rmv_VAS_all = RCSXX(RCSXX.unpleasantVAS == 50 &...
                   RCSXX.painVAS == 50 & RCSXX.worstVAS == 50,:);

    RCSXX_rmv_VAS_any = RCSXX(RCSXX.unpleasantVAS == 50 |...
                   RCSXX.painVAS == 50 | RCSXX.worstVAS == 50,:);




    n_any_VAS_remain = sum(RCSXX.unpleasantVAS ~= 50 &...
                   RCSXX.painVAS ~= 50 & RCSXX.worstVAS ~= 50);

    prop_any_VAS_remain = n_any_VAS_remain /height(RCSXX);


    n_all_VAS_remain = sum(RCSXX.unpleasantVAS ~= 50 |...
                   RCSXX.painVAS ~= 50 | RCSXX.worstVAS ~= 50);

    prop_all_VAS_remain = n_all_VAS_remain /height(RCSXX);


    

%     plotmatrix([RCSXX.mayoNRS, RCSXX.worstNRS, RCSXX.painVAS,...
%         RCSXX.unpleasantVAS, RCSXX.worstVAS, MPQ_som, MPQ_aff]);
% 
%      plotmatrix([RCSXX.mayoNRS, RCSXX.painVAS, MPQ_som]);
% 
%      plotmatrix([MPQ_som,...
%         RCSXX.unpleasantVAS, RCSXX.worstVAS])

%%
    scatter3(RCSXX_trim_VAS_all,'unpleasantVAS','painVAS','MPQsom','filled', ...
        'ColorVariable', 'mayoNRS');


    xlim([1 100]); ylim([1 100]); zlim([0 45]);

    c = colorbar;
    c.Label.String = 'mayoNRS';                 c.Limits = [1,10];

    text(50, 15, 40, ...
        [...
            'Prop(remain w/o ANY VAS = 50s): ', num2str(prop_any_VAS_remain),...
            newline,...
            'Prop(remain w/o ALL VAS = 50): ', num2str(prop_all_VAS_remain)...
        ],...
        'FontSize',14);

    
    scatter3(RCSXX_rmv_VAS_any,'unpleasantVAS','painVAS','MPQsom','filled','MarkerFaceColor',[.7 .7 .7]);

    scatter3(RCSXX_rmv_VAS_all,'unpleasantVAS','painVAS','MPQsom','filled','MarkerFaceColor','k');

  


    format_plot()
  

    

function format_plot()  

    set(gca,'fontSize',14, 'TickLength', [0 0]); 
    grid on;    grid MINOR;      box off
    
end
end