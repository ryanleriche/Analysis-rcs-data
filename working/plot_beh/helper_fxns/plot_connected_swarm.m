function plot_connected_swarm(cfg, ind_by_group)


[xpos, ypos, ~, ~, ~]= UnivarScatter(ind_by_group,... 
    'Compression',50,... distance btwn boxes
    'Width',1,...
    'MarkerFaceColor',[0 0 0],...
    'PointSize', 12);


hold on
% boxplots up, now connect points across the boxplts
h_mean = plot(xpos.',ypos.','Linewidth',1,'Color','k');

% uncomment to make individual subject lines colored
%{
pompadour          = [104  2  63];
deep_sea           = [0  70 60];
deep_pink          = [192  11  111];
persian_rose       = [239   0  150];
persian_green      = [0  160  144];
aqua_marine        = [0  220  181];
carissma           = [255  149  186];
dark_olive         = [61  60   4];
aqua               = [95  255  222 ];
royal_blue         = [0  60  134];
indigo             = [89  10 135];
vivid_purple       = [148   0 230];

cb_colors   = [deep_pink; pompadour ; persian_rose ; persian_green;...
                      aqua_marine ; carissma ; deep_sea; dark_olive ;...
                      aqua ; royal_blue ; indigo ; vivid_purple ]./255;
                  
set(h_mean, {'color'}, num2cell(cb_colors,2));
%}

xticklabels(cfg.labels);
ylabel([cfg.ylabel_txt]);


%ylim([min(ind_by_group,[],'all')*.8, max(ind_by_group,[],'all')*1.2]);

grid on
box off


end