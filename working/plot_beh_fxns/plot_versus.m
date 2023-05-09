function plot_versus(cfg, REDcap)

redcap  = REDcap.(cfg.pt_id);
    

[redcap, date_range] = date_parser(cfg, redcap);


redcap.MPQaff       = sum([redcap.MPQsickening, redcap.MPQfearful, redcap.MPQcruel,  redcap.MPQtiring],2,'omitnan');
redcap.MPQsom       = redcap.MPQtotal - redcap.MPQaff;


% only KEEP if all VAS reports do NOT equal 50

if strcmp(cfg.pt_id, 'RCS06')
    i_trim_VAS_any           = ~(redcap.unpleasantVAS == 50 | redcap.painVAS == 50 | redcap.worstVAS == 50 ...
                                   | redcap.npVAS == 50 | redcap.nocVAS == 50);
end

% only KEEP reports where 3 VAS reports do NOT equal 50, and in cases where
% at least one VAS report equals 50 keep it ONLY if NRS equals 5 (i.e., a true
% neutral report)

i_trim_VAS_all    = ~(redcap.unpleasantVAS == 50 & redcap.painVAS == 50 & redcap.worstVAS == 50)...
                    &...
                        ~(...
                        (redcap.unpleasantVAS == 50 | redcap.painVAS == 50 | redcap.worstVAS == 50)...
                        & redcap.mayoNRS ~= 5);


i_missing          = isnan(redcap.mayoNRS) | isnan(redcap.painVAS);


i_mpq_answered     = ~((redcap.painVAS >=40 | redcap.mayoNRS >=4) & redcap.MPQtotal ==0);

i_kept             = i_trim_VAS_all & ~i_missing & i_mpq_answered;

RCSXX_kept         = redcap(i_kept, :);


%% z-score w/n pain metrics
pain_met_names = redcap.Properties.VariableNames;

% if unpleasantVAS, painVAS, and worstVAS are ALL 50s, ignore those reports

z_RCSXX      = RCSXX_kept;

zscor_xnan   = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));


for i = 2 : length(pain_met_names)

   z_RCSXX.(pain_met_names{i}) = zscor_xnan(z_RCSXX.(pain_met_names{i}));

end


%%
ds  =    datestr(date_range,'dd-mmm-yyyy');

figure('Units', 'Inches', 'Position', [0, 0, 25, 25]);
sgtitle([cfg.pt_id, newline, ds(1,:) ' to ' ds(2,:)], 'Fontsize',24);
hold on

switch cfg.pt_id
    case 'RCS06'
         dropped_cols = {'time', 'MPQthrobbing', 'MPQshooting', 'MPQsharp', ...
                  'MPQcramping', 'MPQgnawing', 'MPQhot_burning', ...
                  'MPQaching', 'MPQheavy', 'MPQtender', 'MPQsplitting',...
                  'MPQtiring', 'MPQsickening', 'MPQcruel',...
                  'MPQstabbing', 'MPQfearful',...
                  ...
                  'MPQaff', 'MPQsom',...
                  'worstNRS', 'worstVAS'};

        parsi_pain_tbl = removevars(RCSXX_kept, dropped_cols);

    case 'RCS07'
         dropped_cols = {'time', 'MPQthrobbing', 'MPQshooting', 'MPQsharp', ...
                  'MPQcramping', 'MPQgnawing', 'MPQhot_burning', ...
                  'MPQaching', 'MPQheavy', 'MPQtender', 'MPQsplitting',...
                  'MPQtiring', 'MPQsickening', 'MPQcruel',...
                  'MPQstabbing', 'MPQfearful',...
                  ...
                  'MPQaff', 'MPQsom',...
                  'worstNRS', 'worstVAS'};

        parsi_pain_tbl = removevars(RCSXX_kept, dropped_cols);

    otherwise
        dropped_cols = {'time', 'MPQthrobbing', 'MPQshooting', 'MPQsharp', ...
                  'MPQcramping', 'MPQgnawing', 'MPQhot_burning', ...
                  'MPQaching', 'MPQheavy', 'MPQtender', 'MPQsplitting',...
                  'MPQtiring', 'MPQsickening', 'MPQcruel',...
                  'MPQstabbing', 'MPQfearful',...
                  ...
                  'MPQaff', 'MPQsom'};


    parsi_pain_tbl = removevars(RCSXX_kept, dropped_cols);
    parsi_pain_tbl = movevars(parsi_pain_tbl, 'worstNRS', 'Before', 'painVAS');
    
end


[~,AX,~,~,~] ...
    ...
    = plotmatrix(parsi_pain_tbl.Variables);
 

% use longest non-NaN rows of pairwise columns to find correlation btwn all
% columns
[corr_coefs, corr_pvals] = ...
    ...
    corr(parsi_pain_tbl.Variables, 'Rows', 'pairwise');

% resize every axis
for i = 1 : length(AX)
    pain_name = parsi_pain_tbl.Properties.VariableNames{i};

    AX(i, 1).YLabel.String   = pain_name;
    AX(i, 1).YLabel.FontSize = 14;

    AX(length(AX), i).XLabel.String = pain_name;
    AX(length(AX), i).XLabel.FontSize = 14;

    % set axes
    switch cfg.pt_id
        case  'RCS06'
            AX(i, 1).XLim = [-1, 11];       AX(i, 2).XLim = [-1, 11];
            AX(i, 3).XLim = [-1, 11];       AX(i, 4).XLim = [-10, 110];
            AX(i, 5).XLim = [-10, 110];     AX(i, 6).XLim = [-10, 110];
            AX(i, 7).XLim = [-10, 110];     AX(i, 8).XLim = [-5, 45];
    
            % same for columns
            AX(1, i).YLim = [-1, 11];       AX(2, i).YLim = [-1, 11];
            AX(3, i).YLim = [-1, 11];       AX(4, i).YLim = [-10, 110];
            AX(5, i).YLim = [-10, 110];     AX(6, i).YLim = [-10, 110];
            AX(7, i).YLim = [-10, 110];     AX(8, i).YLim = [-5, 45];


         case  'RCS07'
            AX(i, 1).XLim = [-1, 11];       AX(i, 2).XLim = [-1, 11];
            AX(i, 3).XLim = [-1, 11];       AX(i, 4).XLim = [-1, 11];
            AX(i, 5).XLim = [-1, 11];       AX(i, 6).XLim = [-10, 110];
            AX(i, 7).XLim = [-10, 110];     AX(i, 8).XLim = [-10, 110];

            AX(i, 9).YLim = [-5, 45];
    
            % same for columns
            AX(1, i).YLim = [-1, 11];       AX(2, i).YLim = [-1, 11];
            AX(3, i).YLim = [-1, 11];       AX(4, i).YLim = [-1, 11];
            AX(5, i).YLim = [-1, 11];       AX(6, i).YLim = [-10, 110];
            AX(7, i).YLim = [-10, 110];     AX(8, i).YLim = [-10, 110];

            AX(9, i).YLim = [-5, 45];
    
        otherwise
            AX(i, 1).XLim = [-1, 11];       AX(i, 2).XLim = [-1, 11];
            AX(i, 3).XLim = [-10, 110];     AX(i, 4).XLim = [-10, 110];
            AX(i, 5).XLim = [-10, 110];     AX(i, 6).XLim = [-5, 45];
    
            % same for columns
            AX(1, i).YLim = [-1, 11];       AX(2, i).YLim = [-1, 11];
            AX(3, i).YLim = [-10, 110];     AX(4, i).YLim = [-10, 110];
            AX(5, i).YLim = [-10, 110];     AX(6, i).YLim = [-5, 45];
    end
end

% add correlation coefficient and p-value
for i = 1 : length(AX)
    for j = 1 : length(AX)

        if corr_coefs(i,j) >= 0.8

            AX(i, j).Title.Color = 'r';

        end

        AX(i, j).Title.String    = ['rho = ', sprintf('%.2f',corr_coefs(i,j)),...
                                     newline,...
                                     'pval = ', sprintf('%.3g',corr_pvals(i,j))];

        AX(i, j).Title.HorizontalAlignment = 'left';
        AX(i, j).Title.FontWeight = "normal";

        % slightly different formatting since RCS07 has more pain metrics
        switch cfg.pt_id
            case 'RCS07'
                AX(i, j).Title.Position(2) = AX(i, j).Title.Position(2) * .6;
            otherwise
                AX(i, j).Title.Position(2) = AX(i, j).Title.Position(2) * .7;
        end

         AX(i, j).Title.Position(1) = AX(i, j).Title.Position(1) * 0.05;
    end
end

i_keep = tril(ones(size(AX)));

for i =1 : size(AX, 1)
    for j = 1: size(AX, 2)

        if i_keep(i,j) ~= 1
            delete(AX(i,j))

        end
    end
end

format_plot();

saveas(gcf, [cd, '/plot_beh/figs/beh_only/', cfg.pt_id, '/plot_versus.png']);



function format_plot()  

    set(gca,'fontSize',14, 'TickLength', [0 0]); 
    grid on;    grid MINOR;      box off; 
    
end
end