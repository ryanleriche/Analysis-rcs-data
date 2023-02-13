function db_RCSXX = calc_INS_lat_per_session(cfg, db_RCSXX)
%% INS to API time latency per session

% cfg                = [];
% cfg.pt_id_side     = 'RCS02R';
% 
% cfg.save_dir       = [github_dir, 'Analysis-rcs-data/working/plot_ephy/aDBS_offline_sessions/'];
% 
% 
% db_RCSXX           = db.(cfg.pt_id_side);

db_RCSXX.INS_API_latency = repmat({''}, height(db_RCSXX),1);


for i=1:height(db_RCSXX)

    if ~isempty(db_RCSXX.path{i})

        TimeSync_filename = findFilesBVQX(db_RCSXX.path{i},'TimeSync.json');
        try
            if ~isempty(TimeSync_filename)
        
                db_RCSXX.INS_API_latency{i}  = INS_API_TimeSync(TimeSync_filename{1});
        
            else
        
                db_RCSXX.INS_API_latency{i}  = 'No entry';
            end
    
        catch
            db_RCSXX.INS_API_latency{i} = 'Error';
    
        end
    else
        db_RCSXX.INS_API_latency{i}  = 'No sess path';
    end
end

%%
i_TimeSync          = find(cellfun(@isstruct, db_RCSXX.INS_API_latency));


INS_API_latency_tbl = [db_RCSXX(i_TimeSync, {'timeStart','path'})...
                        ...
                        struct2table([db_RCSXX.INS_API_latency{i_TimeSync}])];

INS_API_latency_tbl(...
    cellfun(@isempty,INS_API_latency_tbl.timeStart),:) = [];

INS_API_latency_tbl.timeStart =     cellfun(@(x) x(1), INS_API_latency_tbl.timeStart);

%%

INS_API_lat_fig = figure('Units', 'Inches', 'Position', [0, 0, 12, 5]);

stairs(INS_API_latency_tbl.timeStart, INS_API_latency_tbl.mean);

ylabel(['INS Firmware "timestamp"', newline,...
       '–',newline,...
       'API "PacketGenTime"'])

INS_API_latency_tbl.timeStart.TimeZone = 'America/Los_Angeles';


dst_range = datetime({'2020-Mar-08 2:00 AM', '2020-Nov-01 2:00 AM';...
                      '2021-Mar-14 2:00 AM', '2021-Nov-07 2:00 AM';...
                      '2022-Mar-13 2:00 AM', '2022-Nov-06 2:00 AM';...
                      '2023-Mar-12 2:00 AM', '2023-Nov-05 2:00 AM'; ...
                      '2024-Mar-10 2:00 AM', '2024-Nov-03 2:00 AM'; ...
                      '2025-Mar-09 2:00 AM', '2025-Nov-02 2:00 AM'});

dst_range.TimeZone = 'America/Los_Angeles';

y_limits_patch     = [INS_API_lat_fig.CurrentAxes.YLim(1), INS_API_lat_fig.CurrentAxes.YLim(2),...
                      INS_API_lat_fig.CurrentAxes.YLim(2),INS_API_lat_fig.CurrentAxes.YLim(1)];

for i=1 : size(dst_range,1)

    patch(sort([dst_range(i,:),dst_range(i,:)]), y_limits_patch, [0.7, 0.7, 0.7],...
        'FaceAlpha', 0.5, 'EdgeColor', 'none')

end


legend({'Session Start Time', 'Daylight Savings Time'});

title([cfg.pt_id_side, ' | ', 'INS time does NOT account for DST and varies over time'])

set(gca, 'XLim', [min(INS_API_latency_tbl.timeStart), max(INS_API_latency_tbl.timeStart)],...
    'FontSize', 14);

exportgraphics(gcf, [cfg.save_dir, cfg.pt_id_side, '_INS_API_latency.png'])

end
