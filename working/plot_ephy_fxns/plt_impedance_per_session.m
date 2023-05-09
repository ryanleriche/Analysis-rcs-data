function plt_impedance_per_session(cfg, pt_side_id, db)
%%
%{
Per session, return the case to contact impedance as four traces over weeks

%}


%%
eventLog_jsons = vertcat(db.(pt_side_id).eventLogTable{:});

%u_event_names  = unique(eventLog_jsons.EventType);


i_lead_int     = strcmp(eventLog_jsons.EventType, 'Lead Integrity');
tmp_tbl        = eventLog_jsons(i_lead_int, :);


lead_int_tbl = table;

lead_int_tbl.sess_name = tmp_tbl.SessionId;
lead_int_tbl.time      =  datetime(tmp_tbl.HostUnixTime / 1000,...
                                    'ConvertFrom','posixTime','TimeZone','America/Los_Angeles',...
                                    'Format','dd-MMM-yyyy HH:mm:ss.SSS');

i_delimiter            = cellfun(@(x) {regexp(x, '---')}, tmp_tbl.EventSubType);

lead_int_tbl.impedanceBtwn      = cellfun(@(x,y) x(1:y-2), ...
                                     tmp_tbl.EventSubType, i_delimiter, ...
                                     'UniformOut', false);

tmp_int                          = cellfun(@(x,y) str2double(x(y+4:end)), ...
                                        tmp_tbl.EventSubType, i_delimiter,...
                                        'UniformOut', false);

lead_int_tbl.impedanceInKiloOhms = [tmp_int{:}]';

mono_int_tbl = lead_int_tbl(contains(lead_int_tbl.impedanceBtwn, {'16'}), :);

u_mono_int   = unique(mono_int_tbl.impedanceBtwn);
%%%

figure('Units', 'Inches', 'Position', [0, 0, 15, 7]);
            
c_map   = brewermap(NaN, 'Set2');

%tiledlayout(length(i_pain_lbls),1, 'TileSpacing','compact')

for j=1:length(u_mono_int)

    plt_tbl = mono_int_tbl(strcmp(mono_int_tbl.impedanceBtwn, u_mono_int{j}), :);


    plot(plt_tbl.time, plt_tbl.impedanceInKiloOhms, ...
        'o-','Color', c_map(j,:), 'LineWidth', 1);
    
    hold on
    

end

legend(u_mono_int, 'Location', 'northoutside', 'Orientation', 'Horizontal');     

ylabel('Impedance (kΩ)')

grid on; grid minor; legend box off

title(pt_side_id)

set(gca, 'FontSize', 16, 'TickLength', [0,0]);

%%%
cfg.proc_dir = [cfg.proc_dir, '/impedance_per_session/'];

if ~isfolder(cfg.proc_dir);    mkdir(cfg.proc_dir);     end

exportgraphics(gcf, [cfg.proc_dir, pt_side_id, '_per_session_impedance.png'])

end