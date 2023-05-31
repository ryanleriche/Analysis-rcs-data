function ft_form_TD_struct = rcs_TD_peri_survey(cfg,  rcs_pp_dir, pt_side_id, par_db, REDcap)
%{




%}
par_db_RCS0XX   = par_db.(pt_side_id);
redcap          = REDcap.(pt_side_id(1:end-1));

for i_ss  = 1 : height(par_db_RCS0XX)
    i_rcap_oi       = find(...
                            isbetween(redcap.time, ...
                                par_db_RCS0XX.timeStart(i_ss),...
                                par_db_RCS0XX.timeStop(i_ss)));
    if ~isempty(i_rcap_oi)
        par_db_RCS0XX.rcap_indices{i_ss}             = i_rcap_oi;
        par_db_RCS0XX.rcap_latency_timeStart{i_ss}    = redcap.time(i_rcap_oi) - par_db_RCS0XX.timeStart(i_ss);
        par_db_RCS0XX.rcap_N{i_ss}                    = length(i_rcap_oi);
    end
end

N_rcap_assigned = length(unique(vertcat(par_db_RCS0XX.rcap_indices{:})));
N_rcap_possible = sum(isbetween(redcap.time, ...
                            par_db_RCS0XX.timeStart(1), par_db_RCS0XX.timeStop(end)));


fprintf('%s | %g/%g, %.1f%%, of REDcap surveys w/n streaming sessions (**pre-session** rejection)\n',...
    pt_side_id, N_rcap_assigned, N_rcap_possible, N_rcap_assigned/N_rcap_possible*100)

%% save parsed database of interest, the power spectrum per session, frequencies, and channel names
save_dir = fullfile(rcs_pp_dir,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset, cfg.pp_fft,'/');

if ~isfolder(save_dir);     mkdir(save_dir);    end

ft_data_file = [save_dir, pt_side_id, '_ft_form_TD_struct.mat'];

if isfile(ft_data_file) && ~cfg.ignoreold_td_parsing

    fprintf('%s | loading previously ran TD parsing\n', pt_side_id)
    load(ft_data_file);

    return

end
%%
cfg.r_cap_cent_toi   = [duration('0:10:00'), duration('0:00:00')];
cfg.epoch_dur        = duration('0:00:30');
cfg.epoch_nearest_to = 'end';
cfg.min_sess_dur     = cfg.epoch_dur;
%%

%%
dir_ss_pt_side  = fullfile(rcs_pp_dir,'streaming_sessions', pt_side_id, cfg.pp_RCS_TD_subset, '/');

if ~isfolder(dir_ss_pt_side);      error('%s | path to TD data was not found', pt_side_id);    end


save_dir_struct  = dir([dir_ss_pt_side, '*Session*']);

par_db_oi  = par_db_RCS0XX(...
                    ge(par_db_RCS0XX.duration, cfg.min_sess_dur) &...
                    ismember(par_db_RCS0XX.sess_name, {save_dir_struct.name}) &...
                    par_db_RCS0XX.TDsampleRates == 250 ...
                    ,:);

%% circumvent FieldTrip reading functions by matching output of their prepreprocessing
% https://www.fieldtriptoolbox.org/faq/how_can_i_import_my_own_dataformat/
par_db_oi_out     = table;
% FieldTrip formatted time-domain structure
ft_form_TD_struct = struct;
ft_form_TD_struct.dimord = 'chan_time';

cnt = 1;

for i_ss = 1 : height(par_db_oi)
    
    i_rcap_wn_ss = par_db_oi.rcap_indices{i_ss};

    % do NOT parse time-domain if latency < cfg.epoch_dur (accounts for empty cases)
    if all(le(par_db_oi.rcap_latency_timeStart{i_ss}, cfg.epoch_dur))
        continue
    end

    tmp          = load([dir_ss_pt_side, par_db_oi.sess_name{i_ss},'/', 'rcs_streamed.mat']);
    rcs_streamed = tmp.rcs_streamed;

    clear tmp

    for i_rcap = 1 : length(i_rcap_wn_ss)
    
        [rcs_streamed_oi, par_db_oi_out(cnt,:), Fsample_TD] ...
            ...
            =  return_TD_rcs_chunk(...
            ...
        cfg, i_rcap_wn_ss,i_rcap, redcap, rcs_streamed, par_db_oi(i_ss,:));


        if ~isempty(rcs_streamed_oi)
            if cnt == 1
                ch_names = parse_chan_names(pt_side_id, par_db_oi_out(cnt,:))';

                ft_form_TD_struct.label   = ch_names;    % cell-array containing strings, Nchan*1
                ft_form_TD_struct.fsample = Fsample_TD;  % sampling frequency in Hz, single number
               
            end
        
            ft_form_TD_struct.trial{1, cnt}  ...                   % cell-array containing a data matrix for each
                = rcs_streamed_oi{:, compose('Ch_TD%g',0:3)}';     % trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
        
            ft_form_TD_struct.time{1, cnt}  ...                    % cell-array containing a time axis for each
                 = seconds(...
                  rcs_streamed_oi.TimeInPST - ...
                  redcap.time(i_rcap_wn_ss(i_rcap)))';  % trial (1*Ntrial), each time axis is a 1*Nsamples vector                       
        
            cnt = cnt+1;
        end
    end
    clear rcs_streamed
end
%% save outputs
% remove streaming sessions w/o REDcap survey
i_w_sur            = find(~cellfun(@isempty, par_db_oi_out.sess_name));
par_db_oi_out      = par_db_oi_out(i_w_sur,:);

ft_form_TD_struct.rcs = struct;
ft_form_TD_struct.rcs.par_db = par_db_oi_out;
ft_form_TD_struct.rcs.td_cfg = cfg;

save([save_dir, pt_side_id, '_ft_form_TD_struct'], 'ft_form_TD_struct', '-v7.3');
save([save_dir, pt_side_id, '_ch_names'], 'ch_names', '-v7.3');
writetable(par_db_oi_out, [save_dir,pt_side_id, '_parsed_db_oi.xlsx']);


%% plot summary of session durations based on epoching criteria
figure
h = cdfplot(par_db_oi_out.epoch_durationInSecs/60);

title(sprintf('%s | Streaming Sessions', pt_side_id));


epoching_txt = sprintf(['min(session duration), %g min             \n',...
                        'epoch, [%g, %g] min (survey at 0)          \n',...
                         'min(epoch duration), %g min              \n',...
                     ], ...
                     minutes([cfg.min_sess_dur, -cfg.r_cap_cent_toi(1), cfg.epoch_dur])...
                        );

hand_txt = TextLocation(epoching_txt);     hand_txt.FontSize = 12;

h.YData  = 1 -h.YData;


yyaxis right; ylabel('Empirical CDF Probability'); set(gca, 'FontSize', 16)

h.YData   = h.YData * height(par_db_oi_out);

yyaxis left; ylabel('N sessions'); 

xlabel('Session Duration (minutes)');


i_plt  = round(linspace(2, length(h.YData),10));

c = cellfun(@(x) sprintf('%s sessions', num2str(x)), num2cell(h.YData), 'UniformOutput', false);
text(h.XData(i_plt), h.YData(i_plt), c(i_plt), 'Fontsize', 12);

set(gca, 'FontSize', 16);
%%%
anal_dir = fullfile(cfg.anal_dir, cfg.pp_RCS_TD_subset, pt_side_id);

if ~isfolder(anal_dir);       mkdir(anal_dir);      end

exportgraphics(gcf, fullfile(anal_dir, 'session_TD_duration.png'));
%%
N_rcap_assigned = length(unique(par_db_oi_out.rcap_indices));
N_rcap_possible = sum(isbetween(redcap.time, ...
                            par_db_oi.timeStart(1), par_db_oi.timeStop(end)));

fprintf('%s | %g/%g, %.1f%%, of REDcap surveys w/n streaming sessions (**post-session** rejection)\n',...
    pt_side_id, N_rcap_assigned, N_rcap_possible, N_rcap_assigned/N_rcap_possible*100)

end