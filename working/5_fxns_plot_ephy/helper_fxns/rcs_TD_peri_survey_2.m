function ft_form_TD_struct = rcs_TD_peri_survey_2(cfg,  rcs_pp_dir, pt_side_id, vargin)
%{

cfg                 = cfg_rcs;
pt_side_id          = pt_sides{i};
rcs_pp_dir          = dirs.rcs_preproc;

par_db              = par_db_out;



%}
%%
set(0,'DefaultFigureVisible','off')

cfg.r_cap_cent_toi   = [duration('0:20:00'), duration('0:00:00')];
cfg.min_sess_dur     = duration('0:05:00');

cfg.epoch_dur        = duration('0:00:30');


%% save parsed database of interest, the power spectrum per session, frequencies, and channel names
save_dir = fullfile(rcs_pp_dir,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset);

% make folder WITHOUT pt side id
if ~isfolder(save_dir);     mkdir(save_dir);        end

% use pt side id as prefix for all output files from here on forward
save_dir = fullfile(save_dir, pt_side_id);


ft_data_file = fullfile([save_dir, '_ft_form_TD_struct.mat']);



if isfile(ft_data_file) && ~cfg.ignoreold_td_parsing

    fprintf('%s | loading previously ran TD parsing\n', pt_side_id)
    load(ft_data_file);

    load(fullfile([save_dir, '_ch_names.mat']))

    par_db_oi_out = ft_form_TD_struct.rcs.par_db;

else

    
    par_db_RCS0XX   = vargin{1}.(pt_side_id);
    redcap          = vargin{2}.(pt_side_id(1:end-1));
    
    
    for i_ss  = 1 : height(par_db_RCS0XX)
        i_rcap_oi       = find(...
                                isbetween(redcap.time, ...
                                    par_db_RCS0XX.timeStart(i_ss),...
                                    par_db_RCS0XX.timeStop(i_ss)));
        if ~isempty(i_rcap_oi)
            par_db_RCS0XX.rcap_indices{i_ss}              = i_rcap_oi;
            par_db_RCS0XX.rcap_latency_timeStart{i_ss}    = redcap.time(i_rcap_oi) - par_db_RCS0XX.timeStart(i_ss);
            par_db_RCS0XX.rcap_N{i_ss}                    = length(i_rcap_oi);
        end
    end
    
    N_rcap_assigned = length(unique(vertcat(par_db_RCS0XX.rcap_indices{:})));
    N_rcap_possible = sum(isbetween(redcap.time, ...
                                par_db_RCS0XX.timeStart(1), par_db_RCS0XX.timeStop(end)));
    
    fprintf('%s | %g/%g, %.1f%%, of REDcap surveys w/n streaming sessions (**pre-session** rejection)\n',...
        pt_side_id, N_rcap_assigned, N_rcap_possible, N_rcap_assigned/N_rcap_possible*100)

%%


dir_ss_pt_side  = fullfile(rcs_pp_dir,'streaming_sessions', pt_side_id, cfg.pp_RCS_TD_subset, '/');

if ~isfolder(dir_ss_pt_side);      error('%s | path to TD data was not found', pt_side_id);    end

save_dir_struct  = dir([dir_ss_pt_side, '*Session*']);

par_db_oi  = par_db_RCS0XX(...
                    ge(par_db_RCS0XX.duration, cfg.min_sess_dur) &...
                    ismember(par_db_RCS0XX.sess_name, {save_dir_struct.name}) &...
                    par_db_RCS0XX.TDsampleRates == 250 & ...
                    ~cellfun(@isempty, par_db_RCS0XX.rcap_indices)...
                    ,:);


    %% circumvent FieldTrip reading functions by matching output of their prepreprocessing
    % https://www.fieldtriptoolbox.org/faq/how_can_i_import_my_own_dataformat/
    par_db_oi_out     = table;
    
    % FieldTrip formatted time-domain structure
    ft_form_TD_struct        = struct;
    ft_form_TD_struct.dimord = 'chan_time';
    
    cnt = 1;
    
    for i_ss = 1 : height(par_db_oi)
        
        i_rcap_wn_ss              = par_db_oi.rcap_indices{i_ss};
        i_keep_rcap               = ge(par_db_oi.rcap_latency_timeStart{i_ss}, cfg.epoch_dur);
    
        i_rcap_wn_ss              = i_rcap_wn_ss(i_keep_rcap);
    
        %%% --> ONLY parse time-domain if enough (cfg.epoch_dur) TD data before
        %%% REDcap survey
        if isempty(i_rcap_wn_ss);       continue;               end
    
        tmp          = load([dir_ss_pt_side, par_db_oi.sess_name{i_ss},'/', 'rcs_streamed.mat']);
        rcs_streamed = tmp.rcs_streamed;
    
        clear tmp
    
        for i_rcap = 1 : length(i_rcap_wn_ss)
        
            [rcs_streamed_oi, par_db_oi_out(cnt,:), Fsample_TD] ...
                ...
                =  return_TD_rcs_chunk_2(...
                ...
            cfg, i_rcap_wn_ss, i_rcap, redcap, rcs_streamed, par_db_oi(i_ss,:));
    
    
            if isempty(rcs_streamed_oi);     continue;        end
    
            if cnt == 1
                ch_names                  = par_db_oi_out{cnt, compose("Ch%g_chanFullStr", 0:3)}';
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
        clear rcs_streamed
    end
    %% save outputs
    % remove streaming sessions w/o REDcap survey
    i_w_sur                   = find(cellfun(@isempty, par_db_oi_out.sess_name));
    par_db_oi_out(i_w_sur,:) = [];
    
    ft_form_TD_struct.rcs         = struct;
    ft_form_TD_struct.rcs.par_db  = par_db_oi_out;
    ft_form_TD_struct.rcs.td_cfg  = cfg;
    
    save([save_dir, '_ft_form_TD_struct'], 'ft_form_TD_struct', '-v7.3');
    save([save_dir,  '_ch_names'], 'ch_names', '-v7.3');
    
    writetable(par_db_oi_out, [save_dir, '_parsed_db_oi.xlsx']);
    %%%

    N_rcap_assigned = length(unique(par_db_oi_out.rcap_indices));
    N_rcap_possible = sum(isbetween(redcap.time, ...
                            par_db_oi.timeStart(1), par_db_oi.timeStop(end)));
    
    fprintf('%s | %g/%g, %.1f%%, of REDcap surveys w/n streaming sessions (**post-session** rejection)\n',...
    pt_side_id, N_rcap_assigned, N_rcap_possible, N_rcap_assigned/N_rcap_possible*100)

end
%% plot summary of session durations based on epoching criteria

try 
    par_db_oi_out = ft_form_TD_struct.rcs.par_db;

    epoch_min_dur    = par_db_oi_out.epoch_durationInSecs/60;
    epoch_min_dur_TD = (par_db_oi_out.percent_TD_streamed /100).* epoch_min_dur;
    
    
    longest_chunk = epoch_min_dur;
    
    for j = 1 :length(longest_chunk)
    
        longest_chunk(j) = max(par_db_oi_out.td_chunk{j}.length(...
                            strcmp(par_db_oi_out.td_chunk{j}.label, 'TD_chunk')));
    
    
    end
    
    Fsample_TD = unique(par_db_oi_out.TDsampleRates);
    
    longest_chunk = longest_chunk / (60 *Fsample_TD);
    
    figure
        
        h1        = cdfplot(epoch_min_dur); hold on;
        h2        = cdfplot(epoch_min_dur_TD); 
        h3        = cdfplot(longest_chunk);
    
        h1.YData  = 1 -h1.YData;
        h2.YData  = 1 -h2.YData;
        h3.YData  = 1 - h3.YData;
    
        
        title(sprintf('%s | avaiable TD from each streaming session', pt_side_id));
        legend({'TD collected'; 'Cumulative TD w/o packet-loss'; 'Longest TD chunk'}, ...
            'Location','northoutside', 'Orientation','horizontal', 'NumColumns',2)
        
        ylabel('N sessions');
        
        yyaxis right; 
        
            ylabel('Empirical CDF Probability'); set(gca, 'FontSize', 14, 'YColor', 'k')
            
            h1.YData   = h1.YData * height(par_db_oi_out);
            h2.YData   = h2.YData * height(par_db_oi_out);
            h3.YData   = h3.YData * height(par_db_oi_out);
        
        
            epoching_txt = sprintf(['min(session duration), %g min             \n',...
                                    'epoch, [%g, %g] min (survey at 0)          \n',...
                                    ], ...
                                    minutes([cfg.min_sess_dur, -cfg.r_cap_cent_toi])...
                                    );
            
            text(1, max(h1.YData) * .2, epoching_txt, 'FontSize', 14);
           
            xlabel('Minutes prior to REDcap survey');
           
            set(gca, 'FontSize', 14); grid on; grid minor;
    %%%
    anal_dir = fullfile(cfg.anal_dir, cfg.pp_RCS_TD_subset, pt_side_id);
    
    if ~isfolder(anal_dir);       mkdir(anal_dir);      end
    
    exportgraphics(gcf, fullfile(anal_dir, 'session_TD_duration.png'));
catch
end
%%

end