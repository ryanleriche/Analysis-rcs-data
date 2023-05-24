function ft_form_TD_struct = rcs_TD_peri_survey(cfg,  rcs_pp_dir, pt_side_id, par_db, REDcap)
%{


%}
%% save parsed database of interest, the power spectrum per session, frequencies, and channel names
save_dir = fullfile(rcs_pp_dir,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset,cfg.pp_fft,'/');

if ~isfolder(save_dir);     mkdir(save_dir);    end

ft_data_file = [save_dir, pt_side_id, '_ft_form_TD_struct.mat'];

if isfile(ft_data_file) && ~cfg.ignoreold_td_parsing

    fprintf('%s | loading previously ran TD parsing\n', pt_side_id)
    load(ft_data_file);

    return

end
%%
cfg.min_sess_dur    = duration('00:05:00');

%%%
cfg.fftSize         = 1024;             % N samples to use in FFT
cfg.windowLoad      = '100% Hann';      % 100 percent hanning window
cfg.fftinterval     = 500;              % how often in ms a FFT should be calculated


cfg.r_cap_cent_toi = [duration('0:02:00'), duration('0:00:00')];

cfg.epoch_dur      = duration('0:00:30');

cfg.epoch_nearest_to = 'end';
%%

%%
dir_ss_pt_side  = fullfile(rcs_pp_dir,'streaming_sessions', pt_side_id, cfg.pp_RCS_TD_subset, '/');

if ~isfolder(dir_ss_pt_side);      error('%s | path to TD data was not found', pt_side_id);    end


db_RCSXXX        = par_db.(pt_side_id);
save_dir_struct  = dir([dir_ss_pt_side, '*Session*']);

par_db_oi  = db_RCSXXX(...
                    ge(db_RCSXXX.duration, cfg.min_sess_dur) &...
                    ismember(db_RCSXXX.sess_name, {save_dir_struct.name}) &...
                    db_RCSXXX.TDsampleRates == 250 ...
                    ,:);


redcap        = REDcap.(pt_side_id(1:end-1));

%% circumvent FieldTrip reading functions by matching output of their prepreprocessing
% https://www.fieldtriptoolbox.org/faq/how_can_i_import_my_own_dataformat/

par_db_oi_out = table;
% FieldTrip formatted time-domain structure
ft_form_TD_struct = struct;
ft_form_TD_struct.dimord = 'chan_time';

cnt = 1;

for i_sess = 1 : height(par_db_oi)
    
    tmp          = load([dir_ss_pt_side, par_db_oi.sess_name{i_sess},'/', 'rcs_streamed.mat']);
    rcs_streamed = tmp.rcs_streamed;

    clear tmp

    [rcs_streamed_oi, par_db_oi_out(i_sess,:), Fsample_TD] ...
        ...
        =  pull_continuous_TD_rcs_chunk(...
        ...
    cfg, redcap, rcs_streamed, par_db_oi(i_sess,:));

    if ~isempty(rcs_streamed_oi)
        if cnt == 1
            ft_form_TD_struct.fsample = Fsample_TD;              % sampling frequency in Hz, single number
            ft_form_TD_struct.label   = compose('Ch_TD%g',0:3)'; % cell-array containing strings, Nchan*1
        end
    
        ft_form_TD_struct.trial{1, cnt}  ...                   % cell-array containing a data matrix for each
            = rcs_streamed_oi{:, ft_form_TD_struct.label}';    % trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
    
        ft_form_TD_struct.time{1, cnt}  ...                    % cell-array containing a time axis for each
             = seconds(...
                 rcs_streamed_oi.TimeInPST -...
                 rcs_streamed_oi.TimeInPST(1))';                % trial (1*Ntrial), each time axis is a 1*Nsamples vector
    
        ft_form_TD_struct.sampleinfo(cnt, 1:2) = [ft_form_TD_struct.time{1, cnt}(1), ft_form_TD_struct.time{1, cnt}(end)];
        cnt = cnt+1;
    end
    clear rcs_streamed tmp
end
%%
% remove streaming sessions w/o REDcap survey
i_w_sur            = find(~cellfun(@isempty, par_db_oi_out.sess_name));
par_db_oi_out      = par_db_oi_out(i_w_sur,:);

ch_names = parse_chan_names(pt_side_id, par_db_oi_out);


save([save_dir, pt_side_id, '_ft_form_TD_struct'], 'ft_form_TD_struct', '-v7.3');
save([save_dir, pt_side_id, '_ch_names'], 'ch_names', '-v7.3');
writetable(par_db_oi_out, [save_dir,pt_side_id, '_parsed_db_oi.xlsx']);



end