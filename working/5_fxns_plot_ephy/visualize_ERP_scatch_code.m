
%{
* import network-mapping session of interest
    * based on date
    * see if stimLog can be plotted on top of ft_databrowser
        * something like vertical lines w/ stim parameter changes
        * Does this fit within the FieldTrip events motiff?

* see if TMS-evoked potential tutorial in FieldTrip works for sEEG stim-ERPs





%}

par_db_oi   = par_db_out.(pt_sides{i});
nw_dates_oi = {'23-Jun-2023'};



dir_ss_pt_side  = fullfile(dirs.rcs_preproc,'streaming_sessions', pt_sides{i}, sub_cfg.pp_RCS_TD_subset);
nw_dates_oi     = datetime(nw_dates_oi, 'TimeZone','America/Los_Angeles');




[~, i_ss_oi] = min(abs(nw_dates_oi-par_db_oi.timeStart));

nw_ss_par_db = par_db_oi(i_ss_oi, :);


tmp          = load(...
                    fullfile(dir_ss_pt_side, nw_ss_par_db.sess_name{1}, 'rcs_streamed.mat')...
                    );
rcs_streamed = tmp.rcs_streamed;

clear tmp


stim_log          = readtable(fullfile(dir_ss_pt_side, nw_ss_par_db.sess_name{1}, 'session_stim_log.xls'));


%unique(rcs_streamed.TD_samplerate(~isnan(rcs_streamed.TD_samplerate)))
% FieldTrip formatted time-domain structure
ft_form_TD_struct        = struct;
ft_form_TD_struct.dimord = 'chan_time';


ch_names                  = nw_ss_par_db{1, compose("Ch%g_chanOut", 0:3)}';
ft_form_TD_struct.label   = ch_names;                    % cell-array containing strings, Nchan*1
ft_form_TD_struct.fsample = nw_ss_par_db.TDsampleRates;  % sampling frequency in Hz, single number  


ft_form_TD_struct.trial{1, 1}  ...                       % cell-array containing a data matrix for each
    = rcs_streamed{:, compose('Ch_TD%g',0:3)}';          % trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix

ft_form_TD_struct.time{1, 1}  ...                        % cell-array containing a time axis for each
     = seconds(...
      rcs_streamed.TimeInPST - ...
      rcs_streamed.TimeInPST(1))';  % trial (1*Ntrial), each time axis is a 1*Nsamples vector                       
%%

stim_log.duration = [diff(stim_log.stimLogTimeInPT); NaN];
%%
cfg = [];
cfg.continuous = 'yes';


ft_databrowser(cfg, ft_form_TD_struct)




