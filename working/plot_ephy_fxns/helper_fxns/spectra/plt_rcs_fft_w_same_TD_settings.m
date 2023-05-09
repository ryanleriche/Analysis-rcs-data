function plt_rcs_fft_w_same_TD_settings(cfg, pt_side_id, dirs)

%%% load parsed database of interest, the power spectrum per session, and
%%% the frequencies, and channel names
load_dir = [dirs.rcs_preproc,'/spectra_per_sess/', pt_side_id,'/', cfg.pp_RCS_TD_subset,'/'];

% variable names are pwrspectra_by_sess, fft_bins_inHz, ch_names, respectively
load([load_dir, pt_side_id, '_pwrspectra_by_sess.mat']); %#ok<*LOAD> 
load([load_dir, pt_side_id, '_fft_bins_inHz.mat']);      
load([load_dir, pt_side_id, '_ch_names.mat']);           

parsed_db_in = readtable([load_dir, 'parsed_db_oi.xlsx']);


parsed_db_in(cellfun(@isempty, parsed_db_in.activeGroup),...
             'activeGroup') = {'A'};
%%
ch_TD_tbl = table;

for i_ch = 1:4

    tmp_tbl = table;

    hpf_ch  = parsed_db_in.(sprintf('Ch%g_hpf', i_ch-1));

    hpf_ch(hpf_ch == 8.5) = 0.85;

    ch_names(:, i_ch) = cellfun(@(x,y) sprintf('%s HPF-%g', x, y), ...
                        ch_names(:, i_ch), num2cell(hpf_ch),...
                        'UniformOutput', false);
     
    
    [tmp_tbl.chanFullStr, ~, i_chanFullStr] = unique(ch_names(:, i_ch)); %#ok<*USENS> 

    tmp_tbl.N_chanFullStr = hist(i_chanFullStr, unique(i_chanFullStr))';
    
    [N_sess, i_td ] = max(tmp_tbl.N_chanFullStr);

    i_keep                      = strcmp(ch_names(:, i_ch), tmp_tbl.chanFullStr(i_td));

    a_ch_name = sprintf('Ch%g', i_ch -1);

    ch_TD_tbl{a_ch_name, 'every TD settings'} = {tmp_tbl};
    ch_TD_tbl{a_ch_name, 'N_sess'}            = N_sess;
    ch_TD_tbl(a_ch_name, 'TD setting')        = tmp_tbl.chanFullStr(i_td);
    ch_TD_tbl(a_ch_name, 'keep_boolean')      = {i_keep};


    parsed_db_in.(sprintf('%s_chanFullStr', a_ch_name)) =ch_names(:, i_ch);
end

% ensure that every session has therapy Off or at least at 0 mA
amp_oi = nan(height(parsed_db_in), 1);

for i_sess = 1 :height(parsed_db_in)
    amp_oi(i_sess) =  parsed_db_in{i_sess, ...
                                   sprintf('Group%sProg0_ampInMilliamps', parsed_db_in.activeGroup{i_sess})...
                                    };
end

i_off_0mA      = ~parsed_db_in.therapyStatus | amp_oi ==0;
% all Off or 0 mA with the same TD sense settings
i_all_sess_oi  = i_off_0mA & all([ch_TD_tbl.keep_boolean{:}], 2);

% return parsed db, ch_names, and pwrspectra AFTER above filtering
parsed_db_out          = parsed_db_in(i_all_sess_oi, :);
ch_names_out               = ch_names(i_all_sess_oi,:);

pwrspectra_by_sess_out =  pwrspectra_by_sess(:, i_all_sess_oi, :);

%% a "per day" spectrum day is found by the mean of all the sessions that 
% took place on a given day (between 0-3 sessions generally speaking)

save_dir       = [dirs.rcs_preproc,'/spectra_per_sess/', pt_side_id,'/', cfg.pp_RCS_TD_subset,' (pp work up)/'];

if ~isfolder(save_dir);     mkdir(save_dir);     end

cfg.plotted_time     = 'since_inital_session';
cfg.z_score          = true;
cfg.freq             = [0.5, 80];
cfg.y_lbl_txt = ['Days since RCS implant', newline '(implant as day 0)'];


    plt_session_spectra(cfg, fft_bins_inHz, pt_side_id, pwrspectra_by_sess_out, ...
                                     ch_names_out, parsed_db_out , save_dir)


cfg.z_score    = false;
cfg.pwr_limits = [-10, -1];

    plt_session_spectra(cfg, fft_bins_inHz, pt_side_id, pwrspectra_by_sess_out, ...
                                    ch_names_out, parsed_db_out , save_dir);



%% implement manual artifact rejection of spectra

i_ch_sess_rej = plt_rej_sess_spectra(cfg, fft_bins_inHz, pt_side_id, pwrspectra_by_sess_out, ...
                                ch_names_out, parsed_db_out , save_dir);



%% interpolate over artifactual frequnecy bandpeaks
%%% RCS recording artifact--not stimulation

cfg.freq             = [15, 95];
i_keep               = find(i_all_sess_oi);
i_keep               = i_keep(all(~i_ch_sess_rej,2));

pwrspectra_by_sess_out =  pwrspectra_by_sess(:, i_keep, :);
parsed_db_out          = parsed_db_in(i_keep, :);
ch_names_out           = ch_names(i_keep,:);
                            

[pwrspectra_by_sess_out, art_fft_peaks_tbl]...
    ...
    = flatten_sess_spectra(...
    ...
    cfg, fft_bins_inHz, pt_side_id, pwrspectra_by_sess_out, ch_names_out, parsed_db_out , save_dir);



cfg.z_score    = false;
cfg.pwr_limits = [-10, -1];
cfg.freq       = [2, 85];

    plt_session_spectra(cfg, fft_bins_inHz, pt_side_id, pwrspectra_by_sess_out, ...
                                    ch_names_out, parsed_db_out, [save_dir, 'done/']);


i_nan = any(squeeze(isnan(pwrspectra_by_sess_out(:,:,1))));


pwrspectra_by_sess_out =  pwrspectra_by_sess_out(:, ~i_nan, :);
parsed_db_out          = parsed_db_out(~i_nan, :);
ch_names_out           = ch_names_out(~i_nan,:);
%%

ch_names_out = unique(ch_names_out, 'rows', 'stable');

save([save_dir, pt_side_id, '_pwrspectra_by_sess'], ...
    'pwrspectra_by_sess_out', '-v7.3');

save([save_dir, pt_side_id, '_fft_bins_inHz'], ...
    'fft_bins_inHz', '-v7.3');

save([save_dir, pt_side_id, '_ch_names'], ...
    'ch_names_out', '-v7.3');

save([save_dir, pt_side_id, '_art_fft_peaks_tbl'], ...
    'art_fft_peaks_tbl', '-v7.3');

%%% to ensure newest version of parsed database of interest is saved
save_file = [save_dir, pt_side_id, '_parsed_db_oi.xlsx'];

if isfile(save_file);   delete(save_file);  end

writetable(parsed_db_out, save_file);
end