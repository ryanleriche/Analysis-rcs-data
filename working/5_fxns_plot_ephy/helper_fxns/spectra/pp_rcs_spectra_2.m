function pp_rcs_spectra_2(cfg, dirs, pt_side_id)
set(0,'DefaultFigureVisible','off')

%%% load FieldTrip structure containing spectra and RC+S specific meta-data
load_dir = fullfile(dirs.rcs_preproc,'spectra_per_sess', pt_side_id, cfg.pp_RCS_TD_subset, cfg.pp_fft);

load(fullfile(load_dir, [pt_side_id, '_ft_form_FFT_struct.mat'])); %#ok<*LOAD>    

ft_freq_pp = ft_freq_mtm_by_win;

par_db_in = ft_freq_pp.rcs.par_db;
par_db_in(cellfun(@isempty, par_db_in.activeGroup),'activeGroup') = {'A'};
%%
% ensure that every session has therapy Off or at least at 0 mA
amp_oi = nan(height(par_db_in), 1);

for i_sess = 1 :height(par_db_in)
    amp_oi(i_sess) =  par_db_in{i_sess, ...
                                   sprintf('Group%sProg0_ampInMilliamps', par_db_in.activeGroup{i_sess})...
                                    };
end
% all Off or 0 mA with the same TD sense settings
i_off_0mA      = ~par_db_in.therapyStatus | amp_oi ==0;

% return parsed db, and power spectrum
ft_freq_pp.rcs.par_db      = par_db_in(i_off_0mA, :);
ft_freq_pp.powspctrm       = ft_freq_mtm_by_win.powspctrm(i_off_0mA, : ,:);


%% implement manual spectra artifact rejection
save_dir     = fullfile(cfg.anal_dir, 'spectra_per_sess', pt_side_id, [cfg.pp_RCS_TD_subset, ' (pp work up)'], cfg.pp_fft);

if ~isfolder(save_dir);     mkdir(save_dir);    end

cfg.freq            = [1, 95];
cfg.pwr_limits      = [-10, -1];

i_ch_sess_rej = visually_reject_spectra(cfg, save_dir, pt_side_id, ft_freq_pp);

%% interpolate over artifactual frequnecy bandpeaks
%%% RCS recording artifact--not stimulation

cfg.freq               = [15, 95];

i_keep                 = find(i_off_0mA);
i_keep                 = i_keep(all(~i_ch_sess_rej,2));


ft_freq_pp.powspctrm   = ft_freq_pp.powspctrm(i_keep, : ,:);
ft_freq_pp.rcs.par_db  = ft_freq_pp.rcs.par_db(i_keep, :);
                            

ft_freq_pp_2 = inter_over_peaks(cfg, save_dir, pt_side_id, ft_freq_pp);

%%
cfg.z_score         = false;
cfg.freq            = [0.5, 80];
cfg.pwr_limits      = [-10, -1];
cfg.title           = '3_freq_interpolation';

cfg.txt_leg         = ["non-interpolated", "interpolated"];

    plt_rcs_PSDs(cfg, pt_side_id, save_dir, ft_freq_pp, ft_freq_pp_2)

%%
% rename preprocessed FieldTrip structure for simplicity going forward
ft_freq_pp = ft_freq_pp_2;
source_dir = fullfile(dirs.rcs_preproc, 'spectra_per_sess',...
                      pt_side_id, [cfg.pp_RCS_TD_subset, ' (pp work up)'], cfg.pp_fft);

if ~isfolder(source_dir);       mkdir(source_dir);     end

save(fullfile(source_dir, [pt_side_id, '_ft_form_pp_FFT_struct']), 'ft_freq_pp', '-v7.3');


par_db_file = fullfile(source_dir, [pt_side_id, '_par_db.xlsx']);
% RBL noticed odd merging between .xlsx versions if previous was NOT
% deleted first
if isfile(par_db_file);     delete(par_db_file);    end

writetable(ft_freq_pp.rcs.par_db, par_db_file);


end