function plt_daily_spec(cfg, dirs, pt_side_id)


source_dir     = fullfile(dirs.rcs_preproc, 'spectra_per_sess',...
                                    pt_side_id, ...
                                    [cfg.pp_RCS_TD_subset, ' (pp work up)'],...
                                     cfg.pp_fft);

tmp = load(...
            fullfile(source_dir, [pt_side_id, '_ft_form_pp_FFT_struct.mat'])...
            );

ft_freq_pp = tmp.ft_freq_pp;        clear tmp;


%%
i_freq_oi     = find((ft_freq_pp.freq >= cfg.freq(1) & ft_freq_pp.freq <= cfg.freq(2)));
i_ch         = str2double(cfg.rcs_ch_oi(3))+1;
raw_spectrum = log10(...
                    squeeze(ft_freq_pp.powspctrm(:,i_ch, i_freq_oi))...
                    );

i_keep       = find(~any(isnan(raw_spectrum),1));

raw_spectrum = raw_spectrum(:, i_keep);

freq_oi      = ft_freq_pp.freq(i_freq_oi);
freq_oi      = freq_oi(i_keep);

%%%
t_vec_raw         = ft_freq_pp.rcs.par_db.timeStart;

t_start = dateshift(t_vec_raw(1), 'start', 'day');
t_end   = dateshift(t_vec_raw(end), 'end', 'day');

% t_vec   = t_start:duration('24:00:00'):t_end;
% 
% plt_spectrum = nan(length(t_vec), length(freq_oi));
% 
% for i_day = 1 :length(t_vec)-1
% 
%     i_ss = find(...
%                 isbetween(ft_freq_pp.rcs.par_db.timeStart,...
%                           t_vec(i_day), t_vec(i_day+1))...
%                 );
% 
%     plt_spectrum(i_day, :)...
%         =...
%         mean(raw_spectrum(i_ss, :), 1);
% 
% end

plt_t_vec = days(t_vec_raw - t_start);

%%
%%%
figure('Units','inches','Position',[0, 3, 6, 5]);

title_txt      = sprintf('%s | %g sessions over %g days', ...
                            pt_side_id, length(plt_t_vec), days(t_end - t_start));

title(title_txt, 'FontSize', 12)

preset_colormap(@brewermap, "YlOrRd");
h_sur     = surface(freq_oi, plt_t_vec, raw_spectrum);  

format_spectra(cfg, h_sur);

%%% on 2nd y axis plot power spectra on top of each other (rather than in spectrogram)
yyaxis right;         set(gca, 'YColor', 'k'); yticklabels([]);     hold on


[lineout, ~] = stdshade(raw_spectrum, freq_oi, 0.2, 'k'); hold on
            
lineout.LineStyle = '-';              lineout.LineWidth = 3;

%plot(freq_oi, plt_spectrum,  '-k', 'LineWidth', .5, 'Marker','none');

%ylim(cfg.pwr_limits);

exportgraphics(gcf, fullfile(cfg.save_dir, [cfg.fig_name,'.png']))
end