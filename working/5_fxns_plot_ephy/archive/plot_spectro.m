function spectro = plot_spectro(cfg)

if cfg.log10
    cfg.fftoutput = log10(cfg.fftoutput);
end

spectro = pcolor(cfg.fftTime, cfg.fftBins, cfg.fftoutput); shading interp;

caxis(cfg.cb_lim);


cb_info = colorbar;       hold on;
ax = gca;

ax.MinorGridLineStyle = 'none';
ax.GridLineStyle      = 'none';

title(cfg.title, 'FontWeight', 'normal');

xlabel('Time');     ylabel('Frequencies (Hz)'); ylim(cfg.y_lim); xlim(cfg.x_lim);
set(gca, 'FontName', 'Arial','Fontsize',14);


% place, scale, rotate, and move label for color bar
if cfg.log10
    cb_info.Label.String   = 'Log_{10}(mV/Hz)';
else
    cb_info.Label.String   = 'mV / Hz';
end
cb_info.Label.Rotation = 270;
cb_info.Label.Position = [cb_info.Label.Position(1)*1.45, cb_info.Label.Position(2)];
cb_info.Label.FontSize = 14;

end