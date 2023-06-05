function format_spectra(cfg, h_sur)
    hold on 
    c_str = "log_{10}(mV^{2}/Hz)";

    ylabel(cfg.y_lbl_txt);          xlabel('Frequency (Hz)');
   
    h_sur.LineStyle = 'none';

    colormap(preset_colormap);      
    cb_hdl = colorbar;              cb_hdl.Limits          = cfg.pwr_limits; 
    
    set(gca, 'FontSize', 14);       cb_hdl.Label.String    = c_str;
    cb_hdl.FontSize     = 10;       cb_hdl.Label.FontSize  = 18;

    ax = gca;               ax.YDir = 'Reverse';
end