function [nk, rcs] = plt_nk_rcs_spectra(dirs, pt_sides, nk, rcs, fig_name)

    colors = brewermap(NaN, 'Dark2');
    
    
    save_dir = fullfile(dirs.rcs_pia, 'ephy_analysis/staged_spectra_stability/');
    
    if ~isfolder(save_dir);     mkdir(save_dir);    end
    
    figure('Units','inches','Position',[-15, 3, 10, 10])
    
    tiledlayout(3,3, 'Padding','tight')
    sgtitle(sprintf('Periodic power-spectra between inpatient and ambulatory stages'))
    
    
    for i = 1 : length(pt_sides)
    
        rcs_chs = rcs.sense_chan.Properties.RowNames;
        i_rcs   = find(contains(rcs_chs, pt_sides{i}));
    
        nk_chs = nk.sense_chan.Properties.RowNames;
        i_nk   = find(contains(nk_chs, pt_sides{i}));
    
        nexttile
    
        for j = 1 : length(i_rcs)
            plt_spectrum = rcs.sense_chan{i_rcs(j),"pwr_spectra_aperiodic_rmv"}{1};

            plt_spectrum(all(plt_spectrum == 0,2), :) = [];
            
%             avg_spec     = mean(plt_spectrum, 'all');
%             
%             plt_spectrum = (plt_spectrum - avg_spec);

            stdshade(plt_spectrum, rcs.fft_bins_inHz, 0.2, colors(j,:)); hold on
    
            ylim([-.25, 1.25]);
    
            xlabel('Frequency (Hz)');       ylabel('Power (db) above 1/f')
    
            rcs.sense_chan{i_rcs(j),"mean_periodic_spec"} = {mean(plt_spectrum)};

        end
    
        for j = 1 :length(i_nk)
            plt_spectrum = nk.sense_chan{i_nk(j),"pwr_spectra_aperiodic_rmv"}{1};

            plt_spectrum(all(plt_spectrum == 0,2), :) = [];
           
%             avg_spec     = mean(plt_spectrum, 'all');
%             
%             plt_spectrum = (plt_spectrum - avg_spec);

            [lineout, ~] = stdshade(plt_spectrum, nk.fft_bins_inHz, 0.2, colors(j,:)); hold on
    
            lineout.LineStyle = '-.'; lineout.LineWidth = 2.5;

            nk.sense_chan{i_nk(j),"mean_periodic_spec"} = {mean(plt_spectrum)};

        end
    
        title(pt_sides{i});
    
        rcs_only = cellfun(@(x) [x(8:end),  ' (ambulatory)'], rcs_chs(i_rcs), 'UniformOutput',false);
    
        nk_only  = cellfun(@(x) [x(8:end),  ' (inpatient)'], nk_chs(i_nk), 'UniformOutput',false);
    
      
        legend([rcs_only; nk_only]);
        
    end
    
    exportgraphics(gcf, [save_dir, fig_name, '.png']);
end
