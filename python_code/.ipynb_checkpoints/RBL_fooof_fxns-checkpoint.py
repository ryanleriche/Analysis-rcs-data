def setup_fooof_dirs(pt_id, dirs,fooof_param):
    
    import os
    import shutil
    # make nested folders for organized outputs
    ### *PREVIOUS RUNS ARE DELETED* FOR TROUBLESHOOTING PURPOSES AND CLARITY
    base_dir   = dirs['proc']+'/' + pt_id +'/' + dirs['proc_sub'] + ' (' + fooof_param.aperiodic_mode +')' +'/'

    fooof_dirs= {
                'error':        base_dir + 'error_ge_' + str(fooof_param.err_thresh) + '/',
                'aper_fit':     base_dir + 'aperiodic_fits/',
                'pwr_spectra':  base_dir + 'pwr_spectra (aperiodic_fit_removed)/',
                'aper_long':    base_dir + 'aperiodic_longitudinal/',
                'fooof_spec':   base_dir + 'fooofed_spectra (aperiodic_fit_removed)/',
                'fooof_data':   dirs['raw'] +'/'+ pt_id + dirs['raw_sub'] +'/' + 'exported_fooof_data/',
                'base_dir':     base_dir
                }


    for a_dir in fooof_dirs:
        
        tmp_dir = fooof_dirs[a_dir]
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
            
    return  dirs.update(fooof_dirs)



    
### plot/save the aperiodic offset and exponent and find correlations with time
def plt_save_fooof_analysis(pt_hemi, dirs, ch_name, 
                                fooof_param, 
                                t_peri_ephy, fg):

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scipy
    from scipy import stats
    import os
    import shutil
    
    from fooof import FOOOFGroup
    from fooof import FOOOF
    from fooof.plts.spectra import plot_spectrum

    
    ###
    df_vars       = (pt_hemi, ch_name, len(t_peri_ephy), 
                     str(t_peri_ephy[0])[0:10],  str(t_peri_ephy[len(t_peri_ephy)-1]))

    descr_title   = '{} | {} | spectra from {} sessions from {} -> {}'.format(*df_vars)


    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(15, 7))
    
    # using padding
    fig.tight_layout(pad=5.0)

    for i_sess in range(len(fg)):

        fm = fg.get_fooof(ind=i_sess, regenerate=True)
        if i_sess == 0:
            aperiodic_comp               = np.zeros([len(t_peri_ephy), len(fm.freqs)])
            fooof_spectra_aperiodic_rmv  = np.zeros([len(t_peri_ephy), len(fm.freqs)])
            raw_spectra                  = np.zeros([len(t_peri_ephy), len(fm.freqs)])
 
        tmp_spec = np.subtract(fm.fooofed_spectrum_, fm._ap_fit)

        axes[0, 0].plot(fm.freqs, fm.power_spectrum,       False, alpha=0.5)
        axes[0, 1].plot(fm.freqs, fm.fooofed_spectrum_,    False, alpha=0.5)
        axes[0, 2].plot(fm.freqs, fm._ap_fit,              False, alpha=0.5)
        axes[0, 3].plot(fm.freqs, tmp_spec,                False, alpha=0.5)
        
        raw_spectra[i_sess,:]                 += fm.power_spectrum
        aperiodic_comp[i_sess,:]              += fm._ap_fit
        fooof_spectra_aperiodic_rmv[i_sess,:] += tmp_spec
        

    axes[0, 0].set_xlabel('frequency (Hz)')
    axes[0, 1].set_xlabel('frequency (Hz)')
    axes[0, 2].set_xlabel('frequency (Hz)')
    axes[0, 3].set_xlabel('frequency (Hz)')
    
    axes[0, 0].set_ylabel('Raw power (dB)')
    axes[0, 1].set_ylabel('FOOOFed power (dB)')
    axes[0, 2].set_ylabel('Aperiodic power (dB)')
    axes[0, 3].set_ylabel('Periodic power (dB)')
    
    
    ## save raw spectra, flattened fooof spectra, and aperiodic component to xlsx files
    
    pd.DataFrame(raw_spectra).to_excel(dirs['fooof_data'] + ch_name + '_raw_spectra.xlsx', index=False)
    
    pd.DataFrame(fooof_spectra_aperiodic_rmv).to_excel(
        dirs['fooof_data'] + ch_name + '_pwr_spectra_aperiodic_rmv.xlsx', index=False)

    pd.DataFrame(aperiodic_comp).to_excel(dirs['fooof_data'] + ch_name + '_aperiodic_fit.xlsx', index=False)
       
    ###
    plt.suptitle(descr_title)
   

    tmp_params = pd.DataFrame(fg.group_results)

    aperiodic_offset     = [x[0] for x in tmp_params['aperiodic_params']]

    aperiodic_exponent   = [x[1] for x in tmp_params['aperiodic_params']]
    center_freq          = []
    peak_power           = []
    band_width           = []

    for i_sess in range(tmp_params.shape[0]):
        peak_params = []
        peak_params = tmp_params['peak_params'].iloc[i_sess]


        center_freq.append([x[0] for x in peak_params])
        peak_power.append([x[1] for x in peak_params])
        band_width.append([x[2] for x in peak_params])

    values        = [aperiodic_offset, aperiodic_exponent, center_freq, peak_power, band_width, tmp_params['r_squared'], tmp_params['error'], tmp_params['gaussian_params']]

    lbls          = ['aperiodic_offset', 'aperiodic_exponent', 
                     'center_freq', 'peak_power', 'band_width',
                    'r_squared', 'error','gaussian_params' ]

    fooof_params = pd.DataFrame(data = values,  index = lbls).transpose()
    
    fooof_params.to_excel(dirs['fooof_data'] + ch_name + '_fooof_params.xlsx', index=False)

    ###
    days = (t_peri_ephy-t_peri_ephy[0]).dt.total_seconds() / (3600 * 24)
    
    t_exp_corr   = stats.spearmanr(days, fooof_params['aperiodic_exponent'])
    t_off_corr   = stats.spearmanr(days, fooof_params['aperiodic_offset'])
    exp_off_corr = stats.spearmanr(fooof_params['aperiodic_exponent'], fooof_params['aperiodic_offset'])

    plt_corrs = [t_exp_corr[0],t_exp_corr[1], t_off_corr[0], t_off_corr[1], exp_off_corr[0], exp_off_corr[1]]


    set_up_str = 'corr(time, exponent)--> $r_s$ = {:.2f}, p= {:.2g}\n' + \
                 'corr(time, offset)-->$r_s$ = {:.2f}, p= {:.2g}\n' + \
                 'corr(exponent, offset)-->$r_s$ = {:.2f}, p= {:.2g}'

    simple_stat_str   =set_up_str.format(*plt_corrs)


    axes[1,2].scatter(days, fooof_params['aperiodic_offset'],  marker = 'o', alpha = 0.5)
    axes[1,2].set_ylabel('AP Offset (dB) as "●"', fontsize=10)
    axes[1,2].set_xlabel("Days from start", fontsize=10)
    # Instantiate a second axes that shares the same x-axis
    ax2 = axes[1,2].twinx() 
    

    plt.scatter(days, fooof_params['aperiodic_exponent'], marker = 'x', alpha = 0.5)

    plt.gcf().text(0.55, 0.45,simple_stat_str, fontsize=6.5)

    ax2.set_ylabel('AP Exponent "X"', fontsize=10)
   
    #if fooof_param.aperiodic_mode == 'fixed':
        #plt.ylim([1, 3])
        
    axes[1, 0].set_ylabel('Raw pain decoding')
    axes[1, 1].set_ylabel('FOOOFed pain decoding')
    axes[1, 3].set_ylabel('Periodic pain decoding')
    
    axes[1, 0].set_xlabel('frequency (Hz)')
    axes[1, 1].set_xlabel('frequency (Hz)')
    axes[1, 3].set_xlabel('frequency (Hz)')
        
    
    plt.savefig(dirs['base_dir'] + ch_name + '_raw_aperiodic_periodic_spectra.pdf', bbox_inches='tight', dpi = 100)
    plt.close()
    
    ### view erroneuous FOOOF fits
        # could iterate with lower errors for parameter optimization
    to_check    = []
    i_to_check  = []
    
    if os.path.exists(dirs['error']):
        shutil.rmtree(dirs['error'])
    if not os.path.exists(dirs['error']):
        os.makedirs(dirs['error'])

    
    for ind, res in enumerate(fg):
        if res.error > fooof_param.err_thresh:
            to_check.append(fg.get_fooof(ind, regenerate=True))
            i_to_check.append(ind)

    for ind, fm in enumerate(to_check):
        _ = fm.plot()
        fm.save_report(dirs['error'] + ch_name+ '_spectrum_index_'+str(i_to_check[ind]) +'.pdf')
        plt.close()
        
    return aperiodic_comp, raw_spectra, fooof_spectra_aperiodic_rmv, fm.freqs





def plt_save_avg_peaks(fooof_spectra_aperiodic_rmv, aperiodic_comp,freq_oi, pt_id, ch_name, t_peri_ephy, fooof_param, dirs):

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scipy
    from scipy import stats
    import os
    import shutil
    
    from fooof import FOOOFGroup
    from fooof import FOOOF
    
    ###
    df_vars       = (pt_id, ch_name, len(t_peri_ephy), 
                     str(t_peri_ephy[0])[0:10],  str(t_peri_ephy[len(t_peri_ephy)-1]))

    descr_title   = '{} | {} | MEAN of FOOOFed spectra from {} sessions from {} -> {}'.format(*df_vars)

    mean_spectrum =  np.power(10, np.mean(fooof_spectra_aperiodic_rmv, 0) + np.mean(aperiodic_comp, 0))

    fm = FOOOF(
               peak_width_limits  = fooof_param.peak_width_limits,  
               max_n_peaks        = fooof_param.max_n_peaks, 
               min_peak_height    = .1,
               peak_threshold     = fooof_param.peak_threshold,
               aperiodic_mode     = fooof_param.aperiodic_mode
               )

    fm.fit(freq_oi, mean_spectrum, fooof_param.freq_range)

    plt.style.use('seaborn-colorblind')
    fig = plt.figure()
    fig.patch.set_facecolor('white') 

    # Plot an example power spectrum, with a model fit
    fm.plot(plot_peaks='shade', peak_kwargs={'color' : 'green'})

    plt.title(descr_title)
    plt.grid(visible = True, color = 'k', alpha = .5, linestyle = '-')
    # plt.savefig(psd_results + ch_names[i_chan], bbox_inches='tight', dpi = 100)

    plt.savefig(dirs['base_dir'] + ch_name + '_mean_FOOOF.pdf', bbox_inches='tight', dpi = 100)
    plt.close()
    
    ### save to .xlsx
    aperiodic_offset     = fm.aperiodic_params_[0]
    aperiodic_exponent   = fm.aperiodic_params_[1]
    center_freq          = []
    peak_power           = []
    band_width           = []

    peak_params = fm.peak_params_
    center_freq.append([x[0] for x in peak_params])
    peak_power.append([x[1] for x in peak_params])
    band_width.append([x[2] for x in peak_params])

    values        = [aperiodic_offset,
                     aperiodic_exponent, 
                     center_freq[0],
                     peak_power[0], 
                     band_width[0], 
                     fm.r_squared_, 
                     fm.error_, 
                     fm.gaussian_params_[0]]

    lbls          = ['aperiodic_offset',
                     'aperiodic_exponent', 
                     'center_freq', 
                     'peak_power', 
                     'band_width',
                     'r_squared', 
                     'error',
                     'gaussian_params' ]

    fooof_params = pd.DataFrame(data = values,  index = lbls).transpose()
    fooof_params.to_excel(dirs['fooof_data'] + ch_name + '_MEAN_fooof_params.xlsx', index=False)

    return fm