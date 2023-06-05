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
def plt_save_aperiodic_overtime(pt_hemi, dirs, ch_name, 
                                fooof_param, 
                                t_peri_ephy, fg):

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scipy
    from scipy import stats
    
    from fooof import FOOOFGroup
    from fooof import FOOOF
    from fooof.plts.spectra import plot_spectrum

    
    ###
    
    df_vars       = (pt_hemi, ch_name, len(t_peri_ephy), 
                     str(t_peri_ephy[0])[0:10],  str(t_peri_ephy[len(t_peri_ephy)-1]), 
                     fooof_param.freq_range[0], fooof_param.freq_range[1])

    descr_title   = '{} {} | Aperiodic fit spectra from {} sessions from {} -> {}\nfreq range {}–{} Hz'.format(*df_vars)


    _, ax = plt.subplots(figsize=(7, 4))

    for i_sess in range(len(fg)):

        fm = fg.get_fooof(ind=i_sess, regenerate=True)
        if i_sess == 0:
            aperiodic_comp             = np.zeros([len(t_peri_ephy), fm._ap_fit.shape[0]])

        plot_spectrum(fm.freqs, fm._ap_fit, False, alpha=0.5, linestyle='dashed', ax=ax)
        
        aperiodic_comp[i_sess, :] = fm._ap_fit

    plt.title(descr_title)
    plt.xlabel('frequency (Hz)')
    plt.ylabel('log($mV^{2}$/Hz)')
    #plt.ylim([-11, 1])
    plt.savefig(dirs['aper_fit'] + ch_name + '_aperiodic_fit.pdf', bbox_inches='tight', dpi = 100)
    plt.close()

    ## save to xlsx file
    pd.DataFrame(aperiodic_comp).to_excel(dirs['fooof_data'] + ch_name + '_aperiodic_fit.xlsx', index=False)
    
    
    
    ###

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
    t_exp_corr   = stats.spearmanr(t_peri_ephy, fooof_params['aperiodic_exponent'])
    t_off_corr   = stats.spearmanr(t_peri_ephy, fooof_params['aperiodic_offset'])
    exp_off_corr = stats.spearmanr(fooof_params['aperiodic_exponent'], fooof_params['aperiodic_offset'])

    plt_corrs = [t_exp_corr[0],t_exp_corr[1], t_off_corr[0], t_off_corr[1], exp_off_corr[0], exp_off_corr[1]]


    set_up_str = 'Time, aper_exp         | r = {:.3f} p= {:.3g}\n' + \
                 'Time, aper_offset      | r = {:.3f} p= {:.3g}\n' + \
                 'aper_exp, aper_ offset | r = {:.3f} p= {:.3g}'

    simple_stat_str   =set_up_str.format(*plt_corrs)

    df_vars       = (pt_hemi, ch_name, len(t_peri_ephy),
                      str(t_peri_ephy.iloc[0])[0:10],   str(t_peri_ephy[len(t_peri_ephy)-1])[0:10], 
                     )

    descr_title   = '{} {} | Aperiodic offset and exponent from {} sessions from {} -> {}'.format(*df_vars)

    fig, ax1 = plt.subplots(figsize=(12, 6))
   # plt.ylim([-11, 1])
    plt.scatter(t_peri_ephy, 
                fooof_params['aperiodic_offset'],  
                color='blue')

    ax1.set_ylabel("Aperiodic Offset", color='blue', fontsize=14)

    # Instantiate a second axes that shares the same x-axis
    ax2 = ax1.twinx() 
    plt.scatter(t_peri_ephy,
                fooof_params['aperiodic_exponent'],
                color='green',
                )

    plt.title(descr_title, fontsize=14)
    plt.gcf().text(0.5, 0.2,simple_stat_str, fontsize=10)

    ax2.set_ylabel("Aperiodic Exponent", color='green', fontsize=12)
   
    #if fooof_param.aperiodic_mode == 'fixed':
        #plt.ylim([1, 3])

    plt.savefig(dirs['aper_long'] + ch_name + '_aperiodic_longitudinally.pdf', bbox_inches='tight', dpi = 100)
    plt.close()
   

    
def plt_save_periodic(pt_hemi, dirs, ch_name,
                       fooof_param,
                       t_peri_ephy, spectra_oi, fg):
    
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scipy
    from scipy import stats
    
    from fooof import FOOOFGroup
    from fooof import FOOOF
    from fooof.plts.spectra import plot_spectrum
    
    # first, save group report
    fg.save_report(dirs['fooof_spec']  + ch_name + '_fooof_Group_report.pdf')
    
    ### view erroneuous FOOOF fits
        # could iterate with lower errors for parameter optimization
    to_check    = []
    i_to_check  = []
    for ind, res in enumerate(fg):
        if res.error > fooof_param.err_thresh:
            to_check.append(fg.get_fooof(ind, regenerate=True))
            i_to_check.append(ind)

    for ind, fm in enumerate(to_check):
        _ = fm.plot()
        fm.save_report(dirs['error'] + ch_name+ '_spectrum_index_'+str(i_to_check[ind]) +'.pdf')
        plt.close()


    ### save/plot fooofed spectra - aperiodic fit
    # this is the most processe version of the spectra--showing fitted peaks only
   
    _, ax = plt.subplots(figsize=(8, 6))
    for i_sess in range(len(fg)):
        fm = fg.get_fooof(ind=i_sess, regenerate=True)
        
        if i_sess== 0:
            fooof_spectra_aperiodic_rmv  = np.zeros([spectra_oi.shape[0], fm._ap_fit.shape[0]])
 
        tmp_spec = fm.fooofed_spectrum_ - fm._ap_fit

        plot_spectrum(fm.freqs, tmp_spec, False,
                      alpha=0.5, linestyle='dashed', ax=ax)

        fooof_spectra_aperiodic_rmv[i_sess, :] = tmp_spec

    avg_wo_aperiodic = np.mean(fooof_spectra_aperiodic_rmv,0)

    plot_spectrum(fm.freqs, avg_wo_aperiodic , False,
                      color='black', alpha=0.5, linestyle='solid', ax=ax)

    plt.xlabel('frequency (Hz)')
    plt.ylabel('log($mV^{2}$/Hz)')

    plt.ylim([-1.25, 1.25])
    plt.title('Per session FOOOFed power spectra with aperiodic component removed')


    plt.savefig(dirs['fooof_spec'] + ch_name + '_fooof_spectra_aperiodic_rmv.pdf', bbox_inches='tight', dpi = 100)
    plt.close()

    ## save to xlsx file
    pd.DataFrame(fooof_spectra_aperiodic_rmv).to_excel(
        dirs['fooof_data'] + ch_name + '_pwr_spectra_aperiodic_rmv.xlsx', index=False)

    ### repeat, but now for orignal power spectrum minus aperiodic fit 
    ## (no need to save as the raw and aperiodic fit as saved seperately already)

    pwr_spectra_aperiodic_rmv  = np.zeros([spectra_oi.shape[0], fm._ap_fit.shape[0]])

    _, ax = plt.subplots(figsize=(12, 10))
    for i_sess in range(len(fg)):
        fm = fg.get_fooof(ind=i_sess, regenerate=True)
        tmp_spec = fm.power_spectrum - fm._ap_fit

        plot_spectrum(fm.freqs, tmp_spec, False,
                      alpha=0.5, linestyle='dashed', ax=ax)

        pwr_spectra_aperiodic_rmv[i_sess, :] = tmp_spec

    avg_wo_aperiodic = np.mean(pwr_spectra_aperiodic_rmv,0)

    plot_spectrum(fm.freqs, avg_wo_aperiodic , False,
                      color='black', alpha=0.5, linestyle='solid', ax=ax)

    plt.xlabel('frequency (Hz)')
    plt.ylabel('log($mV^{2}$/Hz)')

    plt.ylim([-1.25, 1.25])
    plt.title('Per session raw power spectra with aperiodic component removed')

    plt.savefig(dirs['pwr_spectra'] + ch_name + '_pwr_spectra_aperiodic_rmv.pdf', bbox_inches='tight', dpi = 100)
    plt.close()
    
    
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
            aperiodic_comp               = np.zeros([len(t_peri_ephy), fm._ap_fit.shape[0]])
            fooof_spectra_aperiodic_rmv  = aperiodic_comp 
            raw_spectra                  = aperiodic_comp
 
        tmp_spec = fm.fooofed_spectrum_ - fm._ap_fit

        axes[0, 0].plot(fm.freqs, fm.power_spectrum,       False, alpha=0.5)
        axes[0, 1].plot(fm.freqs, fm.fooofed_spectrum_,    False, alpha=0.5)
        axes[0, 2].plot(fm.freqs, fm._ap_fit,              False, alpha=0.5)
        axes[0, 3].plot(fm.freqs, tmp_spec,                False, alpha=0.5)
        
        aperiodic_comp[i_sess, :]              = fm._ap_fit
        raw_spectra[i_sess, :]                 = fm.power_spectrum
        fooof_spectra_aperiodic_rmv[i_sess, :] = tmp_spec

    axes[0, 0].set_xlabel('frequency (Hz)')
    axes[0, 1].set_xlabel('frequency (Hz)')
    axes[0, 2].set_xlabel('frequency (Hz)')
    axes[0, 3].set_xlabel('frequency (Hz)')
    
    axes[0, 0].set_ylabel('Raw power (dB)')
    axes[0, 1].set_ylabel('FOOOFed power (dB)')
    axes[0, 2].set_ylabel('Aperiodic power (dB)')
    axes[0, 3].set_ylabel('Periodic power (dB)')
    
    

    ## save raw spectra, flattened fooof spectra, and aperiodic component to xlsx files
    
    pd.DataFrame(aperiodic_comp).to_excel(dirs['fooof_data'] + ch_name + '_raw_spectra.xlsx', index=False)
    
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
    
    ###
    # save group report
    fg.save_report(dirs['base_dir']  + ch_name + '_fooof_Group_report.pdf')
    
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
