#### Carefully import .mat of spectra and meta data
# inspect mean spectra as sanity check
def import_FT_RCS_spectra(pt_hemi, dirs):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import mat73

    base_dir= dirs['raw'] +  pt_hemi + dirs['raw_sub'] +  pt_hemi + '_'
    print(base_dir)

    # To load a .mat 7.3 into Python as a dictionary:
    ft_struct = mat73.loadmat(base_dir + 'ft_form_pp_FFT_struct.mat')['ft_freq_pp']
    
    pwrspectra_by_sess = np.reshape(ft_struct['powspctrm'], [ft_struct['powspctrm'].shape[1],
                                    ft_struct['powspctrm'].shape[0], ft_struct['powspctrm'].shape[2]], order='F')
    
    print(np.shape(pwrspectra_by_sess))

    ch_names  = ft_struct['label']
    ch_names  = [i_ch[0] for i_ch in ch_names]
    chs_split = [a_ch.split(' ')[0:2] for a_ch in ch_names]
    ch_names  = [' '.join(a_ch) for a_ch in chs_split]
    print(ch_names)

    fft_bins_inHz = ft_struct['freq']
    print(np.shape(fft_bins_inHz))

    
    parsed_db_oi   = pd.read_excel(base_dir + 'par_db.xlsx')

    #### sanity check: importing from MATLAB has all the same spectra and meta data

    plt.style.use('seaborn-colorblind')
    fig = plt.figure()
    fig.patch.set_facecolor('white')    # nescessary so background of figure is white and NOT transparent
    plt.title(pt_hemi)
    for i_chan in range(len(ch_names)): # loop through channels

        plt.plot(fft_bins_inHz, np.log10(np.mean(pwrspectra_by_sess[i_chan, :, :], 0)))


    plt.ylim([-11,-1])

    plt.xlabel('frequency (Hz)')
    plt.ylabel('log($uV^{2}$/Hz)')

    plt.grid(visible = True, color = 'k', alpha = .5, linestyle = '-')
       # plt.savefig(psd_results + ch_names[i_chan], bbox_inches='tight', dpi = 100)
    plt.legend(ch_names)
    plt.show()
    
    return pwrspectra_by_sess, fft_bins_inHz, ch_names, parsed_db_oi