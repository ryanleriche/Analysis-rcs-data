def  import_NK_spectra(pt_id, dirs):
    
    import numpy as np
    import pandas as pd
    import os
    import re
    import matplotlib.pyplot as plt
    
    base_dir = dirs['raw'] + pt_id + '/'
    
    print(pt_id)
    tmp_labs = [re.split(r'(\d+)', ch_lab[4:]) for ch_lab in np.load(base_dir + 'ch_labs.npy')]
    ch_names = [ch_lab[0].replace(" ", "") + ' ' + ch_lab[1] +'-' + ch_lab[3] for ch_lab in tmp_labs]
    print(ch_names)

    tmp = np.load(base_dir + 'specs.npy')
    dim = (tmp.shape[1], tmp.shape[0], tmp.shape[2])  
    
    raw_pwrspectra = np.reshape(tmp, dim, order='F')
    print(np.shape(raw_pwrspectra)) # spectra look to dimensionally trial by ch by freq 

    fft_bins_inHz = np.load(base_dir + 'freqs.npy')
    print(np.shape(fft_bins_inHz))
    
    #### sanity check: importing from Jeremy's NK TD -> FFT processing
    plt.style.use('seaborn-colorblind')
    fig = plt.figure()
    fig.patch.set_facecolor('white')    # nescessary so background of figure is white and NOT transparent
    plt.title(pt_id)
    for i_chan in range(len(ch_names)): # loop through channels

        plt.plot(fft_bins_inHz, np.log10(np.mean(raw_pwrspectra[i_chan,:, :], 0)))


       # plt.ylim([-11,-1])

        plt.xlabel('frequency (Hz)')
        plt.ylabel('log($uV^{2}$/Hz)')

        plt.grid(visible = True, color = 'k', alpha = .5, linestyle = '-')
       # plt.savefig(psd_results + ch_names[i_chan], bbox_inches='tight', dpi = 100)
    plt.legend(ch_names)
    plt.show()
    
    arm1_df =  pd.read_csv(base_dir + 'arm1_df.csv')
    
    return raw_pwrspectra, fft_bins_inHz, ch_names, arm1_df