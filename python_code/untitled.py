def plt_save_aperiodic_overtime(fg_group_results, ch_name, exp_data_dir, aper_long_dir, parsed_db_oi, pt_side, aperiodic_mode):

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scipy

    tmp_params = pd.DataFrame(fg_group_results)

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
    
    fooof_params.to_excel(exp_data_dir + ch_name + '_fooof_params.xlsx', index=False)

    ###
    t_exp_corr   = stats.spearmanr(parsed_db_oi['timeStart'], fooof_params['aperiodic_exponent'])
    t_off_corr   = stats.spearmanr(parsed_db_oi['timeStart'], fooof_params['aperiodic_offset'])
    exp_off_corr = stats.spearmanr(fooof_params['aperiodic_exponent'], fooof_params['aperiodic_offset'])

    plt_corrs = [t_exp_corr[0],t_exp_corr[1], t_off_corr[0], t_off_corr[1], exp_off_corr[0], exp_off_corr[1]]


    set_up_str = 'Time, aper_exp         | r = {:.3f} p= {:.3g}\n' + \
                 'Time, aper_offset      | r = {:.3f} p= {:.3g}\n' + \
                 'aper_exp, aper_ offset | r = {:.3f} p= {:.3g}'

    simple_stat_str   =set_up_str.format(*plt_corrs)

    df_vars       = (pt_side, ch_name, parsed_db_oi.shape[0],
                      str(parsed_db_oi['timeStart'].iloc[0])[0:10],   str(parsed_db_oi['timeStart'].iloc[-1])[0:10], 
                     )

    descr_title   = '{} {} | Aperiodic offset and exponent from {} sessions from {} -> {}'.format(*df_vars)

    fig, ax1 = plt.subplots(figsize=(15, 6))
    plt.ylim([-11, 1])
    plt.scatter(parsed_db_oi['timeStart'], 
                fooof_params['aperiodic_offset'],  
                color='blue')

    ax1.set_ylabel("Aperiodic Offset", color='blue', fontsize=14)

    # Instantiate a second axes that shares the same x-axis
    ax2 = ax1.twinx() 
    plt.scatter(parsed_db_oi['timeStart'],
                fooof_params['aperiodic_exponent'],
                color='green',
                )

    plt.title(descr_title, fontsize=18)
    plt.gcf().text(0.5, 0.2,simple_stat_str, fontsize=14)

    ax2.set_ylabel("Aperiodic Exponent", color='green', fontsize=14)
   
    if aperiodic_mode == 'fixed':
        plt.ylim([1, 3])

    plt.savefig(aper_long_dir + ch_name + '_aperiodic_longitudinally.pdf', bbox_inches='tight', dpi = 100)
    plt.close()