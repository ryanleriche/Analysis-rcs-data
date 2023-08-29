function [t_vec, state_vec, amp_vec, rate_vec, stream_sess_vec, on_off_vec]...
        ...
        = spon_aDBS_status(...
        ...
        plt_app_oi, start_time, stop_time, step_dur)

t_vec           = start_time:step_dur:stop_time;

state_vec       = nan(length(t_vec),1);
amp_vec         = state_vec;
rate_vec        = state_vec;
stream_sess_vec = state_vec;


 for h = 2:height(plt_app_oi)

    % for t_vec between AppLog entries, flesh out state/stim 
    i_t_vec = find(isbetween(t_vec, ...
                   plt_app_oi.time_INS(h-1), plt_app_oi.time_INS(h))...
                   );

    state_vec(i_t_vec) = plt_app_oi.oldstate(h);

    % 15 is when previous state was therapyStatus Off  
    if plt_app_oi.oldstate(h) == 15
             amp_vec(i_t_vec)       = NaN;
             rate_vec(i_t_vec)      = NaN;

    % stim vector according to amplitude setting of OLD state NOT current state
    else
        amp_vec(i_t_vec)       = plt_app_oi.prog0mA(h);
        rate_vec(i_t_vec)      = plt_app_oi.rateHz(h);
    end

    stream_sess_vec(i_t_vec)    = plt_app_oi.sess_w_same_settings(h);

end

% any amplitude of stim > 0 mA consider aDBS on
on_off_vec                = 100*(amp_vec> 0);

end
