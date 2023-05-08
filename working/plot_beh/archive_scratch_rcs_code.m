



%% parse db to human-readable format
i                   = cellfun(@(x) length(x) == 1, db.RCS02R.duration);
db_RCS02R           = db.RCS02R(i, :);

[~, i_u] = unique(db_RCS02R.sess_name);

db_RCS02R = db_RCS02R(i_u, :);

db_RCS02R.timeStart = cellfun(@(x) x, db_RCS02R.timeStart);
db_RCS02R.timeStop  = cellfun(@(x) x, db_RCS02R.timeStop);
db_RCS02R.duration  = cellfun(@(x) x, db_RCS02R.duration);

% take sessions of useful yet managable duration
i_sess              = ge(db_RCS02R.duration , duration('00:07:00')) &...
                      le(db_RCS02R.duration , duration('01:00:00'));

db_RCS02R           = db_RCS02R(i_sess, :);

i_stim              = cellfun(@(x) ~isempty(x), db_RCS02R.stimSettingsOut);

db_RCS02R.activeGroup(i_stim)  ...
    = cellfun(@(x) x.activeGroup, db_RCS02R.stimSettingsOut(i_stim));


db_RCS02R  = sortrows(db_RCS02R, 'timeStart', 'descend');

% try to simulate certain sessions
sess_oi    = [2, 15, 18:19];
db_RCS02R  = db_RCS02R(sess_oi,:);

nickname = {'i'; 'ii';'iii';'iv'};
db_RCS02R.sess_name = cellfun(@(x,y) [x,'_',y], db_RCS02R.sess_name, nickname, 'UniformOutput', false);

db_RCS02R.per_TD_lost = nan(height(db_RCS02R),1);

for i_sess = 1 : height(db_RCS02R)

    data_dir = db_RCS02R.path{i_sess};
    
    % data using Analysis-rcs-data
    
    [unifiedDerivedTimes,...
        timeDomainData, ~, ~,...
        ~, ~, ~, ...
        PowerData, ~, ~,...
        FFTData, ~, ~,...
        ...
        AdaptiveData, ~, ~, timeDomainSettings, powerSettings,...
        ~, eventLogTable, metaData, stimSettingsOut, stimMetaData,...
        stimLogSettings, DetectorSettings, AdaptiveStimSettings, ...
        AdaptiveEmbeddedRuns_StimSettings, ~] ...
        ...
        = ProcessRCS(data_dir, 3);
    
    dataStreams         = {timeDomainData, PowerData, AdaptiveData, FFTData};

    comb_dt = createCombinedTable(dataStreams, unifiedDerivedTimes, metaData);
    
    [comb_dt_chunks, per_TD_lost]    = chunks_and_gaps(comb_dt);

    db_RCS02R.per_TD_lost(i_sess)    = per_TD_lost;
    db_RCS02R.comb_dt_chunks(i_sess) = comb_dt_chunks;
end
%% simulate LD activity 
set(0,'DefaultFigureVisible','off')
for i_sess = 1 : height(db_RCS02R)

    sim_tbl     = td_to_fft_pb(i_sess, db_RCS02R);
end
set(0,'DefaultFigureVisible','on')
%%
% pull out all sessions w/ FFT channel streamed

fftSettings  = vertcat(db_RCS02R.fftSettings{:});
fftSettings  = struct2table(fftSettings.fftConfig);

% length of 1 implies that it is a channel name rather than 'Disabled'
i_fft_stream = find(cellfun(@(x)  length(x.fftConfig.fftStreamChannel) == 1, ...
                        db_RCS02R.fftSettings));

i_sess       = i_fft_stream(1);

td_to_fft_pb(i_sess, db_RCS02R);

%%

% need aDBS session to see actual LD outputs as validation
% --> use most recent one

i_sess                 = find(strcmp(db_RCS02R.activeGroup, 'D'));
i_sess                 = i_sess(2);



fprintf(['aDBS %s | starting at %s | duration %s (HH:MM:SS.sss)', newline], ...
    db_RCS02R.sess_name{i_sess}, ...
    db_RCS02R.timeStart(i_sess),...
    db_RCS02R.duration(i_sess));







%% unilateral implant AND used INS logs to capture PTM intiated group changes
%{
[REDcap_INSLog.RCS02R, proc_group_changes.RCS02R]  ...
    ...
    = explicit_INS_stim_params(...
    ...
cfg, INS_logs.RCS02R.group_changes, db.RCS02R, REDcap.RCS02, visits.RCS02);



%% add stim params to REDcap based off of EventLog.txt, and DeviceSettings.json files
cfg                     = [];
cfg.stage_dates         = stage_dates{2};
cfg.pt_id               = 'RCS02';

INSLog_w_redcap.RCS02R  = get_INSLog_stim_params(cfg,...
                                     INS_logs.RCS02R.group_changes,...
                                     db.RCS02R,...
                                     REDcap.RCS02,...
                                     visits.RCS02...
                                    );

[wrt_stim_REDcap.RCS02, stimGroups.RCS02] ...
    ...
    = make_stim_groups(...
    ...
'RCS02', [], INSLog_w_redcap.RCS02R, visits.RCS02);


cfg                          = [];
cfg.pt_id                    ='RCS02';
cfg.min_n_reports            = 5;

    plted_stim_groups.RCS02  = plot_stim_groups(cfg, stimGroups.RCS02);

%}


%% organize Streaming Notes, clinic dates, etc

%{
RCS06
    08/17/22 clinic visit

    instructed on Streaming Notes to use HIS timezone 
    
    sliders on pt tablet "jumps"--needs to zoom in to get VAS right where he
    wants it

    CTM occasionally drops during sessions--probably the right side--not
    100%

    for MPQ "sickening" field has been answered--last 3-7 days Crohn's Flare
%}



% experiment w/ stim sweeps parameters
%{
n_runs          = 5;            t_off           = 20;
n_iter          = 2;            t_on            = 60;


t_in_sec        = n_iter * n_runs * (t_off+ t_on);
duty_cycle      = t_on / (t_off + t_on);

t_in_min        = t_in_sec / 60

% logrhythmically spaced frequencies
freq_start       = 10;
freq_stop        = 175;

freq_exp_start   = log10(freq_start);
freq_exp_stop    = log10(freq_stop);

log_freqs        = floor(logspace(freq_exp_start, freq_exp_stop, n_runs))

lin_freqs        = floor(linspace(freq_start, freq_stop, n_runs))

% amp              = repmat(2, 1, n_runs)
pw               = repmat(300, 1, n_runs)

% element-wise multiplication based on percent of time stim-sweep is on
% TEED             = amp.^2 .* lin_freqs .* pw * (duty_cycle)

% need impendance to find mean total electrical energy delievered (TEED) per sec

TEED_now         = 2.^2 * 100 * 300 * 60/(60+20)

amp              = sqrt(TEED_now ./ (lin_freqs.* pw .* duty_cycle))


%}



% n_runs          = 5;            t_off           = 20;
% n_iter          = 2;            t_on            = 60;
% 
% 
% t_in_sec        = n_iter * n_runs * (t_off+ t_on);
% duty_cycle      = t_on / (t_off + t_on);
% 
% t_in_min        = t_in_sec / 60;
% 


%% old code:
% (previously used DeviceSettings.json alone to get stim params)
% --> now likely redundant and less accurate than
% 'align_REDcap_to_stimLog()' fxn

    % stim parameters grouping based off of 'DeviceSettings.json' from 'StimSettingsOut' var in db.RCSXXX
    %{
    %{
    * creates db_beh.RCSXX with timestamps, contacts, amp, PW, freq, cycling,
      group, and laterality (unilateral vs bilateral stim)
    
        -manually inspect laterality to handle edge cases
    
    * align nearest REDcap report as row in db_beh.RCSXX
    
    
    for RCS05 and RCS02 especially look at pain reports within session 5 - 15 minutes
       
    
    %}
    
    cfg                 = [];
    cfg.pt_id           = 'RCS04';
    cfg.pia_raw_dir     = pia_raw_dir;
    cfg.plot_sess_dur   = false;
    
       db_beh.RCS04 = db_sort_beh(cfg, db.RCS04L, db.RCS04R, REDcap.RCS04);
    
    
    cfg.pt_id = 'RCS05';
        db_beh.RCS05 = db_sort_beh(cfg, db.RCS05L, db.RCS05R, REDcap.RCS05, visits.RCS05);
    
    
    cfg.pt_id = 'RCS02';
        db_beh.RCS02 = db_sort_beh(cfg,[], db.RCS02R, REDcap.RCS02);
    
    
    cfg.pt_id = 'RCS06';
        [db_beh.RCS06] = db_sort_beh(cfg,db.RCS06L, db.RCS06R, REDcap.RCS06, visits.RCS06  );
    
    %}


    % validate StimLog.json's w/ db.RCSXXX 'stimSettingsOut' (from DeviceSettings.json)

    %  used for sanity check btwn StimLog.json and DeviceSettings.json)
    %{
    i_stimLog_redcaps = cell2mat(stimLog.i_redcap);
    
    i_emp = cellfun(@isempty, db_beh.RCS05.stimContacts);
    
    i_L = cellfun(@(x) strcmp('L', x(1)), db_beh.RCS05.stimContacts(~i_emp));
    
    db_beh_RCS05L = db_beh.RCS05(~i_emp,:);
    
    db_beh_RCS05L = db_beh_RCS05L(i_L,:);
    
    
    for i = 1 : height(db_beh_RCS05L)
    
       % per db_beh, return index in stimLog if ANY of the REDcap indices are shared
       in_stimLog_at = find(cellfun(@(x) any(ismember(db_beh_RCS05L.i_redcap_near{i}, x)), stimLog.i_redcap));
    
        if  any(~strcmp(db_beh_RCS05L.stimGroup(i), stimLog.activeGroup(in_stimLog_at)))||...
            any(db_beh_RCS05L.stimfreq(i)  ~=  stimLog.stimfreq(in_stimLog_at)) ||...
            any(db_beh_RCS05L.stimPW(i)  ~=  stimLog.stimPW(in_stimLog_at))
        
            disp([num2str(i), ' index of db_beh_RCS05L does not match ', num2str(in_stimLog_at'),...
                ' index of StimLog.json (for "iredcap_near")'])
    
        end
    
         % per db_beh, return index in stimLog if ANY of the REDcap indices are shared
       in_stimLog_at = find(cellfun(@(x) any(ismember(db_beh_RCS05L.i_redcap{i}, x)), stimLog.i_redcap));
    
        if  any(~strcmp(db_beh_RCS05L.stimGroup(i), stimLog.activeGroup(in_stimLog_at)))||...
            any(db_beh_RCS05L.stimfreq(i)  ~=  stimLog.stimfreq(in_stimLog_at)) ||...
            any(db_beh_RCS05L.stimPW(i)  ~=  stimLog.stimPW(in_stimLog_at))
        
            disp([num2str(i), ' index of db_beh_RCS05L does not match ', num2str(in_stimLog_at'),...
                ' index of StimLog.json (for "iredcap")'])
    
        end
    
    end
    
    
    all_db_beh_i_redcap = vertcat(db_beh_RCS05L.i_redcap{:});
    
    if length(all_db_beh_i_redcap) ~= length(unique(all_db_beh_i_redcap)) 
        
        disp('****The same REDcap report(s) "assigned" to different stim settings*****')
    
    end
    
    per_assigned = length(unique(all_db_beh_i_redcap)) ./ height(redcap) * 100;
    
    if  per_assigned ~= 1
    
        disp(['Only ', num2str(per_assigned),'% of REDcap report(s) assigned to stim settings'])
    
    end
    %}


    %% animate pain space visualization

%v = RecordRotation();

%
%{
playRecording(v);
%%
function V2 = RecordRotation()

n = 300;
V = zeros(n,2);
T = timer('period',0.001,'executionmode','fixedrate',...
    'TimerFcn',@captureAzEl,'TasksToExecute',n);
rotate3d on
drawnow;
start(T);
drawnow;
wait(T);
V2 = V;

      function captureAzEl(src,~)
          cnt = get(src,'TasksExecuted');
          [V(cnt,1), V(cnt,2)] = view;  
          drawnow;
      end

end



function playRecording(V)
frame_period = 0.01;
n = length(V);
T = timer('period',frame_period,'executionmode','fixedrate',...
    'TimerFcn',@playRecording,'TasksToExecute',n);
rotate3d on
drawnow;
start(T);
drawnow;
wait(T);

      function playRecording(src,~)
          cnt = get(src,'TasksExecuted');
          view(V(cnt,1), V(cnt,2));  
          drawnow;        
      end
  end

%}
