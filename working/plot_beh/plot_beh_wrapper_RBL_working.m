%% user-inputs
% where RCS files are saved from PIA server
pia_raw_dir     = '/Users/Leriche/pia_server/datastore_spirit/human/rcs_chronic_pain/rcs_device_data/raw/';

% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
github_dir      = '/Users/Leriche/Github/';

% where DropBox desktop is saved locally
dropbox_dir     = ['/Users/Leriche/Dropbox (UCSF Department of Neurological Surgery)/',...
                   'SUBNETS Dropbox/Chronic Pain - Activa and Summit 2.0'];

% application programming interface (API) token which is essentially a
% password to access REDcap remotely, and is unique per researcher per
% study (e.g., Ryan has a unique token for the RCS and PCS studies)

rcs_API_token   = '95FDE91411C10BF91FD77328169F7E1B';
pcs_API_token   = 'DB65F8CB50CFED9CA5A250EFD30F10DB';

% pulls/organizes arms from REDcap (go into fxn to add new arms)
cd([github_dir, 'Analysis-rcs-data/working']);         

addpath(genpath([github_dir, 'Analysis-rcs-data/']));
addpath(genpath([github_dir, 'rcs-simulation/']));


REDcap                = RCS_redcap_painscores(rcs_API_token);

% stage dates and, home/clinic visits for RCS pts 1-7 w/ brief descriptions
[visits, stage_dates] = make_visit_dates;


% need to further distill by pts initals
% fluct                 = RCS_redcap_painscores(rcs_API_token, pcs_API_token, {'FLUCT'});
%% import RCS databases, and INS logs per pt side
%{

* saves RCS session summaries as databases (db) in struct per pt side
* badsessions w/n 'bs' struct

(pg. 12 of the 4NR010 Research Lab Programmer Guide M979053A001)

"Multiple programs can be combined to define a group. Each group and its associated programs can be used
to provide a therapy for specific symptoms or specific patient activities. Pulse width, amplitude, amplitude
limits, and electrode polarity are programmed separately for each program within the group (ie, each
program within the group can have different values). Pulse width limits, rate, rate limits, SoftStart/Stop,
Cycling, and Active Recharge are programmed for each group (ie, each program within the group will have
the same values)."

%}

cfg                    = [];
cfg.load_EventLog      = true;
cfg.ignoreold          = false;
cfg.raw_dir            = pia_raw_dir;

%pt_sides               = {'RCS02R'};

pt_sides               = {'RCS02R','RCS05L','RCS07L'};


%pt_sides               = {'RCS02R','RCS04L','RCS04R','RCS05L','RCS05R','RCS06L','RCS06R','RCS07L','RCS07R'};


for i = 2 : length(pt_sides)

    cfg.pt_id_side                     = pt_sides{i};

    [db.(pt_sides{i}), bs.(pt_sides{i})]...
        ...
        =  makeDatabaseRCS_Ryan(cfg);


    INS_logs.(pt_sides{i})...
        ...
        = RCS_logs(cfg);

end
%%
tt = datetime(1674492396175/1000,...
                    'ConvertFrom','posixTime','TimeZone','America/Los_Angeles',...
                    'Format','dd-MMM-yyyy HH:mm:ss.SSS')


%datetime(1970,1,1,-8,0,0,1674492396175) -datetime(2000,3,1,-8,0,722516198)



%temp_time = 

temp_time.TimeZone = 'America/Los_Angeles';
%%

%%
%{





%}


%clear app_SS_tbl proc_app_log TD_FFT_PB_LD_State_tbl

%pt_sides               = {'RCS02R','RCS04L','RCS04R','RCS05L','RCS05R','RCS06L','RCS06R','RCS07L','RCS07R'};

pt_sides               = {'RCS02R','RCS05L','RCS07L'};



for i=1 %:length(pt_sides)

    [app_SS_tbl.(pt_sides{i}), proc_app_log.(pt_sides{i}), TD_FFT_PB_LD_State_tbl.(pt_sides{i})]...
    ...
        = build_sense_summary_tbl(...
    ...
    db.(pt_sides{i}), INS_logs.(pt_sides{i}).app);

end



%%

%%

cfg                    = [];
cfg.stage_dates        = stage_dates{2};
cfg.pt_id              = 'RCS02';

%%
%
%
%
%


    [REDcap_INSLog.RCS02R, proc_group_changes.RCS02R]  ...
        ...
        = explicit_INS_stim_params(...
        ...
    cfg, INS_logs.RCS02R.group_changes, db.RCS02R, REDcap.RCS02, visits.RCS02);


%%
%
%
%
%










time_API   = datetime(db_RCSXXX.eventLogTable{1809,1}.HostUnixTime /1000,...
                        'ConvertFrom','posixTime',...
                        'TimeZone','America/Los_Angeles',...
                         'Format','dd-MMM-yyyy HH:mm:ss.SSS')


%%
%{
i                   = cellfun(@(x) length(x) == 1, db.RCS04R.duration);
db_RCS04R           = db.RCS04R(i, :);

[~, i_u] = unique(db_RCS04R.sess_name);

db_RCS04R = db_RCS04R(i_u, :);

db_RCS04R.timeStart = cellfun(@(x) x, db_RCS04R.timeStart);
db_RCS04R.timeStop  = cellfun(@(x) x, db_RCS04R.timeStop);
db_RCS04R.duration  = cellfun(@(x) x, db_RCS04R.duration);

% take sessions of useful yet managable duration
i_sess              = ge(db_RCS04R.duration , duration('00:05:00'));

db_RCS04R = db_RCS04R(i_sess,:);

i                   = cellfun(@(x) length(x) == 1, db.RCS04L.duration);
db_RCS04L           = db.RCS04L(i, :);

[~, i_u] = unique(db_RCS04L.sess_name);

db_RCS04L = db_RCS04L(i_u, :);

db_RCS04L.timeStart = cellfun(@(x) x, db_RCS04L.timeStart);
db_RCS04L.timeStop  = cellfun(@(x) x, db_RCS04L.timeStop);
db_RCS04L.duration  = cellfun(@(x) x, db_RCS04L.duration);

% take sessions of useful yet managable duration
i_sess              = ge(db_RCS04L.duration , duration('00:05:00'));

db_RCS04L = db_RCS04L(i_sess,:);



%db_RCSXX = sortrows([db_RCS04L; db_RCS04R], 'timeStart');
%}
    %%
    eventLog_jsons = vertcat(db.RCS04R.eventLogTable{:});

    u_event_names  = unique(eventLog_jsons.EventType);

    i_lead_int     = strcmp(eventLog_jsons.EventType, 'Lead Integrity');

    lead_int_tbl   = eventLog_jsons(i_lead_int, :);
%%
% lead integrity test is likely in kiloohms
ohms = kiloohms / 1000;
amps = milliamps * 1000;


volts = ohms * amps;
%%



%% from databases, parse through StimLog.json files and align to REDcap

for i = 1 : length(pt_sides)
    cfg.pt_id                       = pt_sides{i}(1:end-1);
    cfg.stage_dates                 = stage_dates{str2double(pt_sides{i}(end-1))};

    if ~strcmp(pt_sides{i}, 'RCS02R')

        [stimLog.(pt_sides{i}), REDcap.(pt_sides{i})] ...
            ...
            = align_REDcap_to_stimLog(...
            ...
        cfg, db.(pt_sides{i}), REDcap.(pt_sides{i}(1:end-1)));

    end
end



%% add stim params to REDcap based off of EventLog.txt, and DeviceSettings.json files
cfg                    = [];
cfg.stage_dates        = stage_dates{2};
cfg.pt_id              = 'RCS02';

INSLog_w_redcap.RCS02R  = get_INSLog_stim_params(cfg,...
                                     INS_logs.RCS02R.group_changes,...
                                     db.RCS02R,...
                                     REDcap.RCS02,...
                                     visits.RCS02...
                                    );



%{
% 


* plot of the State of the device overlayed with the mA delivered at that time

    * all the meta data (closest Streaming sessions, date(s) in DD/MM/YYYY, 
      stim contacts, stim freq, sensing contacts and power-bands, LD parameters, etc.,)

* report duty cycle over specified timeframe 

* report min, max, mean, and std of State durations

* visualize distribution of State durations 

* report TEED over specified timeframe

    * use impedance checks calculate TEED per second

%}
%%





%% unilateral implant AND used INS logs to capture PTM intiated group changes
[wrt_stim_REDcap.RCS02, stimGroups.RCS02] ...
    ...
    = make_stim_groups(...
    ...
'RCS02', [], INSLog_w_redcap.RCS02R, visits.RCS02);


cfg                          = [];
cfg.pt_id                    ='RCS02';
cfg.min_n_reports            = 5;

    plted_stim_groups.RCS02  = plot_stim_groups(cfg, stimGroups.RCS02);

%% generate box plots of pain metrics wrt stim parameters
%{

input(s)

___________________________________________________________________________

* w/ aligned REDcap reports, first visualize the parameter space per each
   "major contact pair" (best analgesic contact pairs as assessed clinically--
    heuristically those contacts w/ a 10+ pain reports)

* box plots of NRS, VAS, MPQ, and PC1 (for RCS05 only) based of
  parsimonious amp-PW-freq space *per* contact


___________________________________________________________________________

output(s)
wrt_stim_REDcap.RCS0X | combined left and right stim params PER REDcap survey
stimGroups.RCS0X      | surveys grrouped per major contact pair, stim OFF, @ 0mA, etc


DEFINE baseline pain from pre-Stage 0 fluct data

* evaluate 50% and 30% decrease 
    (see Farrar et al., 2011 for 30% justification based off of patient global
    impression of change (PGIC) of much improved -> very much improved 
    corresponding to 30% decrease across varying pain etiologies, 
    placebo vs pregabalin, ages, btwn, sexes, etc.)


RCS04 stats framework:

see for multivariate approach (predict many pain metrics)
https://www.mathworks.com/help/stats/specify-the-response-and-design-matrices.html


* run Kolmogorov–Smirnov test (see's if data are normal--explore pain
  metric distribtion more generally)

* run Kruskal-Walis test for the pain metrics across each of the contact-pairs 
  (grouping all freq, amp, and PW parameters w/n a contact--lack statistical 
  power for freq/PW/amp)

* for LCaud, RThal, and bilateral stim contacts look run run a Kruskal-Walis 
  test across amplitude w/n contacts 

* follow-up with a Wilconxin signed rank test (?) for signifigant groups

%}

% RCS04, RCS05, RCS06, and RCS07 can be handled together
pts = {'RCS04', 'RCS05', 'RCS06', 'RCS07'};

% pts = {'RCS06', 'RCS07'};

for i = length(pts)

    [wrt_stim_REDcap.(pts{i}), stimGroups.(pts{i})] ...
    ...
    = make_stim_groups(...
    ...
    pts{i}, REDcap.([pts{i}, 'L']), REDcap.([pts{i}, 'R']), visits.(pts{i}));


    cfg                         = [];
    cfg.min_n_reports           = 5;
    cfg.pt_id                   = pts{i};


    plted_stim_groups.(pts{i})  ...
        ...
        = plot_stim_groups(...
        ...
    cfg, stimGroups.(pts{i}));

end


% RCS04 --> stim groups after starting buprenorphine (after July 2022 home visit)

i_epoch          = find(ge(REDcap.RCS04.time, visits.RCS04.dates(11) + duration('24:00:00')));

[~, stimGroups.RCS04_postJul22] ...
    ...
    = make_stim_groups(...
    ...
 'RCS04', REDcap.RCS04L(i_epoch : end,:), REDcap.RCS04R(i_epoch: end, :), visits.RCS04);


cfg                          = [];
cfg.pt_id                    ='RCS04 (Jul 13th to Now)';
cfg.min_n_reports            = 2;

    plted_stim_groups.RCS04_postJul22 = plot_stim_groups(cfg, stimGroups.RCS04_postJul22);

%% distributions of pain metrics and relationship BETWEEN pain metrics
pts = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};

cfg                     = [];
cfg.dates               = 'AllTime';

for i = 4%1:length(pts)
    cfg.pt_id  = pts{i};       

    %plot_hist(cfg, REDcap);          
    plot_versus(cfg, REDcap);
end

%% RCS02: explore NRS and VAS "mismatch" 
cfg.pt_id               = 'RCS02-mismatch';
cfg.stage_dates         = stage_dates{2}; % starts at Stage 1

      plot_timeline(cfg, REDcap);

%{
Takeaways:
    * VAS 50s dominate, but is otherwise bimodal
    * for most pts, MPQ affective subscore is uninformative
    * NRS has normal(ish) distribution
%}

%% Behavioral clustering: identify low/high pain states + reduce dimensionality
%{
* filter-out 50s (a neccesary alteration to identify natural pain subspaces)

* z-score metrics to themselves for easier comparison of magnitude btwn
    metrics

* cluster based off local density while prioritizing distance
    btwn high density points (Rodriguez & Laio, 2014, Science)

* determine first PC to have summary predictive metric for neural analyses

    * assess/finalize stability of inputs to PC and outputs

%}
pts = {'RCS02', 'RCS04', 'RCS05', 'RCS06', 'RCS07'};

pain_space      = [];

cfg             = [];     
cfg.dates       = 'AllTime';
cfg.pca         = false;
cfg.plt_VAS     = true;
cfg.VAS_only    = false;

cfg.CBDP_method = 'top_two';

cfg.source_dir  = ['/Users/Leriche/',...
                   'Dropbox (UCSF Department of Neurological Surgery)/',...
                   'UFlorida_UCSF_RCS_collab/Pain Reports/beh_clustered/'];

cfg.fig_dir     = [github_dir, 'Analysis-rcs-data/working/plot_beh/figs/beh_only/stages1_2_3/'];

set(0,'DefaultFigureVisible','off')

for i =  1:length(pts)

    cfg.pt_id  = pts{i};

    [pain_space.(pts{i})] = plot_pain_space(cfg, REDcap);
end

set(0,'DefaultFigureVisible','on')
%%
%
%% ******************** in-progress versions below ******************** %%
%
%%
% last N days for: 
cfg                     = [];

cfg.pt_id               = 'RCS04';
cfg.dates               = 'PreviousDays';
cfg.ndays               = 5;

cfg.subplot             = true;
cfg.sum_stat_txt        = true;
cfg.stim_parameter      = '';


    plot_timeline(cfg, REDcap);


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

% visualize distribution of PC 1--wait for now--could be useful for distributed closed-loop

%{

figure('Units', 'Inches', 'Position', [0, 0, 15, 10])

sgtitle('Dist of PC 1 per patient')


subplot(311)
title('RCS02');
hold on;
histogram(beh_anl.RCS02_pc.score(:,1)) 

set(gca,'fontSize',14, 'TickLength', [0 0]); 
        grid on;    grid MINOR;      box off


subplot(312)
title('RCS04')
hold on
histogram(beh_anl.RCS04_pc.score(:,1)) 

set(gca,'fontSize',14, 'TickLength', [0 0]); 
        grid on;    grid MINOR;      box off

subplot(313)
title('RCS05');
hold on;
histogram(beh_anl.RCS05_pc.score(:,1)) 

xlabel('PC 1')

 set(gca,'fontSize',14, 'TickLength', [0 0]); 
        grid on;    grid MINOR;      box off
%}

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

%% RCS04 Stim Plan

% 07/13/22; Home Programs
%{
LCaud
    open-loop
        A:     9+11-    1 mA     125 Hz     300 mcs     60s/20s
        B:     C+10-    2 mA     125 Hz     300 mcs     60s/20s  
        C:     C+9-     1 mA     125 Hz     300 mcs     60s/20s

    closed-loop
        D:     C+10-    0-1 mA   100 Hz     300 mcs     aDBS

RThal
    open-loop
        A:     9+11-    1 mA     150 Hz     200 mcs     60s/20s
        B:     C+11-    3 mA     100 Hz     200 mcs     60s/20s

    closed-loop
        D:     C+10-    0-1 mA   100 Hz     300 mcs     aDBS

RACC
    open-loop
        C:     0+3-     2 mA     100 Hz     300 mcs     60s/20sk

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
