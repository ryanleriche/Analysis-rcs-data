%% user-inputs

% where RCS files are saved from PIA server
pia_raw_dir               = '/Users/Leriche/pia_server/datastore_spirit/human/rcs_chronic_pain/rcs_device_data/raw/';

% where 'ryanleriche/Analysis-rcs-data' Github repo is saved locally
github_dir            = '/Users/Leriche/Github/Analysis-rcs-data/';

% application programming interface (API) token which is essentially a
% password to access REDcap remotely, and is unique per researcher per
% study (e.g., Ryan has a unique token for the RCS and PCS studies)

API_token             = '95FDE91411C10BF91FD77328169F7E1B';


% pulls/organizes arms from REDcap (go into fxn to add new arms)
cd([github_dir, 'working']);         addpath(genpath(github_dir));

REDcap                 = RCS_redcap_painscores(API_token);

% Stage start dates I-III for RCS pts 1-6
stage_dates             = {{''}, {'08-Sep-2020'; '31-Jan-2021'; '31-May-2022'},... % RCS02
                           {''}, {'13-May-2021'; '12-Jul-2021'},... % RCS04
                           {'21-Jul-2021';'07-Sep-2021'},... % RCS05
                           {'07-Jun-2022'}}; % RCS06

%% import RCS Databases per pt side

% save RCS session summaries as DataBases w/n 'db' structure and
% badsessions w/n 'bs' structure 

cfg                    = [];
cfg.load_EventLog      = false;
cfg.ignoreold          = false;



% RCS04
[db.RCS04L, bs.RCS04L] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS04L', cfg);

[db.RCS04R, bs.RCS04R] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS04R', cfg);

% RCS05
[db.RCS05L, bs.RCS05L] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS05L',cfg);

[db.RCS05R, bs.RCS05R] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS05R',cfg);

% RCS02
[db.RCS02R, bs.RCS02R] = ...
    makeDatabaseRCS_Ryan(pia_raw_dir, 'RCS02R',cfg);

%% import RCS logs (state changes from INS)
cfg = [];
cfg.pull_adpt_logs      = false;
cfg.pull_event_logs     = false;

[textlog.RCS04L]         = RCS_logs(pia_raw_dir,'RCS04L', cfg);

[textlog.RCS04R]         = RCS_logs(pia_raw_dir,'RCS04R', cfg);


[textlog.RCS05R]         = RCS_logs(pia_raw_dir,'RCS05R', cfg);

%% 

%{
determine when w/n streaming session REDcap survey occured

    using RCS data start and stop times

    go into 'read_adaptive_txt_log' and further deeper fxns for all needed
    fields

    
    
%}



%%
%{
* creates db_beh.RCSXX with timestamps, contacts, amp, PW, freq, cycling,
  group, and laterality (unilateral vs bilateral stim)

    -manually inspect laterality to handle edge cases

* align nearest REDcap report as row in db_beh.RCSXX 

for RCS05 and RCS02 especially look at pain reports within session 5 - 15 minutes
   

%}


cfg       = [];
cfg.pt_id = 'RCS04';

    db_beh.RCS04 = db_sort_beh(cfg, db.RCS04L, db.RCS04R, REDcap.RCS04);


cfg.pt_id = 'RCS05';
    db_beh.RCS05 = db_sort_beh(cfg,db.RCS05L, db.RCS05R, REDcap.RCS05);


cfg.pt_id = 'RCS02';
    db_beh.RCS02 = db_sort_beh(cfg,[], db.RCS02R, REDcap.RCS02);

%%
%{

w/ aligned REDcap reports, first visualize the parameter space per each
"major contact pair" (best analgesic contact pairs as assessed clinically)

    * inspect prop(3hr > SD > 15 min) sessions per contact pair

visualized explored amp-PW-freq space per contact


output(s)
sII

What's are baseline measure to assess "50% reduction in pain" (yet to
determine metric, composite metric yet to ruled out/generated)?

-> Something in Stage 0?

%}
cfg       = [];
cfg.pt_id = 'RCS04';

beh_RCSXX = db_beh.RCS04;



i_kept_sess = le(beh_RCSXX.REDcap_time_diff, '01:00:00') & ...
                            ge(beh_RCSXX.duration, '00:15:00') & ...
                            le(beh_RCSXX.duration, '03:00:00');


beh_RCSXX.i_kept_sess = i_kept_sess;

stim_pairs = unique(beh_RCSXX.stimRegOn(i_kept_sess & beh_RCSXX.stimAmp > 0));

sII     = table();


for i = 1 : length(stim_pairs)

    sII.([stim_pairs{i}]) = {beh_RCSXX(...
                                    strcmp(beh_RCSXX.stimRegOn,stim_pairs{i}) &...
                                    i_kept_sess...
                                    ,:)};
end

% L Caud: 9+11-
figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
sgtitle([cfg.pt_id, ' L Caud: 9+11- Stim Space'], 'Fontsize',16)

subplot(211)
scatter3(sII.('L Caud: 9+11-'){1},...
    'stimAmp','stimPW','stimfreq'...
    ,'filled');

subplot(212)
histogram(sII.('L Caud: 9+11-'){1}.stimAmp, 'BinWidth',0.025)
title('Dist. of Stim Amp (mA)')


% R THAL: 9+11-
figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
sgtitle([cfg.pt_id, ' R THAL: 9+11- Stim Space'], 'Fontsize',16)


subplot(311)
scatter3(sII.('R THAL: 9+11-'){1},...
    'stimAmp','stimPW','stimfreq'...
    ,'filled');

subplot(312)
histogram(sII.('R THAL: 9+11-'){1}.stimAmp, 'BinWidth',0.025)
title('Dist. of Stim Amp (mA)')


subplot(313)
histogram(sII.('R THAL: 9+11-'){1}.stimPW, 'BinWidth',10)
title('Dist. of Stim PW (/mu/s)')


% bilateral stim needs more fine grained alignment given that parameter
% space if doubled (needing to pair session files to respective pain
% reports)
%{ 

figure('Units', 'Inches', 'Position', [0, 0, 15, 10])
sgtitle([cfg.pt_id, ' L Caud: 9+11- & R THAL: 9+11- Parameter Space'], 'Fontsize',16)

subplot(211)
scatter3(sII.('L Caud: 9+11- & R THAL: 9+11-'){1},...
    'stimAmp','stimPW','stimfreq'...
    ,'filled')

subplot(212)
swarmchart(sII.('L Caud: 9+11- & R THAL: 9+11-'){1},...
    'stimAmp','stimfreq'...
    ,'filled')
%}



%*************************************************************************
%% WAIT on pain + stim visualization until stats framework is finalized
%{

DEFINE baseline pain, to evaluate 50% decrease as primary end-point

RCS04 stats framework:

* filter the sessions from 10 min->3 hr 

* run Kolmogorov–Smirnov test (see's if data are normal--explore pain
  metric distribtion more generally)

* run Kruskal-Walis test for the pain metrics across each of the contact-pairs 
  (grouping all freq, amp, and PW parameters w/n a contact--lack statistical 
  power for freq/PW/amp)

* for LCaud, RThal, and bilateral stim contacts look run run a Kruskal-Walis 
  test across amplitude w/n contacts 


* follow-up with a Wilconxin signed rank test (?) for signifigant groups 




for i 1: height(
 eventLogTable               = addRowToTable(newEntry,eventLogTable); 
 cellfun(@x , sII.('L Caud: 9+11-'){1}.REDcap)
%}

% was hoping more detailed API was here:
%{
api_log_root        = [pia_raw_dir,'/','RCS04','/logs/'];
api_log_files       = dir([api_log_root,'*.log']);

for i = 1 : height(api_log_files)

    api_log = readtable([api_log_root, api_log_files(i).name]);

    api_log.Properties.VariableNames = {'time', 'event','d'};

    all_webpage = api_log(cellfun(@(x) contains(x,'WebPage') ,api_log.event),:)

end
%}




%% RCS04 plot daily metrics

% Specify, cfg before calling functions--see below for examples.
cfg                     = [];
cfg.pt_id               = 'RCS04';

cfg.dates               = 'AllTime';
cfg.stage_dates         = stage_dates{4}; % starts at Stage 1
cfg.subplot             = false;

cfg.stim_parameter      = '';
   
    plot_timeline(cfg, REDcap.RCS04, db_beh.RCS04);

%% last 7 days for all pts
cfg                     = [];
cfg.pt_id               = 'RCS04';
cfg.stage_dates         = stage_dates{4}; % starts at Stage 1
cfg.subplot             = true;

cfg.stim_parameter      = 'contacts';

cfg.dates               = 'PreviousDays';
cfg.ndays               = 7;
cfg.subplot             = true;

    plot_timeline(cfg, REDcap.RCS04, db_beh.RCS04);

cfg.pt_id               = 'RCS05';
cfg.stage_dates         = stage_dates{5}; % starts at Stage 1

        plot_timeline(cfg, REDcap.RCS05, db_beh.RCS05);


cfg.pt_id               = 'RCS02';
cfg.stage_dates         = stage_dates{2}; % starts at Stage 1

      plot_timeline(cfg, REDcap.RCS02, db_beh.RCS02);
      

% visually inspect pain metric distributions
cfg             = [];
cfg.pt_id       = 'RCS04';
cfg.dates       = 'AllTime';

    plot_hist(cfg, REDcap.RCS04);

cfg            = [];
cfg.dates      = 'AllTime';
% SEE Table for easy access to common summary statistics
RCS04_sum_stats      = calc_sum_stats(cfg, REDcap.RCS04);

% pain metric A "versus" pain metric B --> see how metrics visually covary
cfg             = [];
cfg.pt_id       = 'RCS04';
cfg.dates       = 'AllTime';

% go into this fxn to change the pain metrics on each axis

% RBL reccomends leaving NRS as the colorbar though, as it is a categorical
% metric--visually creating planes which muddles visualization
    plot_versus(cfg, REDcap.RCS04);



%% ************************************************************************
% clustering behavioral metrics to identify low/high pain states + reduce dimensionality


% build intuition of beh distribution via histograms
cfg             = [];
cfg.dates       = 'AllTime';
cfg.pt_id       = 'RCS02';       

plot_hist(cfg, REDcap.RCS02);

saveas(gcf, [cd, '/plot_beh/figs/RCS02/', 'RCS02_beh_hist.png']);

cfg.pt_id       = 'RCS04';       plot_hist(cfg, REDcap.RCS04);

saveas(gcf, [cd, '/plot_beh/figs/RCS04/', 'RCS04_beh_hist.png']);


cfg.pt_id       = 'RCS05';       plot_hist(cfg, REDcap.RCS05);

saveas(gcf, [cd, '/plot_beh/figs/RCS05/', 'RCS05_beh_hist.png']);


%{
Conclusion:
 VAS 50s dominate, but is otherwise bimodal, MPQ affective is uninformative, 
 and NRS has normal(ish) distribution

Next Steps:

* filter-out 50s (a neccesary alteration to identify natural pain subspaces)

* z-score metrics to themselves for easier comparison of magnitude btwn
    metrics

* cluster based off local density while prioritizing distance
    btwn high density points (Rodriguez & Laio, 2014, Science)

* determine first PC to have summary predictive metric for neural analyses

    * assess/finalize stability of inputs to PC and outputs

%}


cfg             = [];
cfg.pt_id       = 'RCS02';       
cfg.dates       = 'AllTime';
cfg.pca         = false;

[beh_fig.RCS02_filt, beh_fig.RCS02_z_space, beh_fig.RCS02_dec, beh_fig.RCS02_cl, ...
    beh_anl.z_RCS02] ...
    =...
    plot_versus(cfg, REDcap.RCS02);


% rinse & repeat for RCS04
cfg.pt_id       = 'RCS04'; 

[beh_fig.RCS04_filt, beh_fig.RCS04_z_space, beh_fig.RCS04_dec, beh_fig.RCS04_cl,...
    beh_anl.z_RCS04] ...
    =...
    plot_versus(cfg, REDcap.RCS04);


% rinse & repeat for RCS05
cfg.pt_id       = 'RCS05'; 

[beh_fig.RCS05_filt, beh_fig.RCS05_z_space, beh_fig.RCS05_dec, beh_fig.RCS05_cl, ...
    beh_anl.z_RCS05] ...
    =...
    plot_versus(cfg, REDcap.RCS05);

% rinse & repeat for RCS06
cfg.pt_id       = 'RCS06'; 

[beh_fig.RCS06_filt, beh_fig.RCS06_z_space, beh_fig.RCS06_dec, beh_fig.RCS06_cl, ...
    beh_anl.z_RCS06] ...
    =...
    plot_versus(cfg, REDcap.RCS06);

% save([cd,'/plot_beh/beh_anl'],'beh_anl');

save([cd,'/plot_beh/beh_anl_no_pca'],'beh_anl');


fig_names = fieldnames(beh_fig);
for i = 1 :  length(fig_names)

    savefig(beh_fig.(fig_names{i}),[cd,'/plot_beh/figs/', fig_names{i},'.fig'] )

    saveas(beh_fig.(fig_names{i}),...
        [cd,'/plot_beh/figs/', fig_names{i}, '.png']);

end


beh_anl_out     = fieldnames(beh_anl);

n_pts           = length(beh_anl_out);

for i = 1 :  n_pts

    writetable(beh_anl.(beh_anl_out{i}),...
        ['/Users/Leriche/Dropbox (UCSF Department of Neurological Surgery)/UFlorida_UCSF_RCS_collab/Pain Reports/beh_clustered/',...
             beh_anl_out{i}(3:end)], 'FileType','spreadsheet')

end


%%

%{

"ask question and it's opposite"

Chronback's coefficient alpha

%}

openfig([cd,'/plot_beh/figs/RCS02_z_space.fig'])

openfig([cd,'/plot_beh/figs/RCS04_z_space.fig'])

openfig([cd,'/plot_beh/figs/RCS05_z_space.fig'])

% visualize distribution of PC 1--wait for now--could be useful for
% distributed closed-loop

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

%% animate pain space visualization
%{

v = RecordRotation(h);

%%



playRecording(v);
%%
function V2 = RecordRotation(h)

n = 250;
V = zeros(n,2);
T = timer('period',0.05,'executionmode','fixedrate',...
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

%%

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


% Stim Comparisons
%{

07/15/22 - 07/29/22

    RACC
        open-loop
            C:     0+3-     2 mA     100 Hz     300 mcs     60s/20s

    *w/ option to decrease amplitude

    Q: If lowering the amplitude 0.1-0.5 mA does not help after 5-6 hours
    -->
        Lower amplitude further (min of 0.5 mA)?
        Keep waiting?
        Turn off stim?

    Q: Min N of days per contact-pair?

        Q: Within given contact-pair min N of days at given
           amplitude, pulse width, frequency?


Stim Sweeps as therapy (originally IJ's idea)


    
   analgesic effect of:

        *cumulative regular stim sweeps (daily, weekly, etc.)
    
            LCaudB:     C+10-    1 mA     125 Hz     300 mcs     60s/20s 


        "meditation-like" calming by focusing on internal state

            Q: Try double-blinded sham stim sweeps?

        

07/22/22, Stim Sweep
    
    (Program right C for reference)
    RACC
            C:     0+3-     2 mA     100 Hz     300 mcs     60s/20s
    -->
        freq:
        [80, 0, 90, 0, 100, 0, 110, 0, 100, 0, 90, 100, 0, 110]

        PW:
        [0, 200, 0, 250, 0, 300, 0, 200, 0, 300, 0, 250, 0]


        * 2 w/ PW & freq values shuffled while maintaining zeros btwn.








___________________________________________________________________________
bilateral v. unilateral stim

cumalative effect of stim

    regular stim-sweeps

    total energy delieverd as fxn of pain

%}

% Control Comparisons
%{

%}