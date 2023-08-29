function [beh_RCSXX] = db_sort_beh(cfg, db_RCSXXL, db_RCSXXR, redcap, visits_tbl  )
warning("off", "all"); 
% cfg                 = [];
% cfg.pt_id           = 'RCS05';
% cfg.pia_raw_dir     = pia_raw_dir;
% 
% cfg.plot_sess_dur    = false;
% 
% 
% db_RCSXXL           = db.RCS05L;
% db_RCSXXR           = db.RCS05R;
% redcap              = REDcap.RCS05;
% visits_tbl          = visits.RCS05;

%% handle unilateral implants
if isempty(db_RCSXXR)
    beh_RCSXX       = sortrows(db_parse(db_RCSXXL), 'timeStart');
end

if isempty(db_RCSXXL)
    beh_RCSXX       = sortrows(db_parse(db_RCSXXR), 'timeStart');
end

%% aligning bilateral implants
if ~isempty(db_RCSXXL) && ~isempty(db_RCSXXR)

    beh_RCSXXL      = db_parse(db_RCSXXL);
    beh_RCSXXR      = db_parse(db_RCSXXR);
    beh_RCSXX       = sortrows([beh_RCSXXL; beh_RCSXXR], 'timeStart');
end

% look at before and after sess mean

%% finding sessions w/ Stim On

beh_RCSXX.stimRegOn = ...
                repmat({''}, size(beh_RCSXX.stimContacts));

beh_RCSXX.stimRegOn(beh_RCSXX.stimAmp > 0) = ...
    beh_RCSXX.stimContacts(beh_RCSXX.stimAmp > 0);

beh_RCSXX           = beh_RCSXX(ge(beh_RCSXX.duration, '00:03:00'),:);

for j = 1 : height(beh_RCSXX) - 1

   if   le(beh_RCSXX.timeStart(j+1) - beh_RCSXX.timeStart(j), '15:00:00')...
        &&...
            strcmp(beh_RCSXX.stimContacts(j), beh_RCSXX.stimContacts(j+1)) == 0 ...
        &&...
            strcmp(beh_RCSXX.stimSide(j), beh_RCSXX.stimSide(j+1)) == 0 ...
        &&...
            (beh_RCSXX.stimAmp(j) > 0 && beh_RCSXX.stimAmp(j+1) > 0) ...

            
       j_count   = j + 1;


       while j_count < height(beh_RCSXX)

           sorted_str              = sort({beh_RCSXX.stimContacts{j}, ...
                                         beh_RCSXX.stimContacts{j_count}});
    
           beh_RCSXX.stimRegOn{j}  = [sorted_str{1}, ' & ', sorted_str{2}];

           beh_RCSXX.stimRegOn{j_count} = [sorted_str{1}, ' & ', sorted_str{2}];

           if   le(beh_RCSXX.timeStart(j_count) - beh_RCSXX.timeStart(j), '15:00:00')...
                &&...
                    strcmp(beh_RCSXX.stimContacts(j), beh_RCSXX.stimContacts(j_count)) == 0 ...
                &&...
                    strcmp(beh_RCSXX.stimSide(j), beh_RCSXX.stimSide(j_count)) == 0 ...
                &&...
                    beh_RCSXX.stimAmp(j_count) > 0
               
              break


           end

           j_count = j_count + 1;

       end
   end
end

%% visually inspect for edge cases 
if strcmp(cfg.pt_id, 'RCS04')

% when two of the same stim program
% were recorded conseqitively that inital session is falsely "unilateral"
    beh_RCSXX.stimRegOn{404} = 'L Caud: c+11- & R ACC: 0+3-';

    beh_RCSXX.stimRegOn{565} = 'L Caud: c+10- & R THAL: 9+11-';

    beh_RCSXX.stimRegOn([687, 717,801]) = {'L Caud: 9+11- & R THAL: 9+11-'};

% unresolved edge cases
%{
    550, 645 - day before bilateral stim; appears LCaud unilateral was tried
    683, 684 - appears unilateral, but no L session was streamed at all

Conclusion:
    Alignment okay -> use/visualize 'EventLog.txt' to index stim parameters from
    in-btwn sessions and/or during unilateral streaming w/ possible bilateral
    stim.

    so small it's likely moot--PS expects <10%

    pain reports tht are explicitly w/n streaming session

%}

elseif strcmp(cfg.pt_id, 'RCS05')
% visually inspected for edges cases--when two of the same stim program
% were recorded conseqitively that inital session is falsely "unilateral"

    % beh_RCSXX.stimRegOn{} = 

end


% visualizing session durations
if cfg.plot_sess_dur

    prop_ge_10 = length(find(ge(beh_RCSXX.duration, '00:10:00'))) ./ height(beh_RCSXX);
    
    prop_ge_30 =length(find(ge(beh_RCSXX.duration, '00:30:00'))) ./ height(beh_RCSXX);
    
    prop_ge_60 = length(find(ge(beh_RCSXX.duration, '01:00:00'))) ./ height(beh_RCSXX);
    
    
    figure('Units', 'Inches', 'Position', [0, 0, 12, 7])

    histogram(beh_RCSXX.duration,'BinWidth', duration('00:10:00')); xlim(duration({'00:00:00', '03:00:00'}))
    
    title([cfg.pt_id, newline, 'Session Durations'], 'Fontsize',16);
    
    TextLocation(...
        [ ...
            'p(>10 min) : ', num2str(prop_ge_10),...
            newline,...
            'p(>30 min) : ', num2str(prop_ge_30),...
                newline,...
            'p(>60 min) : ', num2str(prop_ge_60) ...
        ],...
            'Location','best');
    
    set(gca,'FontSize',12, 'TickLength', [0 0]); 
end 

%% sanity check of REDcap surverys w/ API logs
%{
api_log_root        = [cfg.pia_raw_dir,'/',cfg.pt_id,'/logs/'];
api_log_files       = dir([api_log_root,'*.log']);

logs = table();

for h = 1: height(api_log_files)

    api_log = readtable([api_log_root, api_log_files(h).name]);

    if size(api_log,2) > 2
 
        n_column = size(api_log,2);
        event = ...
             table2array(mergevars(api_log(:,3: n_column), ...
                         api_log.Properties.VariableNames(3:n_column)));

         api_log = api_log(:, 1:3);
         
         for j = 1: height(event)

            api_log{j,3} = {strjoin(event(j,:))};

         end
    end

    logs.api_log_date{h} = api_log_files(h).name(5:end-4);
    logs.data{h} = api_log;

    if ~isempty(logs.data{h}) && iscellstr(table2array(logs.data{h}(:,2))) %#ok<ISCLSTR> 
    
        i_invok_act = cellfun(@(x) contains(x,'Invoking Action: WebPage'),...
            table2array(logs.data{h}(:,2)));
        
        logs.invoked_act_name{h} = table2array(logs.data{h}(i_invok_act, 2));
        logs.invoked_act_time{h} = table2array(logs.data{h}(i_invok_act, 1));
        
    end
end

% filter out empty logs and then expand logs w/ multiple entries
i_log_act           = cellfun(@(x) ~isempty(x) , logs.invoked_act_name);
temp_WebPage        = logs(i_log_act, [1, 3, 4]);


api_log_WebPage                     = table();
                              
api_log_WebPage.api_log_date        = {''};
api_log_WebPage.invoked_act_name    = {''};
api_log_WebPage.invoked_act_time    = {''};

for i = 1 : height(temp_WebPage)

     tmp_row = temp_WebPage(i, :);

  if length(tmp_row.invoked_act_name{1}) > 1

    for new_row = 1 : length(temp_WebPage.invoked_act_name{i})

         
       api_log_WebPage{end+1,:} = {tmp_row.api_log_date,...
                                   tmp_row.invoked_act_name{1}{new_row},...
                                   tmp_row.invoked_act_time{1}(new_row)};
    end

  else
      api_log_WebPage = [api_log_WebPage; tmp_row];

  end
end

% pull datetime out of cell
% return only 'WebPageOne' (i.e., the "Pain Report" Link)
api_redcap   = api_log_WebPage(...
                cellfun(@(x) contains(x, 'WebPageOne'), api_log_WebPage.invoked_act_name),...
                              :);

api_redcap.invoked_act_time           = cellfun(@(x) x, api_redcap.invoked_act_time);
api_redcap.invoked_act_time.TimeZone  = 'America/Los_Angeles';

% with all WebPage entries expanded, find nearest redcap report and Session
% start and stops for time cross-validation



%%
%{
check for API logs that occured btwn timeStart and timeStop

if one exists, find nearest redcap report --> which ever API log is closest
in time to the redcap report (of course check performance afterwards) treat
that as the true time


see proportion of sessions w/ WebPage one pressed

%}

i_wn_sess = strcmp(beh_RCSXX.log_wn_sess, 'api_log_wn_sess') &...
                beh_RCSXX.timeStop > api_redcap.invoked_act_time(1);

le_10min = beh_RCSXX.redcap_lat(i_wn_sess & ...
                               le(beh_RCSXX.redcap_lat, '0:10:00'));


% visualize API Log to redcap Server Lag
figure('Units', 'Inches', 'Position', [0, 0, 12, 7])

beh_RCSXX           = beh_RCSXX(ge(beh_RCSXX.duration, '00:03:00'),:);


histogram(beh_RCSXX.redcap_lat(ge(beh_RCSXX.redcap_lat, '0:00:01')),...
    'BinWidth', duration({'0:15:00'})); 


title([cfg.pt_id, newline, ...
        'redcap Server to API Log latencies'], 'Fontsize',16);


t = TextLocation(...
    [...
    'prop(API Log "WebPage" w/n Session): ', num2str(sum(i_wn_sess) ./ sum(beh_RCSXX.timeStop > api_redcap.invoked_act_time(1))),...
    newline, newline...
    'prop(latency < 10 min): ', num2str(length(le_10min)./ sum(i_wn_sess)),...
    ],...
    'Location', 'Best');


t.FontSize = 12;

set(gca,'FontSize',12, 'TickLength', [0 0]); 
grid on;    grid MINOR;      box off;

% for RCS02 see prop of API log w/n sess where API logs exist

le_10min = beh_RCSXX.redcap_lat(i_wn_sess & ...
                               le(beh_RCSXX.redcap_lat, '0:10:00'));


% visualize API Log to redcap Server Lag
figure('Units', 'Inches', 'Position', [0, 0, 12, 7])

beh_RCSXX           = beh_RCSXX(ge(beh_RCSXX.duration, '00:03:00'),:);


histogram(beh_RCSXX.redcap_lat(ge(beh_RCSXX.redcap_lat, '0:00:01')),...
    'BinWidth', duration({'0:15:00'})); 


title([cfg.pt_id, newline, ...
        'redcap Server to API Log latencies'], 'Fontsize',16);


t = TextLocation(...
    [...
    'prop(API Log "WebPage" w/n Session): ', num2str(sum(i_wn_sess) ./ height(beh_RCSXX)),...
    newline, newline...
    'prop(latency < 10 min): ', num2str(length(le_10min)./ sum(i_wn_sess)),...
    ],...
   'Location', 'Best');


t.FontSize = 12;

set(gca,'FontSize',12, 'TickLength', [0 0]); 
grid on;    grid MINOR;      box off;

%}


%%
% visualze "nearest" redcap report


for j = 1 : height(beh_RCSXX)

    btwn_start_and_stop = ...
        find(ge(redcap.time, beh_RCSXX.timeStart(j)) & ...
        le(redcap.time, beh_RCSXX.timeStop(j)));


    if ~isempty(btwn_start_and_stop)

        beh_RCSXX.redcap_in_sess{j}     = 'wn_sess';
        beh_RCSXX.redcap_lat_wrt_sess_start{j}...
                                        = duration(redcap.time(btwn_start_and_stop)...
                                            - beh_RCSXX.timeStart(j));

        beh_RCSXX.redcap_lat_wrt_sess_stop{j}...
                                        = duration(beh_RCSXX.timeStop(j)...
                                            - redcap.time(btwn_start_and_stop));
    
        beh_RCSXX.redcap_reports{j}     = redcap(btwn_start_and_stop,:);
        beh_RCSXX.i_redcap(j)           = {btwn_start_and_stop};

     else

        beh_RCSXX.redcap_in_sess{j}     = 'not_wn_sess';

     end


    % repeat with reports *near* streaming sessions
    near_start_and_stop ...
        =find(...
            ge(redcap.time + duration('0:15:00'), beh_RCSXX.timeStart(j)) & ...
            le(redcap.time, beh_RCSXX.timeStop(j)));


    if ~isempty(near_start_and_stop)



    beh_RCSXX.redcap_near_sess{j}     = 'near_sess';
    beh_RCSXX.redcap_lat_near_sess_start{j}...
                                    = duration(redcap.time(near_start_and_stop)...
                                        - beh_RCSXX.timeStart(j));

    beh_RCSXX.redcap_lat_near_sess_stop{j}...
                                    = duration(beh_RCSXX.timeStop(j)...
                                        - redcap.time(near_start_and_stop));

    beh_RCSXX.redcap_near_reports{j}     = redcap(near_start_and_stop,:);
    beh_RCSXX.i_redcap_near(j)           = {near_start_and_stop};

    else

     beh_RCSXX.redcap_near_sess{j}     = 'no_nearby_sess';

    end  

% see if session start occured during inpatient, inclinic, or home testing
    sess_visit_diff =  beh_RCSXX.timeStart(j) - visits_tbl.dates;

    i_visit_day    = find(le(sess_visit_diff, duration('24:00:00')) &...
                          ge(sess_visit_diff, '0:00:00'));

    if ~isempty(i_visit_day)

        beh_RCSXX.visits(j)  = visits_tbl.desc(i_visit_day);
    
    else

        beh_RCSXX.visits(j)  = {''};

    end
    

end


i_1rep_wn_sess       = strcmp(beh_RCSXX.redcap_in_sess, 'wn_sess') &...
                        cellfun(@(x) length(x) == 1, beh_RCSXX.redcap_lat_wrt_sess_start);

i_many_wn_sess       = strcmp(beh_RCSXX.redcap_in_sess, 'wn_sess') &...
                        cellfun(@(x) length(x) > 1, beh_RCSXX.redcap_lat_wrt_sess_start);

i_not_wn_sess        = strcmp(beh_RCSXX.redcap_in_sess, 'not_wn_sess');

n_1rep_wn_sess       = sum(i_1rep_wn_sess);
n_many_wn_sess       = sum(i_many_wn_sess);
n_not_wn_sess        = sum(i_not_wn_sess);


prop_1rep_wn_sess      = n_1rep_wn_sess ./ height(beh_RCSXX);
prop_many_wn_sess      = n_many_wn_sess ./ height(beh_RCSXX);
prop_none_wn_sess      = n_not_wn_sess  ./ height(beh_RCSXX);

lat_sess_start = cellfun(@(x) x , beh_RCSXX.redcap_lat_wrt_sess_start(i_1rep_wn_sess));

i_1rep_near_sess       = strcmp(beh_RCSXX.redcap_near_sess, 'near_sess') &...
                             cellfun(@(x) length(x) == 1, beh_RCSXX.redcap_lat_near_sess_start);

i_many_near_sess       = strcmp(beh_RCSXX.redcap_near_sess, 'near_sess') &...
                             cellfun(@(x) length(x) > 1, beh_RCSXX.redcap_lat_near_sess_start);

i_no_near_sess        = strcmp(beh_RCSXX.redcap_near_sess, 'no_nearby_sess');

n_1rep_near_sess       = sum(i_1rep_near_sess);
n_many_near_sess       = sum(i_many_near_sess);
n_not_near_sess        = sum(i_no_near_sess);


prop_1rep_near_sess      = n_1rep_near_sess ./ height(beh_RCSXX);
prop_many_near_sess        = n_many_near_sess ./ height(beh_RCSXX);
prop_none_near_sess        = n_not_near_sess  ./ height(beh_RCSXX);


i_1rep_near_not_wn         = i_1rep_near_sess & i_not_wn_sess;
n_near_not_wn              = sum(i_1rep_near_not_wn);
prop_near_not_wn           = n_near_not_wn  ./ height(beh_RCSXX);

near_sess_start = cellfun(@(x) x , beh_RCSXX.redcap_lat_near_sess_start(i_1rep_near_not_wn));

near_sess_start = near_sess_start(...
                    le(near_sess_start, '0:00:00'));

%%
figure('Units', 'Inches', 'Position', [0, 0, 12, 7])


histogram(near_sess_start, 'BinWidth', duration({'0:01:00'}));
hold on;
histogram(lat_sess_start, 'BinWidth', duration({'0:01:00'}));

xlim(duration({'-0:30:00', '02:00:00'}));

title([cfg.pt_id, newline, ...
    'REDcap report "time" - RCS Sessision "timeStart"'], 'Fontsize',16);

TextLocation(...
    [...
    'prop( "early" reports ): ', num2str(prop_near_not_wn * 100), ' %', ...
    newline...
    'N reports = ', num2str(n_near_not_wn),...
    newline...
    newline,...
    'prop( w/n Sess ): ', num2str((prop_1rep_wn_sess + prop_many_wn_sess) * 100), ' %', ...
    newline,...
    'N reports = ', num2str(n_1rep_wn_sess+ n_many_wn_sess)...
     newline...
     newline,...
    'prop( w/n Sess + "early" reports ): ', ...
        num2str((prop_1rep_wn_sess + prop_many_wn_sess + prop_near_not_wn) * 100), ' %', ...
    newline,...
    'N reports = ', num2str(n_1rep_wn_sess+ n_many_wn_sess + n_near_not_wn)...
    ],...
    'Location', 'Best');

xlabel('Time "HH:mm:ss"'); ylabel('N of reports')

set(gca,'FontSize',12, 'TickLength', [0 0]); 
grid on;    grid MINOR;    box off;

%% explore so-called missing REDcap surveys
%{
X during visits
    - No, streaming and REDcap survery compliance is excellent during
      visits, but visits to sessions association has been made for RCS05.

* work with EventLogs
    
%}

% RCS05 prop(missing) = 0.1382
prop_missing = 1-(prop_1rep_wn_sess + prop_many_wn_sess + prop_near_not_wn);

i_all_wn_or_near     = i_1rep_wn_sess     | i_many_wn_sess |...
                       i_1rep_near_sess   | i_many_near_sess;


i_visits              = ~strcmp(beh_RCSXX.visits, '');

prop_visits           = sum(i_visits)./ height(beh_RCSXX);

i_miss_during_visits  = i_visits & ~i_all_wn_or_near;


% find(i_miss_during_visits)


assigned_reports = unique(vertcat(beh_RCSXX.i_redcap_near{:}));

%%
figure('Units', 'Inches', 'Position', [0, 0, 12, 7])


near_sess_stop = cellfun(@(x) x , beh_RCSXX.redcap_lat_near_sess_stop(i_1rep_near_not_wn));

histogram(near_sess_stop, 'BinWidth', duration({'0:01:00'}));
hold on

xlim(duration({'-0:20:00', '02:00:00'}));

title([cfg.pt_id, newline, ...
    'RCS Sessision "timeStop" - REDcap report "time"'], 'Fontsize',16);

%%

if sum([prop_1rep_wn_sess, prop_many_wn_sess, prop_none_wn_sess])...
        ~= 1  

   disp('error in indexing REDcap reports wrt Session Starts/Stops')
   disp('code execution paused')
   pause

end

% handling multiple pain reports w/n and/or "near" a session

%{
start_lat_wn    = [];

start_lat_near  = [];
for i = 1 : height(beh_RCSXX.redcap_in_sess)

     tmp_row = beh_RCSXX.redcap_lat_wrt_sess_start{i, :};

  if ~isempty(tmp_row)

    for new_row = 1 : length(temp_WebPage.invoked_act_name{i})

         
       api_log_WebPage{end+1,:} = {tmp_row.api_log_date,...
                                   tmp_row.invoked_act_name{1}{new_row},...
                                   tmp_row.invoked_act_time{1}(new_row)};
    end

  else
      api_log_WebPage = [api_log_WebPage; tmp_row];

  end
end

start_lat       = [beh_RCSXX.redcap_lat_wrt_sess_start{i_wn_sess}]';

stop_lat        = [beh_RCSXX.redcap_lat_wrt_sess_stop{i_wn_sess}]';


%}
%%

% plotting code that does not consider absolute session duration when
% comparing time difference with redcap report

%{
figure; histogram(beh_RCSXX.redcap_near_stop, 'BinWidth', duration({'0:15:00'})); 
xlim(duration({'0:00:00', '06:00:00'}));

title([cfg.pt_id, newline, 'Difference btwn. redcap and mean Session Duration'], 'Fontsize',16);

text(duration('02:00:00'), height(beh_RCSXX)/6,...
    ['p(<30 min) : ', num2str(length(find(le(beh_RCSXX.redcap_near_stop, '00:30:00'))) ./ length(beh_RCSXX.redcap_near_stop)),...
    newline,...
    'p(<60 min) : ', num2str(length(find(le(beh_RCSXX.redcap_near_stop, '01:00:00'))) ./ length(beh_RCSXX.redcap_near_stop)),...
        newline,...
    'p(<90 min) : ', num2str(length(find(le(beh_RCSXX.redcap_near_stop, '01:30:00'))) ./ length(beh_RCSXX.redcap_near_stop))],...
    'FontSize', 14);

%}
% see proportion of time diff btw. redcap report and Session time for
% nearest report

%{
SD_ge_15_le_3hrs = beh_RCSXX.redcap_time_diff((...
    ge(beh_RCSXX.duration, '00:15:00') & le(beh_RCSXX.duration, '03:00:00')));

SD_ge_30_le_3hrs = beh_RCSXX.redcap_time_diff((...
    ge(beh_RCSXX.duration, '00:30:00') & le(beh_RCSXX.duration, '03:00:00')));




figure('Units', 'Inches', 'Position', [0, 0, 12, 7])


histogram(SD_ge_15_le_3hrs, 'BinWidth', duration({'0:10:00'})); 
hold on
histogram(SD_ge_30_le_3hrs, 'BinWidth', duration({'0:10:00'})); 

xlim(duration({'0:00:00', '06:00:00'}));

title([cfg.pt_id, newline, ...
    'Difference btwn. redcap and mean Session Duration'], 'Fontsize',16);

legend({['3hrs > SD > 15 min;   N ', num2str(length(SD_ge_15_le_3hrs))],...
        ['3hrs > SD > 30 min;   N ', num2str(length(SD_ge_30_le_3hrs))]})

grid on;    grid MINOR;   legend boxoff;    box off;


t = TextLocation(...
    [...
    'prop(3hrs > SD > 15 min) w/n 30 min of redcap : ', num2str(sum(le(SD_ge_15_le_3hrs, '00:30:00')) ./ length(SD_ge_15_le_3hrs)),...
    newline,...
    'prop(3hrs > SD > 15 min) w/n 60 min of redcap : ', num2str(sum(le(SD_ge_15_le_3hrs, '01:00:00')) ./ length(SD_ge_15_le_3hrs)),...
    newline, newline,...
    'prop(3hrs > SD > 30 min) w/n 30 min of redcap : ', num2str(sum(le(SD_ge_30_le_3hrs, '00:30:00')) ./ length(SD_ge_30_le_3hrs)),...
    newline,...
    'prop(3hrs > SD > 30 min) w/n 60 min of redcap : ', num2str(sum(le(SD_ge_30_le_3hrs, '01:00:00')) ./ length(SD_ge_30_le_3hrs))...
    ],...
   'Location', 'Best');
t.FontSize = 12;

set(gca,'FontSize',12, 'TickLength', [0 0]); 
%}
%}
%% local fxns

function beh_RCSXXX = db_parse(db_RCSXXX)
% removes any empty or more than 1 timestamps per session
db_RCSXXX               = db_RCSXXX(cellfun(@(x) length(x)==1, db_RCSXXX.timeStart),:);

beh_RCSXXX              = table;

for i_db = 1:height(db_RCSXXX)

    beh_RCSXXX.timeStart(i_db)           = db_RCSXXX.timeStart{i_db};
    
    beh_RCSXXX.timeStop(i_db)            = db_RCSXXX.timeStop{i_db};
    
    beh_RCSXXX.duration(i_db)            = db_RCSXXX.duration{i_db};

    beh_RCSXXX.rec(i_db)                  = db_RCSXXX.rec(i_db);

    if height(db_RCSXXX.stimLogSettings{i_db}) == 1 

        beh_RCSXXX.stimSide{i_db}         = db_RCSXXX.stimReg{i_db,1};

     
        stim_para = strsplit(char(db_RCSXXX.stimLogSettings{i_db}.stimParams_prog1),',');

        beh_RCSXXX.stimContacts{i_db}     = [db_RCSXXX.stimReg{i_db}, ': ', stim_para{1}];
        
        beh_RCSXXX.stimAmp(i_db)          = str2double(stim_para{2}(2:end -2));
        beh_RCSXXX.stimPW(i_db)           = str2double(stim_para{3}(2:end -2));
        beh_RCSXXX.stimfreq(i_db)         = str2double(stim_para{4}(2:end -2));

        beh_RCSXXX.stimGroup{i_db}        = db_RCSXXX.stimSettingsOut{i_db}.activeGroup{1};

        beh_RCSXXX.stimCycleOnTime(i_db)     = db_RCSXXX.stimSettingsOut{i_db}.cycleOnTime{1};
        beh_RCSXXX.stimCycleOffTime(i_db)    = db_RCSXXX.stimSettingsOut{i_db}.cycleOffTime{1};
    

    else

        beh_RCSXXX.stimSide{i_db}              = '';
        beh_RCSXXX.stimContacts{i_db}          = '';

      
        beh_RCSXXX.stimAmp(i_db)               = NaN;
        beh_RCSXXX.stimPW(i_db)                = NaN;
        beh_RCSXXX.stimfreq(i_db)              = NaN;

        beh_RCSXXX.stimCycleOnTime(i_db)       = NaN;
        beh_RCSXXX.stimCycleOffTime(i_db)      = NaN;
        

        beh_RCSXXX.stimGroup{i_db}             = '';
        

    end
end
  
beh_RCSXXX              = table2timetable(beh_RCSXXX);
beh_RCSXXX              = movevars(beh_RCSXXX, {'timeStop', 'duration'},'Before', 'rec');

end
end