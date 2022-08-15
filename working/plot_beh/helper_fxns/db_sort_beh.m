function [beh_RCSXX, api_redcap] = db_sort_beh(cfg, db_RCSXXL, db_RCSXXR, redcap)

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

% beh_RCSXX           = beh_RCSXX(ge(beh_RCSXX.duration, '00:03:00'),:);

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
prop_ge_10 = length(find(ge(beh_RCSXX.duration, '00:10:00'))) ./ height(beh_RCSXX);

prop_ge_30 =length(find(ge(beh_RCSXX.duration, '00:30:00'))) ./ height(beh_RCSXX);

prop_ge_60 = length(find(ge(beh_RCSXX.duration, '01:00:00'))) ./ height(beh_RCSXX);


figure('Units', 'Inches', 'Position', [0, 0, 12, 7])
histogram(beh_RCSXX.duration,'BinWidth', duration('00:10:00')); xlim(duration({'00:00:00', '03:00:00'}))

title([cfg.pt_id, newline, 'Session Durations'], 'Fontsize',16);

t = TextLocation(...
['p(>10 min) : ', num2str(prop_ge_10),...
newline,...
'p(>30 min) : ', num2str(prop_ge_30),...
    newline,...
'p(>60 min) : ', num2str(prop_ge_60)],...
'Location','best');

t.FontSize = 12;

set(gca,'FontSize',12, 'TickLength', [0 0]); 



%%

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

         
       api_log_WebPage{end+1,:} = [{tmp_row.api_log_date,...
                                    tmp_row.invoked_act_name{1}{new_row},...
                                    tmp_row.invoked_act_time{1}(new_row)}];
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

% with all WebPage entries expanded, find nearest REDcap report and Session
% start and stops for time cross-validation



%%
%{
check for API logs that occured btwn timeStart and timeStop

if one exists, find nearest REDcap report --> which ever API log is closest
in time to the REDcap report (of course check performance afterwards) treat
that as the true time


see proportion of sessions w/ WebPage one pressed

%}



for j = 1 : height(beh_RCSXX)

    btwn_start_and_stop = ...
        find(ge(api_redcap.invoked_act_time, beh_RCSXX.timeStart(j)) &...
        le(api_redcap.invoked_act_time, beh_RCSXX.timeStop(j)));
    
    if ~isempty(btwn_start_and_stop)

        beh_RCSXX.log_wn_sess{j} = 'api_log_wn_sess';

        wn_sess_api = api_redcap.invoked_act_time(btwn_start_and_stop);


        % min logic is diff for single log (i.e., scalar by vector
        % subtraction rather than vector by vector subtraction)
        if length(wn_sess_api) == 1

            beh_RCSXX.api_WebPage(j)  =  wn_sess_api;

            % min reflects nearest REDcap report
            [REDcap_lat, i_REDcap] = ...
            min(abs(redcap.time - wn_sess_api));
            
            
            beh_RCSXX.REDcap_lat(j)  =  REDcap_lat;
            beh_RCSXX.REDcap_ind(j)  =  i_REDcap;


        % finding nearest log to REDcap for many logs w/n sess
        else
        
            % min reflects nearest api log to column-wise REDcap reports
            [REDcap_lat, i_log] = ...
                min(abs(redcap.time' - wn_sess_api));
    
    
            % return smallest diff btwn REDcap and logs w/n session
            beh_RCSXX.REDcap_lat(j)  =  min(REDcap_lat); 
            
    
            % return timestamp of nearest log to REDcap report w/n session
            i_near_log                = i_log(min(REDcap_lat) == REDcap_lat);
            beh_RCSXX.api_WebPage(j)  =  wn_sess_api(i_near_log(1)); % in-case duplicate log timestamps
            
            % min reflects nearest api log to column-wise REDcap reports
            [~, i_REDcap] = ...
                min(abs(wn_sess_api' - redcap.time));
    
            beh_RCSXX.REDcap_ind(j) = i_REDcap(1);

        end

    else
        beh_RCSXX.log_wn_sess{j} = 'no_api_log_wn_sess';

    end
end


i_wn_sess = strcmp(beh_RCSXX.log_wn_sess, 'api_log_wn_sess');

le_5min = beh_RCSXX.REDcap_lat(i_wn_sess & ...
                               le(beh_RCSXX.REDcap_lat, '0:10:00'));

ge_55min_le_65min = beh_RCSXX.REDcap_lat(i_wn_sess & ...
                               ge(beh_RCSXX.REDcap_lat, '0:50:00') &...
                               le(beh_RCSXX.REDcap_lat, '1:10:00'));

% visualize API Log to REDcap Server Lag
figure('Units', 'Inches', 'Position', [0, 0, 12, 7])

beh_RCSXX           = beh_RCSXX(ge(beh_RCSXX.duration, '00:03:00'),:);


histogram(beh_RCSXX.REDcap_lat(ge(beh_RCSXX.REDcap_lat, '0:00:01')),...
    'BinWidth', duration({'0:15:00'})); 


title([cfg.pt_id, newline, ...
        'REDcap Server to API Log latencies'], 'Fontsize',16);


t = TextLocation(...
    [...
    'prop(API Log "WebPage" w/n Session): ', num2str(sum(i_wn_sess) ./ height(beh_RCSXX)),...
    newline, newline...
    'prop(latency < 5 min): ', num2str(length(le_5min)./ sum(i_wn_sess)),...
    newline,...
    'prop(55 min < latency < 65 min): ', num2str(length(ge_55min_le_65min)./ sum(i_wn_sess))],...
   'Location', 'Best');


t.FontSize = 12;

set(gca,'FontSize',12, 'TickLength', [0 0]); 
grid on;    grid MINOR;      box off;

%%
% mean session time as better estimate of pain report

% visualze "nearest" REDcap report

%{
for j = 1 : height(beh_RCSXX)

    btwn_start_and_stop = ...
        find(ge(REDcap.time, beh_RCSXX.timeStart(j)) & le(REDcap.time, beh_RCSXX.timeStop(j)));
    
    
    [time_diff, i_rep] = min(abs(beh_RCSXX.timeStop(j) - REDcap.time));
    
    beh_RCSXX.i_REDcap_near_stop(j)  = i_rep;

    beh_RCSXX.REDcap_near_stop(j)    = time_diff;




    if ~isempty(btwn_start_and_stop)

        beh_RCSXX.REDcap_in_sess{j}     = 'within_session';
        beh_RCSXX.REDcap_lat_wrt_sess_start{j}...
                                        = duration(REDcap.time(btwn_start_and_stop)...
                                            - beh_RCSXX.timeStart(j));

        beh_RCSXX.REDcap_lat_wrt_sess_stop{j}...
                                        = duration(beh_RCSXX.timeStop(j)...
                                            - REDcap.time(btwn_start_and_stop));

        
        beh_RCSXX.REDcap_reports{j}     = REDcap(btwn_start_and_stop,:);
        beh_RCSXX.i_REDcap(j)           = {btwn_start_and_stop};


    else

         beh_RCSXX.REDcap_in_sess{j}     = 'no_report_within_session';

    end
end

i_wn_sess       = strcmp(beh_RCSXX.REDcap_in_sess, 'within_session') &...
                    cellfun(@(x) length(x)==1, beh_RCSXX.REDcap_lat_wrt_sess_start);

n_wn_sess       = sum(i_wn_sess);

prop_wn_sess    = n_wn_sess ./ height(beh_RCSXX);

start_lat       = [beh_RCSXX.REDcap_lat_wrt_sess_start{i_wn_sess}]';

stop_lat        =  [beh_RCSXX.REDcap_lat_wrt_sess_stop{i_wn_sess}]';


alignment_check = sum(...
    beh_RCSXX.duration(i_wn_sess) == start_lat + stop_lat);

if alignment_check ~= n_wn_sess

    disp('possible error in REDcap to RCS session file alignment')
end


figure('Units', 'Inches', 'Position', [0, 0, 12, 7])


histogram(start_lat, 'BinWidth', duration({'0:10:00'})); 
hold on
histogram(stop_lat, 'BinWidth', duration({'0:10:00'})); 

xlim(duration({'0:00:00', '04:00:00'}));

title([cfg.pt_id, newline, ...
    'REDcap latencies'], 'Fontsize',16);

t = TextLocation(...
    [...
    'prop(REDcap w/n Start/Stop): ', num2str(prop_wn_sess)],...
   'Location', 'Best');

legend({'REDcap wrt Session Start',...
        'REDcap wrt Session Stop'})

t.FontSize = 12;

set(gca,'FontSize',12, 'TickLength', [0 0]); 
grid on;    grid MINOR;   legend boxoff;    box off;


% plotting code that does not consider absolute session duration when
% comparing time difference with REDcap report


figure; histogram(beh_RCSXX.REDcap_near_stop, 'BinWidth', duration({'0:15:00'})); 
xlim(duration({'0:00:00', '06:00:00'}));

title([cfg.pt_id, newline, 'Difference btwn. REDcap and mean Session Duration'], 'Fontsize',16);

text(duration('02:00:00'), height(beh_RCSXX)/6,...
    ['p(<30 min) : ', num2str(length(find(le(beh_RCSXX.REDcap_near_stop, '00:30:00'))) ./ length(beh_RCSXX.REDcap_near_stop)),...
    newline,...
    'p(<60 min) : ', num2str(length(find(le(beh_RCSXX.REDcap_near_stop, '01:00:00'))) ./ length(beh_RCSXX.REDcap_near_stop)),...
        newline,...
    'p(<90 min) : ', num2str(length(find(le(beh_RCSXX.REDcap_near_stop, '01:30:00'))) ./ length(beh_RCSXX.REDcap_near_stop))],...
    'FontSize', 14);


% see proportion of time diff btw. REDcap report and Session time for
% nearest report

%{
SD_ge_15_le_3hrs = beh_RCSXX.REDcap_time_diff((...
    ge(beh_RCSXX.duration, '00:15:00') & le(beh_RCSXX.duration, '03:00:00')));

SD_ge_30_le_3hrs = beh_RCSXX.REDcap_time_diff((...
    ge(beh_RCSXX.duration, '00:30:00') & le(beh_RCSXX.duration, '03:00:00')));




figure('Units', 'Inches', 'Position', [0, 0, 12, 7])


histogram(SD_ge_15_le_3hrs, 'BinWidth', duration({'0:10:00'})); 
hold on
histogram(SD_ge_30_le_3hrs, 'BinWidth', duration({'0:10:00'})); 

xlim(duration({'0:00:00', '06:00:00'}));

title([cfg.pt_id, newline, ...
    'Difference btwn. REDcap and mean Session Duration'], 'Fontsize',16);

legend({['3hrs > SD > 15 min;   N ', num2str(length(SD_ge_15_le_3hrs))],...
        ['3hrs > SD > 30 min;   N ', num2str(length(SD_ge_30_le_3hrs))]})

grid on;    grid MINOR;   legend boxoff;    box off;


t = TextLocation(...
    [...
    'prop(3hrs > SD > 15 min) w/n 30 min of REDcap : ', num2str(sum(le(SD_ge_15_le_3hrs, '00:30:00')) ./ length(SD_ge_15_le_3hrs)),...
    newline,...
    'prop(3hrs > SD > 15 min) w/n 60 min of REDcap : ', num2str(sum(le(SD_ge_15_le_3hrs, '01:00:00')) ./ length(SD_ge_15_le_3hrs)),...
    newline, newline,...
    'prop(3hrs > SD > 30 min) w/n 30 min of REDcap : ', num2str(sum(le(SD_ge_30_le_3hrs, '00:30:00')) ./ length(SD_ge_30_le_3hrs)),...
    newline,...
    'prop(3hrs > SD > 30 min) w/n 60 min of REDcap : ', num2str(sum(le(SD_ge_30_le_3hrs, '01:00:00')) ./ length(SD_ge_30_le_3hrs))...
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

        beh_RCSXXX.stimSide{i_db}         = db_RCSXXX.stimReg{i_db,1}(1);

     
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