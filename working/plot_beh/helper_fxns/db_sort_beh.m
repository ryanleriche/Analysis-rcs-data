function beh_RCSXX = db_sort_beh(cfg, db_RCSXXL, db_RCSXXR, REDcap)

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


% mean session time as better estimate of pain report
sess_time_mean = mean([beh_RCSXX.timeStart, beh_RCSXX.timeStop] ,2);

% visualze "nearest" REDcap report
for j = 1 : height(beh_RCSXX)

    diff_vec = REDcap.time - beh_RCSXX.timeStart(j); % + duration('0:05:00');

    diff_vec = diff_vec(ge(diff_vec, 0));

    [time_diff, j_near]  = min(diff_vec);

    nearest_REDcap_time  = REDcap.time(find(diff_vec == time_diff));

    if ge(nearest_REDcap_time, beh_RCSXX.timeStop(j))
        

    end
    % beh_RCSXX.REDcap{j}             = REDcap(j_near,:);
    beh_RCSXX.REDcap_after_start(j)   = time_diff;
    beh_RCSXX.i_REDcap(j)           = j_near;


end
% plotting code that does not consider absolute session duration when
% comparing time difference with REDcap report
%{
figure; histogram(beh_RCSXX.REDcap_time_diff, 'BinWidth', duration({'0:15:00'})); 
xlim(duration({'0:00:00', '06:00:00'}));

title([cfg.pt_id, newline, 'Difference btwn. REDcap and mean Session Duration'], 'Fontsize',16);

text(duration('02:00:00'), height(beh_RCSXX)/6,...
    ['p(<30 min) : ', num2str(length(find(le(beh_RCSXX.REDcap_time_diff, '00:30:00'))) ./ length(beh_RCSXX.REDcap_time_diff)),...
    newline,...
    'p(<60 min) : ', num2str(length(find(le(beh_RCSXX.REDcap_time_diff, '01:00:00'))) ./ length(beh_RCSXX.REDcap_time_diff)),...
        newline,...
    'p(<90 min) : ', num2str(length(find(le(beh_RCSXX.REDcap_time_diff, '01:30:00'))) ./ length(beh_RCSXX.REDcap_time_diff))],...
    'FontSize', 14);
%}

% see proportion of time diff btw. REDcap report and Session time for
% sessions of different min duration

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


%% local fxns


function beh_RCSXXX = db_parse(db_RCSXXX)
% removes any empty or more than 1 timestamps per session
db_RCSXXX               = db_RCSXXX(cellfun(@(x) length(x)==1,db_RCSXXX.timeStart),:);

beh_RCSXXX              = table;

for i = 1:height(db_RCSXXX)

    beh_RCSXXX.timeStart(i)           = db_RCSXXX.timeStart{i};
    
    beh_RCSXXX.timeStop(i)            = db_RCSXXX.timeStop{i};
    
    beh_RCSXXX.duration(i)            = db_RCSXXX.duration{i};

    beh_RCSXXX.rec(i)                  = db_RCSXXX.rec(i);

    if height(db_RCSXXX.stimLogSettings{i}) == 1

        beh_RCSXXX.stimSide{i}         = db_RCSXXX.stimReg{i,1}(1);

     
        stim_para = strsplit(char(db_RCSXXX.stimLogSettings{i}.stimParams_prog1),',');

        beh_RCSXXX.stimContacts{i}     = [db_RCSXXX.stimReg{i}, ': ', stim_para{1}];
        
        beh_RCSXXX.stimAmp(i)          = str2double(stim_para{2}(2:end -2));
        beh_RCSXXX.stimPW(i)           = str2double(stim_para{3}(2:end -2));
        beh_RCSXXX.stimfreq(i)         = str2double(stim_para{4}(2:end -2));

        beh_RCSXXX.stimGroup{i}        = db_RCSXXX.stimSettingsOut{i}.activeGroup{1};

        beh_RCSXXX.stimCycleOnTime(i)     = db_RCSXXX.stimSettingsOut{i}.cycleOnTime{1};
        beh_RCSXXX.stimCycleOffTime(i)    = db_RCSXXX.stimSettingsOut{i}.cycleOffTime{1};
    

    else

        beh_RCSXXX.stimSide{i}              = '';
        beh_RCSXXX.stimContacts{i}          = '';

      
        beh_RCSXXX.stimAmp(i)               = NaN;
        beh_RCSXXX.stimPW(i)                = NaN;
        beh_RCSXXX.stimfreq(i)              = NaN;

        beh_RCSXXX.stimCycleOnTime(i)       = NaN;
        beh_RCSXXX.stimCycleOffTime(i)      = NaN;
        

        beh_RCSXXX.stimGroup{i}             = '';
        

    end
end
  
beh_RCSXXX              = table2timetable(beh_RCSXXX);
beh_RCSXXX              = movevars(beh_RCSXXX, {'timeStop', 'duration'},'Before', 'rec');

end
end