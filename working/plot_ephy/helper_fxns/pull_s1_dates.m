function [s1_start, s1_end] ...
          = pull_s1_dates(parsed_db_in)


% ensure that every session has therapy Off or at least at 0 mA
amp_oi = nan(height(parsed_db_in), 1);

parsed_db_in(cellfun(@isempty, parsed_db_in.activeGroup),...
             'activeGroup') = {'A'};

for i_sess = 1 :height(parsed_db_in)
    amp_oi(i_sess) =  parsed_db_in{i_sess, ...
                                   sprintf('Group%sProg0_ampInMilliamps', parsed_db_in.activeGroup{i_sess})...
                                    };
end

i_on     = find(parsed_db_in.therapyStatus & amp_oi > 0 & amp_oi ~= 8.5);


% return FIRST session ever - > session PRIOR to first session with sitm ON above 0 mA
s1_start = dateshift(parsed_db_in.timeStart((1)), 'start', 'day');
  
post_s1      = s1_start + duration('120:00:00');
i_on_peri_s1 = find(le(parsed_db_in.timeStart(i_on), post_s1));

if ~isempty(i_on_peri_s1)
    i_on_at_home = i_on_peri_s1(end) +1;

else
    i_on_at_home =1;
end

s1_end   = dateshift(parsed_db_in.timeStart(i_on(i_on_at_home)-1), 'start', 'day');
end