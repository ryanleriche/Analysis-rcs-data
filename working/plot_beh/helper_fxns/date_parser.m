function [db_beh_RCSXX, RCSXX, date_range] = date_parser(cfg, RCSXX, db_beh_RCSXX)


    if strcmp(cfg.dates, 'AllTime') == 1

        date_range           = [RCSXX.time(1), RCSXX.time(end)];

    % takes current time and finds the nearest report 
    elseif strcmp(cfg.dates, 'PreviousDays') == 1
        
        [~, i_date_start] = min(abs(RCSXX.time ...
            - (datetime('now','TimeZone','America/Los_Angeles') - days(cfg.ndays))));

        [~, db_beh_date_start] = min(abs(...
                                mean([db_beh_RCSXX.timeStart,db_beh_RCSXX.timeStop],2) ...
                        - (datetime('now','TimeZone','America/Los_Angeles') - days(cfg.ndays))));


        
        date_range         = [RCSXX.time(i_date_start), RCSXX.time(end)];

        if length(unique(date_range)) == 1
            i_date_start = i_date_start - 1;
        end

        date_range         = [RCSXX.time(i_date_start), RCSXX.time(end)];

        RCSXX              = RCSXX(i_date_start:end, :);

        db_beh_RCSXX       = db_beh_RCSXX(db_beh_date_start: end, :);

     % takes 'cfg.date_range' and finds the nearest report 
     elseif strcmp(cfg.dates, 'DateRange') == 1

        [~, i_date_start] = min(abs(RCSXX.time - datetime(cfg.date_range(1),'TimeZone','America/Los_Angeles')));
       
        [~, i_date_end] = min(abs(RCSXX.time - datetime(cfg.date_range(2),'TimeZone','America/Los_Angeles')));

        [~, db_beh_date_start] = min(abs(...
             mean([db_beh_RCSXX.timeStart,db_beh_RCSXX.timeStop],2) ...
            - datetime(cfg.date_range(1),'TimeZone','America/Los_Angeles')));
       
        [~, db_beh_date_end] = min(abs(...
             mean([db_beh_RCSXX.timeStart,db_beh_RCSXX.timeStop],2) ...
            - datetime(cfg.date_range(2),'TimeZone','America/Los_Angeles')));
       
        
        date_range         = [RCSXX.time(i_date_start), RCSXX.time(i_date_end)]; 

        RCSXX              = RCSXX(i_date_start:i_date_end, :);


        db_beh_RCSXX       = db_beh_RCSXX(db_beh_date_start: db_beh_date_end, :);

    end
end