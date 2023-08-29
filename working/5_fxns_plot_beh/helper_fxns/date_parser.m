function [redcap, date_range] = date_parser(cfg, redcap)
 
    
    if strcmp(cfg.dates, 'AllTime') == 1

        date_range           = [redcap.time(1), redcap.time(end)];

    % takes current time and finds the nearest report 
    elseif strcmp(cfg.dates, 'PreviousDays') == 1
        
        [~, i_date_start] = min(abs(redcap.time ...
            - (datetime('now','TimeZone','America/Los_Angeles') - days(cfg.ndays))));

        
        date_range         = [redcap.time(i_date_start), redcap.time(end)];

        if length(unique(date_range)) == 1
            i_date_start = i_date_start - 1;
        end

        date_range         = [redcap.time(i_date_start), redcap.time(end)];

        redcap              = redcap(i_date_start:end, :);

     % takes 'cfg.date_range' and finds the nearest report 
     elseif strcmp(cfg.dates, 'DateRange') == 1

         i_btwn = find(...
                     isbetween(...
                     redcap.time,  ...
                     datetime(cfg.date_range(1),'TimeZone','America/Los_Angeles'),...
                     datetime(cfg.date_range(2),'TimeZone','America/Los_Angeles')))
        
        date_range         = [redcap.time(i_btwn(1)), redcap.time(i_btwn(end))]; 
        
    end
end