function s0_REDcap    = append_pain_BD(s0_map_dir, s0_REDcap, visits_tbl)
% 
% visits_tbl = visits;
% 
% s0_map_dir  = [dropbox_dir, ...
%             'Stage 0 ALL PATIENTS Redcap Records/PBD_metrics/'];


s0_maps_tbl = readtable([s0_map_dir, 'ALL pain metrics.csv']);
s0_maps_tbl = removevars(s0_maps_tbl, {'Index', 'hoursum', 'hourdiff', 'nrsdiff'...
                'intensitydiff','mood_vasdiff', 'unpleasant_vasdiff'});

pt_ids      = fieldnames(s0_REDcap);
pt_ids      = cellfun(@(x) x(7:end), pt_ids, 'UniformOutput', false);


%find(s0_maps_tbl.avgpainintensity == min(s0_maps_tbl.avgpainintensity))

for i = 1 : length(pt_ids)

    redcap  = s0_REDcap.(['stage0', pt_ids{i}]);
    tmp_tbl = s0_maps_tbl(strcmp(s0_maps_tbl.ID, pt_ids{i}), :);

    
    % only pain maps that occured DURING pt's Stage 0 (rather than DEMO
    % day)
    i_s0_start = contains(visits_tbl.(pt_ids{i}).desc, 's0_implant');
    i_s0_stop  = contains(visits_tbl.(pt_ids{i}).desc, 's0_explant');
    
    tmp_tbl.arm1_timestamp.TimeZone = 'America/Los_Angeles';


    if any(le(tmp_tbl.arm1_timestamp, datetime('01-Jan-2019','TimeZone','America/Los_Angeles')))
    
        % if every but NOT all datetimes are shifted -> toss error
        if ~all(le(tmp_tbl.arm1_timestamp, datetime('01-Jan-2019','TimeZone','America/Los_Angeles')))
        
            error('RBL: inconsistent datetime formatting w/n %s .csv file', pt_id)
        end

        tmp_tbl.arm1_timestamp = tmp_tbl.arm1_timestamp + calyears(2000);


    end
    
    



    % from the entirety of implant day to the entiretyof the explant day
    i_s0       = ge(tmp_tbl.arm1_timestamp, visits_tbl.(pt_ids{i}).dates(i_s0_start)) &...
                 le(tmp_tbl.arm1_timestamp, visits_tbl.(pt_ids{i}).dates(i_s0_stop) + duration('24:00:00'));
    
    tmp_tbl   = tmp_tbl(i_s0, :);



    % define pain body diagram (PBD) metrics
    redcap.PBD_sum    = nan(height(redcap), 1);
    redcap.PBD_coverage    = nan(height(redcap), 1);
    redcap.PBD_mean    = nan(height(redcap), 1);

    % add timezone and round
    tmp_time = redcap.scales_vasmpq_timestamp;
    tmp_time = dateshift(tmp_time, 'start', 'minute');

    for j = 1 : height(tmp_tbl)

        i_map = find(eq(tmp_time, tmp_tbl.arm1_timestamp(j)));

        redcap.PBD_sum(i_map)       = tmp_tbl.pixelsum(j);
        redcap.PBD_coverage(i_map)       = tmp_tbl.painspread(j);
        redcap.PBD_mean(i_map)       = tmp_tbl.avgpainintensity(j);
    
    end

    if length(find(~isnan(redcap.PBD_sum))) == ...
            length(find(~isnan(tmp_tbl.pixelsum)))

        disp(['stage0', pt_ids{i}, ' | pain maps aligned to Arm 1 timestamp'])

    else 
        error(['stage0', pt_ids{i}, ' | see for timing disprepancy btwn',...
            ' pain maps and pain surveys'])
    end

    % percent of pixels w/ any hue
    redcap.PBD_coverage       = redcap.PBD_coverage * 100;


    % based on sex of pts, male or female body maps have different N pixels 
    % (and shape which is NOT explicitly considered here)

%     if any(strcmp(pt_ids{i}, {'RCS02', 'RCS04', 'RCS07'}))
% 
% 
%         redcap.PBD_sum   = redcap.PBD_sum / (820452);
% 
%     elseif any(strcmp(pt_ids{i}, {'RCS05', 'RCS06'}))
% 
%         redcap.PBD_sum   = redcap.PBD_sum / (724608);
% 
% 
%     end
%     

    s0_REDcap.(['stage0',pt_ids{i}]) = redcap;
end
end