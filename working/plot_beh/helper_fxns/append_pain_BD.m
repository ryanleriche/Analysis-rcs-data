function s0_REDcap    = append_pain_BD(s0_map_dir, s0_REDcap, visits_tbl)

% visits_tbl = visits;
% 
% s0_map_dir  = [dropbox_dir, ...
%             '/DATA ANALYSIS/Stage 0 ALL PATIENTS Redcap Records/pain_map_metrics/'];
% 

s0_maps_tbl = readtable([s0_map_dir, 'ALL pain metrics.csv']);
s0_maps_tbl = removevars(s0_maps_tbl, {'Index', 'hoursum', 'hourdiff', 'nrsdiff'...
                'intensitydiff','mood_vasdiff', 'unpleasant_vasdiff'});

pt_ids      = fieldnames(s0_REDcap);
pt_ids      = cellfun(@(x) x(7:end), pt_ids, 'UniformOutput', false);

for i = 1 : length(pt_ids)

    redcap  = s0_REDcap.(['stage0', pt_ids{i}]);
    tmp_tbl = s0_maps_tbl(strcmp(s0_maps_tbl.ID, pt_ids{i}), :);


    % only pain maps that occured DURING pt's Stage 0 (rather than DEMO
    % day)
    i_s0_start = contains(visits_tbl.(pt_ids{i}).desc, 's0_implant');
    i_s0_stop  = contains(visits_tbl.(pt_ids{i}).desc, 's0_explant');
    
    tmp_tbl.arm1_timestamp.TimeZone = 'America/Los_Angeles';

    % from the entirety of implant day to the entiretyof the explant day
    i_s0       = ge(tmp_tbl.arm1_timestamp, visits_tbl.(pt_ids{i}).dates(i_s0_start)) &...
                 le(tmp_tbl.arm1_timestamp, visits_tbl.(pt_ids{i}).dates(i_s0_stop) + duration('24:00:00'));
    
    tmp_tbl   = tmp_tbl(i_s0, :);



    % define pain body diagram (PBD) metrics
    redcap.PBD_sumpixelval    = nan(height(redcap), 1);
    redcap.PBD_pixelspread    = nan(height(redcap), 1);
    redcap.PBD_avgpixelval    = nan(height(redcap), 1);

    % add timezone and round
    tmp_time = redcap.scales_vasmpq_timestamp;
    tmp_time = dateshift(tmp_time, 'start', 'minute');

    for j = 1 : height(tmp_tbl)

        i_map = find(eq(tmp_time, tmp_tbl.arm1_timestamp(j)));

        redcap.PBD_sumpixelval(i_map)       = tmp_tbl.pixelvalue(j);
        redcap.PBD_pixelspread(i_map)       = tmp_tbl.painspread(j);
        redcap.PBD_avgpixelval(i_map)       = tmp_tbl.avgpainintensity(j);
    
    end

    if length(find(~isnan(redcap.PBD_sumpixelval))) == ...
            length(find(~isnan(tmp_tbl.pixelvalue)))

        disp(['stage0', pt_ids{i}, ' | pain maps aligned to Arm 1 timestamp'])

    else 
        error(['stage0', pt_ids{i}, ' | see for timing disprepancy btwn',...
            ' pain maps and pain surveys'])
    end

    % percent of pixels w/ any hue
    redcap.PBD_pixelspread       = redcap.PBD_pixelspread * 100;


    % based on sex of pts, male or female body maps have different N pixels 
    % (and shape which is NOT explicitly considered here)

    if any(strcmp(pt_ids{i}, {'RCS02', 'RCS04', 'RCS07'}))


        redcap.PBD_sumpixelval   = redcap.PBD_sumpixelval / (820452);

    elseif any(strcmp(pt_ids{i}, {'RCS05', 'RCS06'}))

        redcap.PBD_sumpixelval   = redcap.PBD_sumpixelval / (724608);


    end

    s0_REDcap.(['stage0',pt_ids{i}]) = redcap;
end
end