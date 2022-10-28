% ARRANGE THE REDCAP REPORS AND STIM GROUPS BASED ON DATA BASE ONLY 


try
    [stimLog_w_redcap.([cfg.pt_id 'L'])] = align_REDcap_to_stimLog(cfg, db.([cfg.pt_id 'L']), REDcap.(cfg.pt_id));
     fprintf(' - OK done %s \n',[cfg.pt_id 'L'])
catch 
    fprintf(' - NO stimlog for %s \n',[cfg.pt_id 'L'])

end 

try 
    [stimLog_w_redcap.([cfg.pt_id 'R'])] = align_REDcap_to_stimLog(cfg, db.([cfg.pt_id 'R']), REDcap.(cfg.pt_id));
    fprintf(' - OK done %s \n',[cfg.pt_id 'R'])
catch 
    fprintf(' - NO stimlog for %s \n',[cfg.pt_id 'R'])

end


% MAKE STIM GROUPS
fprintf('Making stim groups:\n')
if contains(cfg.pt_id,{'01','02'})  %if 01 or 02, make L side empty
        stimLogLEFT   = [];
else
        stimLogLEFT     = stimLog_w_redcap.([cfg.pt_id 'L']);
end

stimLogRIGHT     = stimLog_w_redcap.([cfg.pt_id 'R']);


[redcap, stimGroups,sortedmetrics] =  make_stim_groups(...
     cfg.pt_id, stimLogLEFT, stimLogRIGHT, REDcap.(cfg.pt_id), cfg.visits.(cfg.pt_id));

