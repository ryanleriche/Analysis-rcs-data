function proc_g_changes = determine_current_group_status(proc_g_changes)


proc_g_changes(...
    contains(proc_g_changes.event, 'LeadIntegrityTest'), :) = [];


i_grp    = contains(proc_g_changes.event, 'Group');
i_status = find(~i_grp);
i_grp    = find(i_grp);

proc_g_changes.Group_Status = cell(height(proc_g_changes),1);


for i_event = 1 : height(proc_g_changes)
    if contains(proc_g_changes.event(i_event), {'On', 'Off'})

        ind_diff = i_grp - i_event;
        N_diff   = max(ind_diff(ind_diff < 0));

        if isempty(N_diff) % inital INS Logs prior to first Group defintion
            
            proc_g_changes.Group_Status{i_event} =proc_g_changes.event{i_event};

        else               % append previous group to current therapy status 
            i_near_prev = i_grp(ind_diff == N_diff);

            proc_g_changes.Group_Status{i_event} ...
                ...
                =...
            [proc_g_changes.event{i_near_prev}, '_',proc_g_changes.event{i_event}];
        end

    else % append previous **therapy status** to current group

        ind_diff   = i_status - i_event;
        N_diff     = max(ind_diff(ind_diff < 0));

        i_near_prev = i_status(ind_diff == N_diff);

        proc_g_changes.Group_Status{i_event} ...
            ...
            =...
        [proc_g_changes.event{i_event}, '_', proc_g_changes.event{i_near_prev}];
    end
end

% return therapyStatus boolean for ease of parsing therapy "Off"s
proc_g_changes.therapyStatus = contains(proc_g_changes.Group_Status, 'On');

% identical times were actually loaded in order, but missed due to INS 1 Hz timestamping
% -> prevent reordering by adding in half second in between

end
