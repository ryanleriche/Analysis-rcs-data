function ch_names = parse_chan_names(pt_side_id, par_db_oi_out)

switch pt_side_id
    case 'RCS02R';  stimReg    = [{'RACC ', ["0","1","2","3"]}; {'RThal ', ["8","9","10","11"]}];
             
    case 'RCS04L';  stimReg    = [{'LACC ', ["0","1","2","3"]}; {'LCaud ', ["8","9","10","11"]}];
    case 'RCS04R';  stimReg    = [{'RACC ', ["0","1","2","3"]}; {'RThal ', ["8","9","10","11"]}];
        
        
    case 'RCS05L';  stimReg    = [{'LCaud ', ["0","1","2","3"]}; {'LACC ', ["8","9","10","11"]}];
    case 'RCS05R';  stimReg    = [{'RThal ', ["0","1","2","3"]}; {'RIFG ', ["8","9","10","11"]}];

        
    case 'RCS06L';  stimReg    = [{'LACC ',  ["0","1","2","3"]}; {'LCaud ', ["8","9","10","11"]}];
    case 'RCS06R';  stimReg    = [{'RThal ', ["0","1","2","3"]}; {'RSFG ', ["8","9","10","11"]}];
     
    % SGC (not ACC) was bilaterally implanted (RCS data is wrong)
    case 'RCS07L';  stimReg    = [{'LGPi ',  ["0","1","2","3"]}; {'LSGC ', ["8","9","10","11"]}];
    case 'RCS07R';  stimReg    = [{'RThal ', ["0","1","2","3"]}; {'RSGC ', ["8","9","10","11"]}];    
end


ch_tbl_var   = compose("Ch%g_plusInput", 0:3);
ch_pos_input = par_db_oi_out{:,ch_tbl_var};


ch_tbl_var   = compose("Ch%g_chanFullStr", 0:3);
ch_names     = par_db_oi_out{:,ch_tbl_var};


for j_ch = 1  :  size(ch_names, 2)

    tf_0_3   = cellfun(@(x) any(strcmp(x, stimReg{1,2})), ch_pos_input(:,j_ch));

    tf_8_11  = cellfun(@(x) any(strcmp(x, stimReg{2,2})), ch_pos_input(:,j_ch));

    if all(tf_0_3) && ~all(tf_8_11)
        ch_names(:, j_ch)  = cellfun(@(x) [stimReg{1,1}, x], ch_names(:, j_ch), 'UniformOutput',false);
        
    elseif ~all(tf_0_3) && all(tf_8_11)
         ch_names(:, j_ch) = cellfun(@(x) [stimReg{2,1}, x], ch_names(:, j_ch), 'UniformOutput',false);
    end
end
end