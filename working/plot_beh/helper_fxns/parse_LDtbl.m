function LD_tbl = parse_LDtbl(input_LD_tbl)  

    LD_tbl         = struct2table(input_LD_tbl);
    tmp_ld_feat_tbl = table();
    
    for d=1:4
        LD_feat = struct2table(cellfun(@(x) x(d), LD_tbl.features));

        tmp_LD  = varfun(@(x) double(typecast(int64(x),'int32')), LD_feat);

       
        LD_feat = renamevars(tmp_LD(1:2:height(tmp_LD),:),...
                   {'Fun_normalizationMultiplyVector', 'Fun_normalizationSubtractVector','Fun_weightVector'},...
                   {['normMultiply', num2str(d-1)],['normSubtract', num2str(d-1)],['normWeight', num2str(d-1)]});
    
    
        tmp_ld_feat_tbl = [tmp_ld_feat_tbl, LD_feat]; %#ok<AGROW> 
    
    end

    pb_input    = repmat({'Disabled'}, height(LD_tbl), 8);

    for d=1:8

        i_pb_input_d = find(cellfun(@(x) length(x), LD_tbl.powerband_Inputs) >=d);

        tmp_pb = cellfun(@(x) x(d), LD_tbl.powerband_Inputs(i_pb_input_d));
        pb_input(i_pb_input_d, d) = tmp_pb;

    end

     pb_input_tbl = cell2table(pb_input, "VariableNames", ...
        cellfun(@(x) ['PowerBand_input_',num2str(x)], num2cell(0:7), 'UniformOutput', false));

    
    tmp_b0 = cellfun(@(x) ...
                double(typecast(int64(x(1)),'int32')),...
                LD_tbl.biasTerm, 'UniformOutput',false);

    tmp_ld_feat_tbl.biasTerm0 = cellfun(@(x) x(1), tmp_b0);

    
    tmp_b1 = cellfun(@(x) ...
                double(typecast(int64(x(2)),'int32')),...
                LD_tbl.biasTerm, 'UniformOutput',false);

    tmp_ld_feat_tbl.biasTerm1 = cellfun(@(x) x(1), tmp_b1);
    


    LD_tbl = [tmp_ld_feat_tbl,  LD_tbl,pb_input_tbl];
    LD_tbl = removevars(LD_tbl, {'biasTerm', 'features','powerband_Inputs'});

end
    