function PR_r_cap = import_presidio_csvs(raw_dir)


PR_r_cap           = struct;
csv_PR_tbl         = struct2table(dir(fullfile(raw_dir,'/*.csv')));

for i = 1: height(csv_PR_tbl)

    file_name    = csv_PR_tbl{i, 'name'}{1};

    path_oi      = fullfile(csv_PR_tbl{i, 'folder'}{1}, file_name);

    pr_beh   = readtable(path_oi);

    vars         = pr_beh.Properties.VariableNames;
    classes      = varfun(@class, pr_beh, 'OutputFormat','cell');
       
    i_double     = find(strcmp(classes, 'double'));
    i_rmv        = all(isnan(pr_beh{:,i_double}));

    all_nan_vars = vars(i_double(i_rmv));

    pr_beh   = removevars(pr_beh, all_nan_vars);


    classes      = varfun(@class, pr_beh, 'OutputFormat','cell');


    pr_beh = renamevars(pr_beh, vars(strcmp(classes, 'datetime')), 'time_hamd');

    
    PR_r_cap.(file_name(1:end-4)) = pr_beh;

end

end
