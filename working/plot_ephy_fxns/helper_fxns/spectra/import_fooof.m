
function [sense_chan, fft_bins_inHz] = import_fooof(pp_dir, pt_ids, pp_subdir)


sense_chan = table;
freq_loaded = false;


for i = 1 : length(pt_ids)

    files_tbl = struct2table(...
                        dir(...
                        fullfile(pp_dir, pt_ids{i}, pp_subdir, '/exported_fooof_data/*.xlsx')...
                        ));

    files_tbl.path = cellfun(@(x, y) [x,'/',y], ...
                                files_tbl.folder, files_tbl.name,...
                                'UniformOutput', false);

    chan_name = unique(cellfun(@(x) ...
                            extractBefore(x , '_'), files_tbl.name,...
                            'UniformOutput', false));

    i_fooof = contains(chan_name, 'fooof');
    if any(i_fooof) && ~freq_loaded

         file_freq     = files_tbl.path{contains(files_tbl.path, 'fooof_freqs')};
         fft_bins_inHz = readmatrix(file_freq, 'Range', 2);                    
         freq_loaded   = true;

    end
    chan_name = chan_name(~i_fooof);



    for i_chan = 1 : length(chan_name)


        path     = files_tbl.path(...
                                contains(files_tbl.path, chan_name(i_chan)) ...
                                );

        file_aper     = path{contains(path, 'aperiodic_fit')};
        aperiodic_fit = readmatrix(file_aper, 'Range', 2);


        file_fooof    = path{contains(path, 'fooof_params')};
        fooof_params  = readtable(file_fooof);

        % format as cell containg double rather than character array
        var_class    = varfun(@class,fooof_params,'OutputFormat','cell');
        for i_var = 1 : length(var_class)

            if strcmp(var_class(i_var), 'cell')

                fooof_params{:,i_var} =cellfun(@(x) ...
                        sscanf(x(2:end-2),'%f,'), fooof_params{:,i_var}, ...
                        'UniformOutput', false);
             

            end
        end
        

        file_pwr                  = path{contains(path, 'pwr_spectra_aperiodic_rmv')};
        pwr_spectra_aperiodic_rmv = readmatrix(file_pwr, 'Range', 2);

        %%% now save hemisphere channel as row

        % since NK spectra are not saved out per hemisphere--add tag for
        % ease of comparision between the RCS and NK
        if ~endsWith(pt_ids{i}, {'R', 'L'})
            hemi = [pt_ids{i}, chan_name{i_chan}(1)];
        else
            hemi = pt_ids{i};
        end

        hemi_chan = [hemi,'_', chan_name{i_chan}];


        var_names = {'path', 'fooof_params',...
                     'pwr_spectra_aperiodic_rmv', 'aperiodic_fit'};


        sense_chan(hemi_chan,var_names)    ...
            ...
            = [{path}, {fooof_params},...
               {pwr_spectra_aperiodic_rmv}, {aperiodic_fit}];

    end
end