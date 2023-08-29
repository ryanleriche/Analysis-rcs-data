
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


        file_fooof    = path(contains(path, 'fooof_params'));

        for i_fooof = 1 : length(file_fooof)
            tmp_params  = readtable(file_fooof{i_fooof});
    
            % format as cell containg double rather than character array
            var_class    = varfun(@class,tmp_params,'OutputFormat','cell');
            for i_var = 1 : length(var_class)
    
                if strcmp(var_class(i_var), 'cell')
    
                    tmp_params{:,i_var} =cellfun(@(x) ...
                            sscanf(x(2:end-2),'%f,'), tmp_params{:,i_var}, ...
                            'UniformOutput', false);
     
                end

                if contains(file_fooof{i_fooof}, 'MEAN')
                    mean_fooof_params = tmp_params;
                else
                    fooof_params = tmp_params;
                end
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


% 
%        var_names = {'path', 'fooof_params',...
%                      'pwr_spectra_aperiodic_rmv', 'aperiodic_fit'};
% 
% 
%         sense_chan(hemi_chan,var_names)    ...
%             ...
%             = [{path}, {fooof_params},...
%                {pwr_spectra_aperiodic_rmv}, {aperiodic_fit}];


        var_names = {'path', 'fooof_params',...
                     'pwr_spectra_aperiodic_rmv', 'aperiodic_fit', 'mean_fooof_params'};


        sense_chan(hemi_chan,var_names)    ...
            ...
            = [{path}, {fooof_params},...
               {pwr_spectra_aperiodic_rmv}, {aperiodic_fit}, {mean_fooof_params}];

        

        %%

    end
   

end

sense_chan = sortrows(sense_chan, 'RowNames');


%%
row_names = sense_chan.Properties.RowNames;

N_sess    = cellfun(@(x) ...
                   height(x), sense_chan.pwr_spectra_aperiodic_rmv, ...
                  'UniformOutput',false);
N_sess    = [N_sess{:}];

ALL_bands = table;

pb_range.theta_freq     = [3,7];
pb_range.alpha_freq     = [8,12];
pb_range.b_low_freq  = [13, 29];
pb_range.b_high_freq = [30, 40];


i_pt_hemi = [1, cumsum(N_sess)];

cnt = 1;
for j = 1 : height(sense_chan)

    if j == 1
        ALL_bands.theta_freq   = nan(sum(N_sess),1);
        ALL_bands.theta_pwr    = nan(sum(N_sess),1);
        ALL_bands.theta_width  = nan(sum(N_sess),1);


        ALL_bands.alpha_freq   = nan(sum(N_sess),1);
        ALL_bands.alpha_pwr    = nan(sum(N_sess),1);
        ALL_bands.alpha_width  = nan(sum(N_sess),1);

        ALL_bands.b_low_freq  = nan(sum(N_sess),1);
        ALL_bands.b_low_pwr  = nan(sum(N_sess),1);
        ALL_bands.b_low_width  = nan(sum(N_sess),1);

        ALL_bands.b_high_freq = nan(sum(N_sess),1);
        ALL_bands.b_high_pwr = nan(sum(N_sess),1);
        ALL_bands.b_high_width = nan(sum(N_sess),1);

    end

     ALL_bands.pt_chan(i_pt_hemi(j):i_pt_hemi(j+1)) = row_names(j);

     ALL_bands.pt_hemi(i_pt_hemi(j):i_pt_hemi(j+1)) = {row_names{j}(1:6)};


    for h = 1 : height(sense_chan.fooof_params{j})

        tmp_freq   = sense_chan.fooof_params{j}.center_freq{h};
        
        i_theta    = find(tmp_freq >= pb_range.theta_freq(1) & tmp_freq  <= pb_range.theta_freq(2));
        i_alpha    = find(tmp_freq >= pb_range.alpha_freq(1) & tmp_freq  <= pb_range.alpha_freq(2));
        i_b_low    = find(tmp_freq >= pb_range.b_low_freq(1) & tmp_freq  <= pb_range.b_low_freq(2));
        i_b_high   = find(tmp_freq >= pb_range.b_high_freq(1) & tmp_freq  <= pb_range.b_high_freq(2));

        tmp_pwr    = sense_chan.fooof_params{j}.peak_power{h};
        tmp_width  = sense_chan.fooof_params{j}.band_width{h};

        if ~isempty(i_theta)
            ALL_bands.theta_freq(cnt)      = tmp_freq(i_theta(1));
            ALL_bands.theta_pwr(cnt)       = tmp_pwr(i_theta(1));
            ALL_bands.theta_width(cnt)     = tmp_width(i_theta(1));
        end

        if ~isempty(i_alpha)
            ALL_bands.alpha_freq(cnt)      = tmp_freq(i_alpha(1));
            ALL_bands.alpha_pwr(cnt)       = tmp_pwr(i_alpha(1));
            ALL_bands.alpha_width(cnt)     = tmp_width(i_alpha(1));
        end

        if ~isempty(i_b_low)
            ALL_bands.b_low_freq(cnt)      = tmp_freq(i_b_low(1));
            ALL_bands.b_low_pwr(cnt)       = tmp_pwr(i_b_low(1));
            ALL_bands.b_low_width(cnt)     = tmp_width(i_b_low(1));
        end
    
        if ~isempty(i_b_high)
            ALL_bands.b_high_freq(cnt)     = tmp_freq(i_b_high(1));
            ALL_bands.b_high_pwr(cnt)      = tmp_pwr(i_b_high(1));
            ALL_bands.b_high_width(cnt)    = tmp_width(i_b_high(1));
        end
        cnt = cnt+1;
    end


end

%%
writetable(ALL_bands, fullfile(pp_dir, 'ALL_bands.xlsx'))


