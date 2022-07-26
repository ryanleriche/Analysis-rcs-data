function db_RCSXXX_out = parse_db(db_RCSXXX)
    % db_RCSXXX_out = removevars(db_RCSXXX, {'contacts', 'amp', 'PW','freq','stimName'});
    

    for i_sess = 1 : height(db_RCSXXX)
        
    
        if height(db_RCSXXX.stimparams{i_sess}) <= 1
    
            i_para = strsplit(char(db_RCSXXX_out.stimparams{i_sess}),',');
        
            if ~isempty(i_para{1})
        
                db_RCSXXX_out.contacts{i_sess}    = i_para{1};
                db_RCSXXX_out.amp(i_sess)         = str2double(i_para{2}(2:end -2));
                db_RCSXXX_out.PW(i_sess)          = str2double(i_para{3}(2:end -2));
                db_RCSXXX_out.freq(i_sess)        = str2double(i_para{4}(2:end -2));
        
            else % without stim parameters 
                db_RCSXXX_out.contacts{i_sess}    = '';
                db_RCSXXX_out.amp(i_sess)         = NaN;
                db_RCSXXX_out.PW(i_sess)          = NaN;
                db_RCSXXX_out.freq(i_sess)        = NaN;
    
            end
    
        else
            db_RCSXXX_out.contacts{i_sess}    = 'MANY';
            db_RCSXXX_out.amp(i_sess)         = NaN;
            db_RCSXXX_out.PW(i_sess)          = NaN;
            db_RCSXXX_out.freq(i_sess)        = NaN;
        end

        if ~isempty(db_RCSXXX.stimName{i_sess})
        
            db_RCSXXX_out.stimName{i_sess}    = db_RCSXXX.stimName{i_sess};
        else
        
            db_RCSXXX_out.stimName{i_sess}    = '';
        end
        
    end
end