function db_RCSXXX = parse_db(db_RCSXXX)

    for i_sess = 1 : height(db_RCSXXX)
    
        if length(db_RCSXXX.stimparams{i_sess}) <= 1
    
            i_para = strsplit(char(db_RCSXXX.stimparams{i_sess}),',');
        
            if ~isempty(i_para{1})
        
                db_RCSXXX.contacts{i_sess}    = i_para{1};
                db_RCSXXX.amp(i_sess)         = str2double(i_para{2}(2:end -2));
                db_RCSXXX.PW(i_sess)          = str2double(i_para{3}(2:end -2));
                db_RCSXXX.freq(i_sess)        = str2double(i_para{4}(2:end -2));
        
            else % without stim parameters 
                db_RCSXXX.contacts{i_sess}    = '';
                db_RCSXXX.amp(i_sess)         = NaN;
                db_RCSXXX.PW(i_sess)          = NaN;
                db_RCSXXX.freq(i_sess)        = NaN;
    
            end
    
        else
            db_RCSXXX.contacts{i_sess}    = 'MANY';
            db_RCSXXX.amp(i_sess)         = NaN;
            db_RCSXXX.PW(i_sess)          = NaN;
            db_RCSXXX.freq(i_sess)        = NaN;
        end
    end
end