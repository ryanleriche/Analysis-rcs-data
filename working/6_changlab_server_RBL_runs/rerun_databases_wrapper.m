% call CONFIG_XXX to load directories, and REDcap data 
CONFIG_changlab_server;

%% import RCS databases, and INS logs per pt side as structures
%{
*   saves RCS session summaries as databases (db) from constellation of
    .jsons saved during streaming sessions

*   saves INS logs (with AppLog.txt, and EventLog.txt changes) as INS_logs
    importantly contains Adaptive state changes/stim defintions and Group changes

*   from db created a parsed database (par_db) allowing for .xlsx
    exportable and human readable one line summaries per streaming session

    -> Note LAST stim/sense settings are returned for parsimonious report
       of aDBS settings between streaming sessions

*   from db create a stimLog containing every change in stim parameter
    during a streaming session ("misses" offline PTM intiated changed)

Inital run takes hours for running multiple pts w/ 1000s of streaming
sessions.
%}

% option to load previous database for efficient processing
sub_cfg.ignoreold_db                = true;
sub_cfg.ignoreold_INS_logs          = true;
sub_cfg.ignoreold_par_db            = true;    % <-- keep as true to avoid version issues

% specify patient hemispheres
pt_sides        = {'RCS02R','RCS04R','RCS04L','RCS05R', 'RCS05L','RCS06R','RCS06L','RCS07L', 'RCS07R'};
%pt_sides         = {'RCS04R'};

for i = 1  : length(pt_sides)
    %%% process RCS .jsons into searchable database
    [db.(pt_sides{i}), bs.(pt_sides{i})] ...
        ...
        = makeDatabaseRCS_Ryan(...
        ...
        sub_cfg, pt_sides{i});

        %%% unpack all sense, LD, and stimulation settings as own variable in table
        % allows for programmatic discernment of unique RC+S settings
    [par_db.(pt_sides{i}), ~]...
        ...
        = makeParsedDatabaseRCS(...
        ...
        sub_cfg , pt_sides{i}, db);

end
% %%
% for i = 1  : length(pt_sides)
%     %%% process INS logs .txts based on unique entries only
%         % (INS logs have mostly repeating entries)
%     INS_logs.(pt_sides{i})  ...
%         ...
%         = RCS_logs( ...
%         ...
%          sub_cfg, pt_sides{i});
% end