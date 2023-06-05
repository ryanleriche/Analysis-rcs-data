% Stage start dates, home, and clinic visits for RCS pts 1-7
stage_dates             = {{''}, {'08-Sep-2020'; '31-Jan-2021'; '31-May-2022'},... % RCS02
                           {''}, {'13-May-2021'; '12-Jul-2021'},... % RCS04
                           {'21-Jul-2021';'07-Sep-2021'},... % RCS05
                           {'07-Jun-2022','18-Aug-2022'},...
                           {'19-Oct-2022'}}; % RCS07


pt_meta.RCS02        = table;
pt_meta.RCS02.dates  = datetime(...
   {'27-Jul-2020'; '06-Aug-2020';...
    '08-Sep-2020'; '09-Sep-2020'; '10-Sep-2020'; '11-Sep-2020';...
    '13-Oct-2020';...
    ...
    '18-Oct-2020'; ...
    ...
    '02-Nov-2021'; '09-Nov-2021';...
    '01-Feb-2021'; '02-Feb-2021';  '03-Feb-2021';...
    '13-Apr-2021'; '14-Apr-2021';...
    '27-Sep-2021'; '28-Sep-2021';...
    ...
    '31-May-2022';...
    '27-Jun-2022';...
    '28-Jun-2022'; '29-Jun-2022';...
    '02-Jul-2022';...
    '21-Jul-2022'; '24-Jul-2022';...
    '06-Sep-2022'; '12-Sep-2022';...
    ...
    '14-Nov-2022'}, ...
    ...
    'TimeZone', 'America/Los_Angeles');

pt_meta.RCS02.desc   = ...
    {'s0_implant';     's0_explant';...
     's1_implant';     's1_inpatient';    's1_inclinic_d1';  's1_inclinic_d2';...
     'remove_hardware_inpatient';...
     ...
     's2_start';
     ...
     'washout_testing'; 'washout_testing';
     's2_inclinic_d1'; 's2_inclinic_d2' ; 's2_inclinic_d3';...
     's2_home_d1'    ; 's2_home_d2';...
     's2_inclinic_d1'; 's2_inclinic_d2' ;...
     ...
     's3_start_ol_versus_sham';...
     's3_break_start';...
     's3_inclinic_d1'; 's3_inclinic_d2';...
     's3_break_stop';...
     's3_break_start'; 's3_break_stop';...
     's3_break_start'; 's3_break_stop';...
     ...
     's3_stop_ol_versus_sham'};

%% RCS04
pt_meta.RCS04        = table;
pt_meta.RCS04.dates  = datetime(...
   {'09-Mar-2021'; '19-Mar-2021';...
    '13-May-2021'; '14-May-2021'; '15-May-2021'; '16-May-2021';...   
    '12-Jul-2021'; '13-Jul-2021'; '14-Jul-2021';...
    '17-Aug-2021'; '18-Aug-2021'; '19-Aug-2021';...
    '13-Jul-2022'; ...
    '12-Apr-2023';  '13-Apr-2023'},...
    ...
    'TimeZone', 'America/Los_Angeles');

pt_meta.RCS04.desc   = ...
    {'s0_implant';      's0_explant';...
     's1_implant';      's1_inpatient_d1';  's1_inpatient_d2'; 's1_inpatient_d3'; ...
     's2_inclinic_d1';  's2_inclinic_d2' ;  's2_inclinic_d3';...
     's2_home_d1'    ;  's2_home_d2';       's2_home_d2';...
     's2_home_d1'    ;...
     's2_home_d1'    ;  's2_home_d2';...
     };
%% RCS05

pt_meta.RCS05        = table;
pt_meta.RCS05.dates  = [datetime(...
   {'08-Jun-2021';   '18-Jun-2021'; ...
    '21-Jul-2021';   '22-Jul-2021';    '23-Jul-2021';...
    '07-Sep-2021';   '08-Sep-2021';    '09-Sep-2021';...
    '01-Dec-2021';   '02-Dec-2021';...
    '10-Aug-2022';   '11-Aug-2022';...
    '03-Mar-2023';    },...
    ...
    'TimeZone', 'America/Los_Angeles');...
    ...
    ...
    datetime('today', 'TimeZone', 'America/Los_Angeles')];

pt_meta.RCS05.desc   = ...
    {'s0_implant';       's0_explant';...
     's1_implant';       's1_inpatient_d1';   's1_inpatient_d2';...
     's2_inclinic_d1';   's2_inclinic_d2' ;   's2_inclinic_d3';...
     's2_home_d1'    ;   's2_home_d2';...
     's2_inclinic_d1';          's2_inclinic_d2';...
     'blinded_testing_start';   'blinded_testing_stop'};

%% RCS06
pt_meta.RCS06        = table;
pt_meta.RCS06.dates  = datetime(...
   {'01-Mar-2022';  '11-Mar-2022';...
    '07-Jun-2022';  '08-Jun-2022'; '09-Jun-2022';...
    '17-Aug-2022';  '18-Aug-2022'}, ...
    ...
    'TimeZone', 'America/Los_Angeles');

pt_meta.RCS06.desc   = ...
    {'s0_implant';       's0_explant';...
     's1_implant';       's1_inpatient_d1';  's1_inpatient_d2';...
     's2_inclinic_d1';   's2_inclinic_d2'};

%% RCS07
pt_meta.RCS07        = table;
pt_meta.RCS07.dates  = datetime(...
   {'20-Sep-2022'; '30-Sep-2022'; ...
    '19-Oct-2022'; '20-Oct-2022'; '21-Oct-2022';...
    '13-Dec-2022'; '14-Dec-2022'},...
    'TimeZone', 'America/Los_Angeles');

pt_meta.RCS07.desc   = ...
    {'s0_implant';       's0_explant';...
    's1_implant';        's1_inpatient_d1';  's1_inpatient_d2';...
    's2_inclinic_d1';   's2_inclinic_d2'};