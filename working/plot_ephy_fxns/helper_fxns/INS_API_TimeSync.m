%% by-session INS time to API time latency summary
function INS_API_latency  = INS_API_TimeSync(TimeSync_filename)

% TimeSync_filename = ['/Users/Leriche/pia_server/datastore_spirit/human/rcs_chronic_pain/rcs_device_data/',...
%             'raw/RCS02/SummitData/SummitContinuousBilateralStreaming/RCS02R/Session1673885676526/',...
%             'DeviceNPC700457H/TimeSync.json'];



% load .json as struct -> convert into table for clarity
TimeSync_sess = deserializeJSON(TimeSync_filename);

TimeSyncDataHeader     = struct2table([TimeSync_sess.TimeSyncData.Header]);
TimeSyncData           = removevars(struct2table([TimeSync_sess.TimeSyncData]), 'Header');


TimeSyncData           = [TimeSyncDataHeader, TimeSyncData];
TimeSyncData.timestamp = [TimeSyncData.timestamp.seconds]';

% 'timestamp':
% "Elapsed number of seconds since March 1, 2000, in units of seconds.
% Implemented in INS firmware" (Sellers et al 2021 Front. Hum. Neurosci.)

time_INS               = datetime(2000, 3, 1, 0, 0, TimeSyncData.timestamp);
time_INS.TimeZone      = 'America/Los_Angeles';

% 'PacketGenTime':
% "API estimate of when the packet was created on the INS. Unix time with
% resolution to millisecond" (Sellers et al 2021)

time_API = datetime(TimeSyncData.PacketGenTime / 1000,...
                    'ConvertFrom','posixTime','TimeZone','America/Los_Angeles',...
                    'Format','dd-MMM-yyyy HH:mm:ss.SSS');

t_INS_sub_API     = time_INS - time_API;


% simple summary statistics should be sufficient
INS_API_latency.mean = mean(t_INS_sub_API);
INS_API_latency.std  = std(t_INS_sub_API);

INS_API_latency.max  = max(t_INS_sub_API);
INS_API_latency.min  = min(t_INS_sub_API);


% for trouble-shooting uncomment below and see each TimeSync.json entry
% INS_API_Sync_tbl = table(time_INS, time_API,t_INS_sub_API);

end