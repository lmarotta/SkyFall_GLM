%% ParseProbeOutput
% Inputs:
% optional
% YYYY_start, MM_start, DD_start, HH_start, MIN_start - data cutoff start
% in UTC timestamp (+5h from Chicago time)
% YYYY_end, MM_end, DD_end, HH_end, MIN_end  - data cutoff end in UTC
% timestamp (+5h from Chicago time)
%
% Loads a user-selected .txt file containing exported data from Fallnet
% and current model class_params file
%
% saves labels structure with parsed data into FallProbe_TestData.mat file

function ParseProbeOutput(YYYY_start, MM_start, DD_start, HH_start, MIN_start, YYYY_end, MM_end, DD_end, HH_end, MIN_end)

if nargin ~= 0 &&  nargin ~= 10
  error('the function might not take any input or you may specify cutoff start and end')
end

% load file
[FileName,PathName,~] = uigetfile('*.txt');
PayloadTable=readtable([PathName '/' FileName],'Delimiter','\t');
Phone = PayloadTable.UserID;
Probe = PayloadTable.Probe;
Payload=PayloadTable.Payload;

%load current fall model file from github
[FileName,PathName,~] = uigetfile('*.mat');
model=load([PathName '/' FileName]);
Nfeatures = length(model.b);

if nargin == 10
    cutoff_start = datetime(YYYY_start,MM_start,DD_start,HH_start,MIN_start,0);  %This is in UTC (+5h from Chicago time)
    cutoff_end = datetime(YYYY_end,MM_end,DD_end,HH_end,MIN_end,0);  %This is in UTC (+5h from Chicago time)
    
    timestr=cellfun(@(x) strsplit(x(strfind(x,'"TIMESTAMP"'):strfind(x,'"TIMESTAMP"')+100),{':' ','}),Payload,'UniformOutput',false);
    timestamp=cellfun(@(x) str2double(x{:,2}),timestr);
    timestamp=datetime(1970,1,1,0,0,timestamp);

    filtered_times=timestamp>cutoff_start & timestamp<cutoff_end;
    Phone = Phone(filtered_times);
    Probe = Probe(filtered_times);
    Payload = Payload(filtered_times);
end

%indices of correct and failure records
correct = cellfun(@(x) isempty(strfind(x,'Failure')) && ~isempty(strfind(x,'FallNet')),Probe);
ind_correct = find(correct);
failure = cellfun(@(x) ~isempty(strfind(x,'Failure')),Probe);
ind_failures = find(failure);

%Indices of Activity Detection records
act_inds = find(cellfun(@(x) ~isempty(strfind(x,'ActivityDetection')),Probe));

%Initialize structure contents (each element corresponds to a window)
winsize = -ones(length(ind_correct),1);
timestampSTART_END=-ones(length(ind_correct),2); %start end of each window 
evalstart = -ones(length(ind_correct),1);  %timestamp when the model started to evaluate  
sensor_counts =-ones(length(ind_correct),3); % # of samples per sensor
duration = -ones(length(ind_correct),3);  %evaluation time for each record (preparation, verification, evaluation)  
acce = cell(length(ind_correct),1);      %sensor data for each window
gyro = cell(length(ind_correct),1);
baro = cell(length(ind_correct),1);
%FAILURES
failure_timestampSTART_END = -ones(length(ind_failures),2); %start and end timestamps for sensor data in the window
failure_reason = cell(length(ind_failures),1);
failure_value = -ones(length(ind_failures),1);
failure_duration = -ones(length(ind_failures),3); 
failure_evalstart = -ones(length(ind_failures),1); %timestamp when the model started to evaluate 

for ind=1:length(ind_correct)
    
    i_correct = ind_correct(ind); %index
    temp_correct=strsplit(Payload{i_correct},{'{' '"' ':' ',' ' ' '[' ']' '}'}); 
 
    x=find(strcmp(temp_correct,'EVALUATION_WINDOW_SIZE'));
    winsize(ind) = str2double(temp_correct(x+1));
    
    %Additional Start and End timestamps for the clip
    xs = find(strcmp(temp_correct,'MIN_ABSOLUTE_TIMESTAMP'));
    xe = find(strcmp(temp_correct,'MAX_ABSOLUTE_TIMESTAMP'));
    timestampSTART_END(ind,:) = [str2double(temp_correct(xs+1)) str2double(temp_correct(xe+1))];
        
    eval_start_ts = find(strcmp(temp_correct,'EVALUATION_START'));
    evalstart(ind) = str2double(temp_correct(eval_start_ts+1));
    
    ac = find(strcmp(temp_correct,'ACCELEROMETER_READING_COUNT'));
    gc = find(strcmp(temp_correct,'GYROSCOPE_READING_COUNT'));
    bc = find(strcmp(temp_correct,'BAROMETER_READING_COUNT'));
    sensor_counts(ind,:)=[str2double(temp_correct(ac+1)) str2double(temp_correct(gc+1)) str2double(temp_correct(bc+1))];
    
    % parse duration times
    prepd = find(strcmp(temp_correct,'PREPARE_DURATION'));
    verifd = find(strcmp(temp_correct,'VERIFY_DURATION')); 
    predd = find(strcmp(temp_correct,'PREDICT_DURATION'));
    duration(ind,:)=[str2double(temp_correct(prepd+1)) str2double(temp_correct(verifd+1)) str2double(temp_correct(predd+1))];
    
    % parse raw sensors data
    accstart = find(strcmp(temp_correct,'ACCELEROMETER_SAMPLES'));
    gyrstart = find(strcmp(temp_correct,'GYROSCOPE_SAMPLES'));
    barstart = find(strcmp(temp_correct,'BAROMETER_SAMPLES'));
    sensortimestamps = find(strcmp(temp_correct,'TIMESTAMPS'));
    ix = find(strcmp(temp_correct,'X'));
    iy = find(strcmp(temp_correct,'Y'));
    iz = find(strcmp(temp_correct,'Z'));
    
    ixacc = ix(find(ix > accstart,1));
    iyacc = iy(find(iy > accstart,1));
    izacc = iz(find(iz > accstart,1));
    stsacc = sensortimestamps(find(sensortimestamps > accstart,1));
    Nsamplesacc = sort([ixacc iyacc izacc stsacc]);
    Nsamplesacc = Nsamplesacc(2)-Nsamplesacc(1)-1;
    acce{ind} = [str2double(temp_correct(stsacc+1:stsacc+Nsamplesacc))' str2double(temp_correct(ixacc+1:ixacc+Nsamplesacc))' ...
        str2double(temp_correct(iyacc+1:iyacc+Nsamplesacc))' str2double(temp_correct(izacc+1:izacc+Nsamplesacc))'];
    
    ixgyr = ix(find(ix > gyrstart,1));
    iygyr = iy(find(iy > gyrstart,1));
    izgyr = iz(find(iz > gyrstart,1));
    stsgyr = sensortimestamps(find(sensortimestamps > gyrstart,1));
    Nsamplesgyr = sort([ixgyr iygyr izgyr stsgyr]);
    Nsamplesgyr = Nsamplesgyr(2)-Nsamplesgyr(1)-1;
    gyro{ind} = [str2double(temp_correct(stsgyr+1:stsgyr+Nsamplesgyr))' str2double(temp_correct(ixgyr+1:ixgyr+Nsamplesgyr))' ...
        str2double(temp_correct(iygyr+1:iygyr+Nsamplesgyr))' str2double(temp_correct(izgyr+1:izgyr+Nsamplesgyr))'];
    
    ipbar = find(strcmp(temp_correct,'PRESSURE'));
    iabar = find(strcmp(temp_correct,'ALTITUDE'));
    stsbar = sensortimestamps(find(sensortimestamps > barstart,1));
    Nsamplesbar = sort([ipbar iabar stsbar]);
    Nsamplesbar = Nsamplesbar(2)-Nsamplesbar(1)-1;
    baro{ind} = [str2double(temp_correct(stsbar+1:stsbar+Nsamplesbar))' str2double(temp_correct(ipbar+1:ipbar+Nsamplesbar))' ...
        str2double(temp_correct(iabar+1:iabar+Nsamplesbar))'];
        
end

labels.winsize = winsize;
labels.timestampSTART_END = timestampSTART_END; %fallnet model features
labels.evalstart = evalstart;  %timestamp when the model started to evaluate  
labels.sensor_counts = sensor_counts; % # of samples per sensor
labels.duration = duration;  %evaluation time for each record (preparation, verification, evaluation)  
labels.acce = acce;      %sensor data for each window
labels.gyro = gyro;
labels.baro = baro;

labels.phone = Phone(ind_correct);

for ind=1:length(ind_failures)
    i = ind_failures(ind); %index
    temp=strsplit(Payload{i},{'{' '"' ':' ',' ' ' '[' ']' '}'});
    
    %Additional Start and End timestamps for the clip
    xs = find(strcmp(temp,'MIN_ABSOLUTE_TIMESTAMP'));
    xe = find(strcmp(temp,'MAX_ABSOLUTE_TIMESTAMP'));
    failure_timestampSTART_END(ind,:) = [str2double(temp(xs+1)) str2double(temp(xe+1))];
        
    eval_start_ts = find(strcmp(temp,'EVALUATION_START'));
    failure_evalstart(ind) = str2double(temp(eval_start_ts+1));
    
    % parse duration times
    prepd = find(strcmp(temp,'PREPARE_DURATION'));
    verifd = find(strcmp(temp,'VERIFY_DURATION')); 
    predd = find(strcmp(temp,'PREDICT_DURATION'));
    failure_duration(ind,:)=[str2double(temp(prepd+1)) str2double(temp(verifd+1)) str2double(temp(predd+1))];
    
    failind = find(strcmp(temp,'FAILURE_REASON'));
    failure_reason(ind) = temp(failind+1);
    
    if ~isempty(strfind(temp{failind+1},'INSUFFICIENT'))
        fvind = find(strcmp(temp,'SAMPLES_COUNT'));
        failure_value(ind) = str2double(temp(fvind+1));
    elseif ~isempty(strfind(temp{failind+1},'GAP'))
        fvind = find(strcmp(temp,'GAP_SIZE'));
        failure_value(ind) = str2double(temp(fvind+1));    
    end
end

labels.failure.timestampSTART_END = failure_timestampSTART_END; %start and end timestamps for sensor data in the window
labels.failure.reason = failure_reason;
labels.failure.value = failure_value;
labels.failure.duration = failure_duration; 
labels.failure.evalstart = failure_evalstart; %timestamp when the model started to evaluate 
labels.failure.phone = Phone(ind_failures);

act_type_all={};
act_conf_all={};
timestamps={};

for i=1:length(act_inds)
    ind = act_inds(i);
    temp = strsplit(Payload{ind},{'{' '"' ':' ',' ' ' '[' ']' '}'});
    
    act_counts = temp{find(strcmp(temp,'ACTIVITY_COUNT'))+1};
    
    act_type={};
    act_conf=[];
    
    type_inds = find(strcmp(temp,'ACTIVITY_TYPE'));
    conf_inds = find(strcmp(temp,'ACTIVITY_CONFIDENCE'));
    
    for j=1:str2double(act_counts)
        act_type = [act_type temp{type_inds(j)+1}];
        act_conf = [act_conf str2double(temp{conf_inds(j)+1})]; 
    end
    
    act_type_all = [act_type_all {act_type}];
    act_conf_all = [act_conf_all act_conf];
    timestamps = [timestamps temp(find(strcmp(temp,'TIMESTAMP'))+1)];
end

activity_detection = struct('Type',act_type_all,'Conf',act_conf_all,'Timestamp',timestamps);

save FallProbe_TestData labels activity_detection -v7.3
