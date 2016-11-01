% Loads a user-selected .txt file containing exported data from Fallnet
% and sensor probes and separates probe data fro analysis

clear all

% load file
[FileName,PathName,~] = uigetfile('*.txt');
PayloadTable=readtable([PathName '/' FileName],'Delimiter','\t');
Probe = PayloadTable.Probe;
Payload=PayloadTable.Payload;

%load current fall model file from github
[FileName,PathName,~] = uigetfile('class_params*.mat');
model=load([PathName '/' FileName]);
Nfeatures = length(model.b);

fall_probe=[];
acc_probe=[];
gyr_probe=[];
bar_probe=[];
fall_probe_failure_counts = [];
fall_probe_failure_gap = [];
fall_probe_failure_interp =[];
fall_probe_failure_other =[];

cutoff_start = datetime(2016,10,31,15,0,0,0);  %This is in UTC (+5h from Chicago time)
cutoff_end = datetime(2016,10,31,19,0,0,0);  %This is in UTC (+5h from Chicago time)

timestr=cellfun(@(x) strsplit(x(strfind(x,'"TIMESTAMP"'):strfind(x,'"TIMESTAMP"')+100),{':' ','}),Payload,'UniformOutput',false);
timestamp=cellfun(@(x) str2double(x{:,2}),timestr);
timestamp=datetime(1970,1,1,0,0,timestamp);

filtered_times=timestamp>cutoff_start & timestamp<cutoff_end;

%indices of correct and failure records
correct = cellfun(@(x) isempty(strfind(x,'Failure')),Probe);
ind_correct = find(correct & filtered_times);
ind_failures = find(~correct & filtered_times);

%Initialize structure contents (each element corresponds to a window)
winsize = -ones(length(ind_correct),1);
features = -ones(length(ind_correct),Nfeatures); %fallnet model features
timestampSTART_END=-ones(length(ind_correct),2); %start end of each window 
evalstart = -ones(length(ind_correct),1);  %timestamp when the model started to evaluate  
sensor_counts =-ones(length(ind_correct),3); % # of samples per sensor
duration = -ones(length(ind_correct),3);  %evaluation time for each record (preparation, verification, evaluation)  
acce = cell(length(ind_correct),1);      %sensor data for each window
gyro = cell(length(ind_correct),1);
baro = cell(length(ind_correct),1);
acce_inerp = cell(length(ind_correct),1);      %inerpolated sensor data for each window
gyro_inerp = cell(length(ind_correct),1);
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
    
    v=find(strcmp(temp_correct,'FALL_VALUES'));
    features(ind,:) = str2double(temp_correct(v+1:v+Nfeatures));
    
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
    
    % parse interpolated sensors data
    acc_inerp_start = find(strcmp(temp_correct,'INTERPOLATED_ACCELEROMETER_SAMPLES'));
    gyr_inerp_start = find(strcmp(temp_correct,'INTERPOLATED_GYROSCOPE_SAMPLES'));
    
    ix_int_acc = ix(find(ix > acc_inerp_start,1));
    iy_int_acc = iy(find(iy > acc_inerp_start,1));
    iz_int_acc = iz(find(iz > acc_inerp_start,1));
    sts_int_acc = sensortimestamps(find(sensortimestamps > acc_inerp_start,1));
    Nsamplesacc_int = sort([ix_int_acc iy_int_acc iz_int_acc sts_int_acc]);
    Nsamplesacc_int = Nsamplesacc_int(2)-Nsamplesacc_int(1)-1;
    acce_inerp{ind} = [str2double(temp_correct(sts_int_acc+1:sts_int_acc+Nsamplesacc_int))' str2double(temp_correct(ix_int_acc+1:ix_int_acc+Nsamplesacc_int))' ...
        str2double(temp_correct(iy_int_acc+1:iy_int_acc+Nsamplesacc_int))' str2double(temp_correct(iz_int_acc+1:iz_int_acc+Nsamplesacc_int))'];
    
    ix_int_gyr = ix(find(ix > gyr_inerp_start,1));
    iy_int_gyr = iy(find(iy > gyr_inerp_start,1));
    iz_int_gyr = iz(find(iz > gyr_inerp_start,1));
    sts_int_gyr = sensortimestamps(find(sensortimestamps > gyr_inerp_start,1));
    Nsamplesgyr_int = sort([ix_int_gyr iy_int_gyr iz_int_gyr sts_int_gyr]);
    Nsamplesgyr_int = Nsamplesgyr_int(2)-Nsamplesgyr_int(1)-1;
    gyro_inerp{ind} = [str2double(temp_correct(sts_int_gyr+1:sts_int_gyr+Nsamplesgyr_int))' str2double(temp_correct(ix_int_gyr+1:ix_int_gyr+Nsamplesgyr_int))' ...
        str2double(temp_correct(iy_int_gyr+1:iy_int_gyr+Nsamplesgyr_int))' str2double(temp_correct(iz_int_gyr+1:iz_int_gyr+Nsamplesgyr_int))'];
end

labels.winsize = winsize;
labels.features = features; %fallnet model features
labels.timestampSTART_END = timestampSTART_END; %fallnet model features
labels.evalstart = evalstart;  %timestamp when the model started to evaluate  
labels.sensor_counts = sensor_counts; % # of samples per sensor
labels.duration = duration;  %evaluation time for each record (preparation, verification, evaluation)  
labels.acce = acce;      %sensor data for each window
labels.gyro = gyro;
labels.baro = baro;
labels.acce_inerp = acce_inerp;      %interpolated sensor data for each window
labels.gyro_inerp = gyro_inerp;


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

save FallProbe_TestData labels
