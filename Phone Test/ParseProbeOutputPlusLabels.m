% Loads a user-selected .txt file containing exported data from Fallnet
% and sensor probes and separates probe data fro analysis

clear all

% flags
filter_by_time = 0;

% load file
[FileName,PathName,~] = uigetfile('*.txt');
PayloadTable=readtable([PathName '/' FileName],'Delimiter','\t');
Probe = PayloadTable.Probe;
Payload=PayloadTable.Payload;
PhoneID = PayloadTable.UserID;

%load current fall model file from github
[FileName,PathName,~] = uigetfile('class_params*.mat');
model=load([PathName '/' FileName]);
Nfeatures = length(model.b);

if filter_by_time
    cutoff_start = datetime(2016,10,17,15,40,0,0);  %This is in UTC (+5h from Chicago time)
    cutoff_end = datetime(2016,10,17,17,10,0,0);  %This is in UTC (+5h from Chicago time)

    timestr=cellfun(@(x) strsplit(x(strfind(x,'"TIMESTAMP"'):strfind(x,'"TIMESTAMP"')+100),{':' ','}),Payload,'UniformOutput',false);
    timestamp=cellfun(@(x) str2double(x{:,2}),timestr);
    timestamp=datetime(1970,1,1,0,0,timestamp);

    filtered_times=timestamp>cutoff_start & timestamp<cutoff_end;
    Probe = Probe(filtered_times);
    Payload = Payload(filtered_times);
    PhoneID = PhoneID(filtered_times);
end    

% split labels and data
fall_labels_ind = cellfun(@(x) strcmp(x,'falls'),Probe);
fall_labels = Payload(fall_labels_ind);

activity_labels_ind = cellfun(@(x) strcmp(x,'activities'),Probe);
activity_labels = Payload(activity_labels_ind);

labels_phone_id = PhoneID(fall_labels_ind);
labels_phone_id = labels_phone_id{1};

data_ind = cellfun(@(x) ~strcmp(x, labels_phone_id),PhoneID);
data_fallnet = Payload(data_ind);
data_probe = Probe(data_ind);

% parse falls labels
fall_types = cell(length(fall_labels),1);
fall_subject = cell(length(fall_labels),1);
fall_location = cell(length(fall_labels),1);
fall_start_end = -ones(length(fall_labels),2);

for i=1:length(fall_labels)
    temp=strsplit(fall_labels{i},{'{' '"' ':' ',' '[' ']' '}'});
    
    ts_start = find(strcmp(temp,'start'));
    ts_end = find(strcmp(temp,'end'));
    fall_start_end(i,:) = [str2double(temp(ts_start+1)) str2double(temp(ts_end+1))];
    
    type = find(strcmp(temp,'name'));
    fall_types(i) = temp(type+2);
    
    subj = find(strcmp(temp,'Subject ID'));
    fall_subject(i) = temp(subj+2);
    
    loc = find(strcmp(temp,'phone location'));
    fall_location(i) = temp(loc+2);
end

% parse activities labels
activity_types = cell(length(activity_labels),1);
activity_subject = cell(length(activity_labels),1);
activity_location = cell(length(activity_labels),1);
activity_start_end = -ones(length(activity_labels),2);

for i=1:length(activity_labels)
    temp=strsplit(activity_labels{i},{'{' '"' ':' ',' ' ' '[' ']' '}'});
    
    ts_start = find(strcmp(temp,'start'));
    ts_end = find(strcmp(temp,'end'));
    activity_start_end(i,:) = [str2double(temp(ts_start+1)) str2double(temp(ts_end+1))];
    
    type = find(strcmp(temp,'name'));
    activity_types(i) = temp(type+2);
    
    subj = find(strcmp(temp,'Subject ID'));
    activity_subject(i) = temp(subj+2);
    
    loc = find(strcmp(temp,'phone location'));
    activity_location(i) = temp(loc+2);
end


%indices of correct and failure records
correct = cellfun(@(x) isempty(strfind(x,'Failure')),data_probe);
ind_correct = find(correct);
ind_failures = find(~correct);

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
%FAILURES
failure_timestampSTART_END = -ones(length(ind_failures),2); %start and end timestamps for sensor data in the window
failure_reason = cell(length(ind_failures),1);
failure_value = -ones(length(ind_failures),1);
failure_duration = -ones(length(ind_failures),3); 
failure_evalstart = -ones(length(ind_failures),1); %timestamp when the model started to evaluate 

for ind=1:length(ind_correct)
    
    i_correct = ind_correct(ind); %index
    temp_correct=strsplit(data_fallnet{i_correct},{'{' '"' ':' ',' ' ' '[' ']' '}'}); 
 
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
labels.features = features; %fallnet model features
labels.timestampSTART_END = timestampSTART_END; %fallnet model features
labels.evalstart = evalstart;  %timestamp when the model started to evaluate  
labels.sensor_counts = sensor_counts; % # of samples per sensor
labels.duration = duration;  %evaluation time for each record (preparation, verification, evaluation)  
labels.acce = acce;      %sensor data for each window
labels.gyro = gyro;
labels.baro = baro;


for ind=1:length(ind_failures)
    i = ind_failures(ind); %index
    temp=strsplit(data_fallnet{i},{'{' '"' ':' ',' ' ' '[' ']' '}'});
    
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

% label activities data
if ~isempty(fall_labels)
    
    keep_ind = zeros(length(labels.winsize),1);
    location = {};
    subject = {};
    type = [];
    type_str = {};
    
    keep_failure = zeros(length(ind_failures));
    
    curr_data_ind = 1;
    for i=1:length(fall_labels)
        % current labels data
        fall_start = fall_start_end(i,1);
        fall_end = fall_start_end(i,2);
        curr_type = fall_types(i);
        curr_loc = fall_location(i);
        curr_subj = fall_subject(i);        
      
        clip_found = 0;
        for j=curr_data_ind:length(ind_correct)
            
            if ~clip_found
                % current sensors data
                curr_clip_start = labels.evalstart(j);        
                curr_clip_end = curr_clip_start + labels.winsize(j)*1000;

                if fall_start > curr_clip_start && fall_start < curr_clip_end

                    keep_ind(j) = 1;
                    type_str = [type_str; curr_type];
                    location = [location; curr_loc];
                    subject = [subject; curr_subj];                

                    clip_found = 1;
                    curr_data_ind = j+1;

                    if fall_end > curr_clip_end
                        curr_data_ind = j+2;

                        keep_ind(j+1) = 1;
                        type_str = [type_str; curr_type];
                        location = [location; curr_loc];
                        subject = [subject; curr_subj];                    

                    end

                end
            else
                break;
            end 
        end
        
        
        % look for failures in the labeled range (aasumes there are not a
        % lot of failures)     
        for j=1:length(ind_failures)
            curr_clip_start = labels.failure.evalstart(j);
            curr_clip_end = curr_clip_start + 5000;
            if fall_start > curr_clip_start && fall_start < curr_clip_end
                keep_failure(j) = 1;
            end
        end
        
        
    end
end

keep_ind = logical(keep_ind);
keep_failure = logical(keep_failure);

data.winsize = labels.winsize(keep_ind);
data.features = labels.features(keep_ind); %fallnet model features
data.timestampSTART_END = labels.timestampSTART_END(keep_ind); %fallnet model features
data.evalstart = labels.evalstart(keep_ind);  %timestamp when the model started to evaluate  
data.sensor_counts = labels.sensor_counts(keep_ind); % # of samples per sensor
data.duration = labels.duration(keep_ind);  %evaluation time for each record (preparation, verification, evaluation)  
data.acce = labels.acce(keep_ind);      %sensor data for each window
data.gyro = labels.gyro(keep_ind);
data.baro = labels.baro(keep_ind);

data.type_str = type_str;
data.subject = subject;
data.location = location;

data.failure.timestampSTART_END = labels.failure.timestampSTART_END(keep_failure); %start and end timestamps for sensor data in the window
data.failure.reason = labels.failure.reason(keep_failure);
data.failure.value = labels.failure.value(keep_failure);
data.failure.duration = labels.failure.duration(keep_failure); 
data.failure.evalstart = labels.failure.evalstart(keep_failure); %timestamp when the model started to evaluate 

save falls_data data



% label activities data
%{
keep_ind = zeros(length(labels.winsize),1);

curr_act_start = activity_start_end(1,1);
curr_act_end = activity_start_end(1,2);
for i=1:length(labels.winsize)
    
    
    
    
end
%}






