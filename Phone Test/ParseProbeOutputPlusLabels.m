%%  ParseProbeOutputPlusLabels
% Inputs:
% plot_data = 1; - flag to plot parsed data
% labels_offset = 0; offset for data from 10-24-2016: 396000; 0 otherwise
% optional:
% YYYY_start, MM_start, DD_start, HH_start, MIN_start - data cutoff start
% in UTC timestamp (+5h from Chicago time)
% YYYY_end, MM_end, DD_end, HH_end, MIN_end  - data cutoff end in UTC
% timestamp (+5h from Chicago time)
%
% Loads a user-selected .txt file containing exported data from Fallnet
% and sensor probes and separates probe data for analysis

function ParseProbeOutputPlusLabels(plot_data, labels_offset, YYYY_start, MM_start, DD_start, HH_start, MIN_start, YYYY_end, MM_end, DD_end, HH_end, MIN_end)
% clear all

if nargin == 0
    labels_offset = 0; % offset for data from 10-24-2016: 396000; 0 otherwise
    plot_data = 1; % flag to plot parsed data
elseif nargin ~= 2 && nargin ~= 12
    error('invalid number of input arguments. Should be 0, 2 or 12');
end

% load file
[FileName,PathName,~] = uigetfile('*.txt');
PayloadTable=readtable([PathName '/' FileName],'Delimiter','\t');
Probe = PayloadTable.Probe;
Payload=PayloadTable.Payload;
PhoneID = PayloadTable.UserID;

date_of_data_collection = PayloadTable.Logged{1};
date_of_data_collection = strcat(date_of_data_collection(3:4), date_of_data_collection(6:7), date_of_data_collection(9:10));

%load current fall model file from github
[FileName,PathName,~] = uigetfile('class_params*.mat');
model=load([PathName '/' FileName]);
Nfeatures = length(model.b);

if nargin == 12
    cutoff_start = datetime(YYYY_start,MM_start,DD_start,HH_start,MIN_start,0,0);  %This is in UTC (+5h from Chicago time)
    cutoff_end = datetime(YYYY_end,MM_end,DD_end,HH_end,MIN_end,0,0);  %This is in UTC (+5h from Chicago time)

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
fall_labels = Payload(fall_labels_ind);  % if there are mislabeled falls, delete from this array

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

fall_start_end = fall_start_end - labels_offset;

falllabels.types = fall_types;
falllabels.start_end_marked = fall_start_end;
falllabels.location = fall_location;
falllabels.subject = fall_subject;

save falllabels falllabels


% parse activities labels
activity_types = cell(length(activity_labels),1);
activity_subject = cell(length(activity_labels),1);
activity_location = cell(length(activity_labels),1);
activity_start_end = -ones(length(activity_labels),2);

for i=1:length(activity_labels)
    temp=strsplit(activity_labels{i},{'{' '"' ':' ',' '[' ']' '}'});
    
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

activity_start_end=activity_start_end-labels_offset;

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
    
    % parse interpolated sensors data
    %{
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
    %}
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

save fall_data_unlabeled labels

%% label activities data
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
                %curr_clip_start = labels.evalstart(j);        
                %curr_clip_end = curr_clip_start + labels.winsize(j)*1000;
                curr_clip_start = labels.timestampSTART_END(j,1);        
                curr_clip_end = labels.timestampSTART_END(j,2);

                if fall_start > curr_clip_start && fall_start < curr_clip_end
                    
                    % pull previous clip 
                    %{
                    keep_ind(j-1) = 1;
                    type_str = [type_str; curr_type];
                    location = [location; curr_loc];
                    subject = [subject; curr_subj];
                    %}

                    %pull clip where start of the fall is marked
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
                        
                        if plot_data
                            % plot 1 fall split into 2 clips for pouch location
                            if strcmp(curr_loc, 'pouch')
                                clip0 = labels.gyro{j-1};
                                clip1 = labels.gyro{j};
                                clip2 = labels.gyro{j+1}; 
                                acc_clip0 = labels.acce{j-1};
                                acc_clip1 = labels.acce{j};
                                acc_clip2 = labels.acce{j+1}; 
                                %t1 = 1:size(clip1,1);
                                %t2 = 1:size(clip2,1);
                                prev_clip_start = labels.timestampSTART_END(j-1,1);
                                next_clip_start = labels.timestampSTART_END(j+1,1);
                                t0 = clip0(:,1)*1000 + prev_clip_start;
                                t1 = clip1(:,1)*1000 + curr_clip_start;
                                t2 = clip2(:,1)*1000 + next_clip_start;

                                acc_t0 = acc_clip0(:,1)*1000 + prev_clip_start;
                                acc_t1 = acc_clip1(:,1)*1000 + curr_clip_start;
                                acc_t2 = acc_clip2(:,1)*1000 + next_clip_start;

                                plot_title = strcat(curr_type,{' '},curr_loc);

                                figure, subplot(3,2,1), plot(t0,clip0(:,2), t0,clip0(:,3), t0,clip0(:,4)), legend('X','Y','Z')
                                title(plot_title)
                                subplot(3,2,3), plot(t1,clip1(:,2), t1,clip1(:,3), t1,clip1(:,4))                            
                                y1=get(gca,'ylim'); hold on, plot([fall_start fall_start],y1)
                                subplot(3,2,5), plot(t2,clip2(:,2), t2,clip2(:,3), t2,clip2(:,4))
                                y1=get(gca,'ylim'); hold on, plot([fall_end fall_end],y1)

                                subplot(3,2,2), plot(acc_t0,acc_clip0(:,2), acc_t0,acc_clip0(:,3), acc_t0,acc_clip0(:,4))
                                subplot(3,2,4), plot(acc_t1,acc_clip1(:,2), acc_t1,acc_clip1(:,3), acc_t1,acc_clip1(:,4))                            
                                y1=get(gca,'ylim'); hold on, plot([fall_start fall_start],y1)
                                subplot(3,2,6), plot(acc_t2,acc_clip2(:,2), acc_t2,acc_clip2(:,3), acc_t2,acc_clip2(:,4))
                                y1=get(gca,'ylim'); hold on, plot([fall_end fall_end],y1)
                            end
                        end
                    else
                        if plot_data
                            if strcmp(curr_loc, 'pouch')
                                clip0 = labels.gyro{j-1};
                                clip = labels.gyro{j};
                                acc_clip0 = labels.acce{j-1};
                                acc_clip = labels.acce{j};
                                prev_clip_start = labels.timestampSTART_END(j-1,1);
                                t0 = clip0(:,1)*1000 + prev_clip_start;
                                t = clip(:,1)*1000 + curr_clip_start;

                                acc_t0 = acc_clip0(:,1)*1000 + prev_clip_start;
                                acc_t = acc_clip(:,1)*1000 + curr_clip_start;
                                plot_title = strcat(curr_type,{' '},curr_loc);

                                figure, subplot(2,2,1), plot(t0,clip0(:,2), t0,clip0(:,3), t0,clip0(:,4)), legend('X','Y','Z')
                                title(plot_title)
                                subplot(2,2,3), plot(t,clip(:,2), t,clip(:,3), t,clip(:,4))                           
                                y1=get(gca,'ylim'); hold on, plot([fall_start fall_start],y1)
                                y1=get(gca,'ylim'); hold on, plot([fall_end fall_end],y1)

                                subplot(2,2,2), plot(acc_t0,acc_clip0(:,2), acc_t0,acc_clip0(:,3), acc_t0,acc_clip0(:,4))
                                subplot(2,2,4), plot(acc_t,acc_clip(:,2), acc_t,acc_clip(:,3), acc_t,acc_clip(:,4))                           
                                y1=get(gca,'ylim'); hold on, plot([fall_start fall_start],y1)
                                y1=get(gca,'ylim'); hold on, plot([fall_end fall_end],y1)
                            end
                        end
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

%%
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

% Values field for data
data.value=get_value(data.type_str);

save falls_data data

%% Identify activities data to include with falls
%close all

falls_size=sum(keep_ind);
subject={};

subjs=unique(activity_subject);
subj_counts=zeros(1,length(subjs));
for indSubj=1:length(subjs)
    Subj_inds=strcmp(subjs{indSubj},activity_subject);
    start_ind=find(Subj_inds,1);
    end_ind=find(Subj_inds,1,'last');
    start_time=activity_start_end(start_ind,1);
    end_time=activity_start_end(end_ind,2);
    
    keep_ind=keep_ind | (start_time<labels.timestampSTART_END(:,2) & end_time>labels.timestampSTART_END(:,1));
    subj_counts(indSubj)=sum(keep_ind)-falls_size-sum(subj_counts(1:indSubj-1));
end

data.winsize = labels.winsize(keep_ind);
data.features = labels.features(keep_ind); %fallnet model features
data.timestampSTART_END = labels.timestampSTART_END(keep_ind); %fallnet model features
data.evalstart = labels.evalstart(keep_ind);  %timestamp when the model started to evaluate  
data.sensor_counts = labels.sensor_counts(keep_ind); % # of samples per sensor
data.duration = labels.duration(keep_ind);  %evaluation time for each record (preparation, verification, evaluation)  
data.acce = labels.acce(keep_ind);      %sensor data for each window
data.gyro = labels.gyro(keep_ind);
data.baro = labels.baro(keep_ind);

act_size=sum(keep_ind)-falls_size;

data.subject = [data.subject; cell(act_size,1)];
for i=1:length(subjs)
    data.subject(falls_size+1+sum(subj_counts(1:i-1)):falls_size+sum(subj_counts(1:i)))=subjs(i);
end
data.type_str(falls_size+1:falls_size+act_size) = ({'Non-Fall'});
data.location(falls_size+1:falls_size+act_size) = ({'NA'});

data.value=[data.value; repmat(9,[act_size 1])];

save falls_act_data data

%% label activities data
if ~isempty(activity_labels)
    
    keep_ind = zeros(length(labels.winsize),1);
    location = {};
    subject = {};
    type = [];
    type_str = {};
    
    keep_failure = zeros(length(ind_failures));
    
    curr_data_ind = 1;
    for i=1:length(activity_labels)
        % current labels data
        curr_act_start = activity_start_end(i,1);
        curr_act_end = activity_start_end(i,2);
        curr_type = activity_types(i);
        curr_loc = activity_location(i);
        curr_subj = activity_subject(i);        
      
        act_found = 0;
        for j=curr_data_ind:length(ind_correct)
            
            if ~act_found
                % current sensors data
                %curr_clip_start = labels.evalstart(j);        
                %curr_clip_end = curr_clip_start + labels.winsize(j)*1000;
                curr_clip_start = labels.timestampSTART_END(j,1);        
                curr_clip_end = labels.timestampSTART_END(j,2);

                if curr_act_start > curr_clip_start && curr_act_start < curr_clip_end
                    % beginning of the activity
                    ind_act_start = j+1; 
                end
                
                if curr_act_end > curr_clip_start && curr_act_end < curr_clip_end
                    % end of the activity
                    ind_act_end = j-1;
                    
                    keep_ind(ind_act_start:ind_act_end) = 1;
                    num_clips = ind_act_end - ind_act_start + 1;
                    ct = cell(num_clips,1);
                    ct(:) = {curr_type};
                    cl = cell(num_clips,1);
                    cl(:) = {curr_loc};
                    cs = cell(num_clips,1);
                    cs(:) = {curr_subj};
                    
                    type_str = [type_str; ct];
                    location = [location; cl];
                    subject = [subject; cs];
                    
                    curr_data_ind = j;
                    act_found = 1;
                    
                    % plot activity
                    if plot_data
                        data_to_plot = labels.acce(ind_act_start:ind_act_end);
                        % convert timestamps to absolute
                        for k=1:length(data_to_plot)
                            data_to_plot{k}(:,1) = data_to_plot{k}(:,1)*1000 + labels.timestampSTART_END(ind_act_start+k,1);
                        end                
                        data_to_plot = cell2mat(data_to_plot);
                        t = data_to_plot(:,1);
                        plot_title = strcat(curr_type,{' '},curr_loc);
                        figure, plot(t,data_to_plot(:,2), t,data_to_plot(:,3), t,data_to_plot(:,4)), legend('X','Y','Z')
                        title(plot_title);
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
            if curr_clip_start > curr_act_start && curr_clip_start < curr_act_end
                keep_failure(j) = 1;
            end
        end
        
        
    end
end
%%
keep_ind = logical(keep_ind);
keep_failure = logical(keep_failure);

actdata.winsize = labels.winsize(keep_ind);
actdata.features = labels.features(keep_ind); %fallnet model features
actdata.timestampSTART_END = labels.timestampSTART_END(keep_ind); %fallnet model features
actdata.evalstart = labels.evalstart(keep_ind);  %timestamp when the model started to evaluate  
actdata.sensor_counts = labels.sensor_counts(keep_ind); % # of samples per sensor
actdata.duration = labels.duration(keep_ind);  %evaluation time for each record (preparation, verification, evaluation)  
actdata.acce = labels.acce(keep_ind);      %sensor data for each window
actdata.gyro = labels.gyro(keep_ind);
actdata.baro = labels.baro(keep_ind);

actdata.type_str = type_str;
actdata.subject = subject;
actdata.location = location;

actdata.failure.timestampSTART_END = labels.failure.timestampSTART_END(keep_failure); %start and end timestamps for sensor data in the window
actdata.failure.reason = labels.failure.reason(keep_failure);
actdata.failure.value = labels.failure.value(keep_failure);
actdata.failure.duration = labels.failure.duration(keep_failure); 
actdata.failure.evalstart = labels.failure.evalstart(keep_failure); %timestamp when the model started to evaluate 

actdata.value=get_value(actdata.type_str);

save activities_data actdata

combine_data_into_10sec_clips(date_of_data_collection, plot_data)







