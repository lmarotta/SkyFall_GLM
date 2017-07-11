%%  ParseProbeOutputPlusLabels
% Inputs:
% plot_data = 1; - flag to plot parsed data
% labels_offset = 0; offset for time discrepancy between labels and data
%
% Loads a user-selected .txt file containing exported data from Fallnet
% and sensor probes and separates probe data for analysis

function ParseProbeOutputPlusLabels(plot_data, labels_offset)
% clear all

if nargin == 0
    labels_offset = 0; % offset for data from 10-24-2016: 396000; 0 otherwise
    plot_data = 1; % flag to plot parsed data
elseif nargin ~= 2
    error('invalid number of input arguments. Should be 0 or 2');
end

% load file
[FileName,PathName,~] = uigetfile('*.txt');
PayloadTable=readtable([PathName '/' FileName],'Delimiter','\t');
Probe = PayloadTable.Probe;
Payload=PayloadTable.Payload;
PhoneID = PayloadTable.UserID;

date_of_data_collection = PayloadTable.Logged{1};
date_of_data_collection = strcat(date_of_data_collection(3:4), date_of_data_collection(6:7), date_of_data_collection(9:10)); 

% split labels and data
fall_labels_ind = cellfun(@(x) ~isempty(strfind(x,'falls')),Payload);
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
fall_laytime = -ones(length(fall_labels),1);

for i=1:length(fall_labels)
    temp=strsplit(fall_labels{i},{'{' '"' ':' ',' '[' ']' '}'});
    
    ts_start = find(strcmp(temp,'start'));
    ts_end = find(strcmp(temp,'end'));
    fall_start_end(i,:) = [str2double(temp(ts_start+1)) str2double(temp(ts_end+1))];
    
    type = find(strcmp(temp,'name'));
    fall_types(i) = temp(type+2);
    
    subj = find(strcmp(temp,'SubID'));
    fall_subject(i) = temp(subj+2);
    
    loc = find(strcmp(temp,'location'));
    fall_location(i) = temp(loc+2);
    
    laytime = find(strcmp(temp,'lietime'));
    fall_laytime(i) = str2double(temp(laytime+2));
end

fall_start_end = fall_start_end - labels_offset;

falllabels.types = fall_types;
falllabels.start_end_marked = fall_start_end;
falllabels.location = fall_location;
falllabels.subject = fall_subject;
falllabels.laytime = fall_laytime;

%indices of correct and failure records
correct = cellfun(@(x) isempty(strfind(x,'Failure')) && ~isempty(strfind(x,'FallNet')),data_probe);
ind_correct = find(correct);
failure = cellfun(@(x) ~isempty(strfind(x,'Failure')),data_probe);
ind_failures = find(failure);

%Indices of Activity Detection records
act_inds = find(cellfun(@(x) ~isempty(strfind(x,'ActivityDetection')),Probe));

%Initialize structure contents (each element corresponds to a window)
winsize = -ones(length(ind_correct),1);
timestampSTART_END=-ones(length(ind_correct),2); %start end of each window 
sensor_counts =-ones(length(ind_correct),3); % # of samples per sensor 
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
    
    %Additional Start and End timestamps for the clip
    xs = find(strcmp(temp_correct,'MIN_ABSOLUTE_TIMESTAMP'));
    xe = find(strcmp(temp_correct,'MAX_ABSOLUTE_TIMESTAMP'));
    timestampSTART_END(ind,:) = [str2double(temp_correct(xs+1)) str2double(temp_correct(xe+1))];
    
    ac = find(strcmp(temp_correct,'ACCELEROMETER_READING_COUNT'));
    gc = find(strcmp(temp_correct,'GYROSCOPE_READING_COUNT'));
    bc = find(strcmp(temp_correct,'BAROMETER_READING_COUNT'));
    sensor_counts(ind,:)=[str2double(temp_correct(ac+1)) str2double(temp_correct(gc+1)) str2double(temp_correct(bc+1))];
    
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
labels.sensor_counts = sensor_counts; % # of samples per sensor 
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

%% label falls data
if ~isempty(fall_labels)
    
    keep_ind = zeros(length(labels.winsize),1);
    location = {};
    subject = {};
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
data.timestampSTART_END = labels.timestampSTART_END(keep_ind, true(1,2)); %fallnet model features
data.sensor_counts = labels.sensor_counts(keep_ind, true(1,3)); % # of samples per sensor
data.acce = labels.acce(keep_ind);      %sensor data for each window
data.gyro = labels.gyro(keep_ind);
data.baro = labels.baro(keep_ind);

data.type_str = type_str;
data.subject = subject;
data.location = location;

data.failure.timestampSTART_END = labels.failure.timestampSTART_END(keep_failure, true(1,2)); %start and end timestamps for sensor data in the window
data.failure.reason = labels.failure.reason(keep_failure);
data.failure.value = labels.failure.value(keep_failure);
data.failure.duration = labels.failure.duration(keep_failure, true(1,3)); 
data.failure.evalstart = labels.failure.evalstart(keep_failure); %timestamp when the model started to evaluate 

% Values field for data
data.value=get_value(data.type_str);

save(['falls_data_' date_of_data_collection], 'data')


%% Activity Detection Probe
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

save(['ActData_' date_of_data_collection], 'activity_detection')

combine_data_into_10sec_clips(date_of_data_collection, plot_data)







