% Loads a user-selected .txt file containing exported data from Fallnet
% and sensor probes and separates probe data fro analysis

clear all

% load file
[FileName,PathName,~] = uigetfile('*.txt');
Payload=readtable([PathName '\' FileName],'Delimiter','\t');
Probe = Payload.Probe;
Payload=Payload.Payload;

fall_probe=[];
acc_probe=[];
gyr_probe=[];
bar_probe=[];
fall_probe_failure_counts = [];
fall_probe_failure_gap = [];
fall_probe_failure_interp =[];
fall_probe_failure_other =[];

for i=1:length(Payload)
    % skip junk data
    if length(Payload{i})<9
        continue
    end
    
    temp=strsplit(Payload{i},{'{' '"' ':' ',' ' ' '[' ']' '}'});
    % for files exported after 07/12/16
    % if strcmp(Payload{i}(3:9),'fallnet')
    % for files exported after 07/12/16
    if strcmp(Probe{i}(end-11:end),'FallNetProbe')
        % only stores fall_probe data of the same size
        % prevents errors if data contains probe data for
        % multiple sets of model parameters
        if size(temp,2)==size(fall_probe,2) || size(fall_probe,2)==0
            fall_probe=[fall_probe; temp];
        end
    elseif strcmp(Probe{i}(60:68),'Gyroscope')
        gyr_probe=[gyr_probe; temp];
    elseif strcmp(Probe{i}(60:72),'Accelerometer')
        acc_probe=[acc_probe; temp];
    elseif strcmp(Probe{i}(60:67),'Pressure')
        bar_probe=[bar_probe; temp];
    elseif strcmp(Probe{i}(end-6:end),'Failure')
        if find(strcmp(temp,'SAMPLES_EXPECTED'))
            fall_probe_failure_counts = [fall_probe_failure_counts; temp];
        elseif find(strcmp(temp,'GAP_SIZE'))
            fall_probe_failure_gap = [fall_probe_failure_gap; temp]; %'GAP_SIZE'
        elseif find(strcmp(temp,'INTERPOLATOR_EXCEPTION'))
            fall_probe_failure_interp = [fall_probe_failure_interp; temp]; % 'interpolator failure'
        else
            fall_probe_failure_other = [fall_probe_failure_other; temp];
        end
        
    else
        dummy=1;
    end
    
end

% combine data from each sensor into single vectors

if ~isempty(acc_probe) && ~isempty(gyr_probe)
Data_type='accelerometer';
y=find(strcmp(acc_probe(1,:),'Y'));
x=find(strcmp(acc_probe(1,:),'X'));
s=find(strcmp(gyr_probe(1,:),'SENSOR_TIMESTAMP'));
e=find(strcmp(acc_probe(1,:),'EVENT_TIMESTAMP'));
len=x-y-1;
X=[];
Y=[];
Z=[];
NormalizedTimestamp=[];
eTimestamp=[];
for i=1:size(acc_probe,1)
    Y=[Y; str2double(acc_probe(i,y+1:y+len)).'];
    X=[X; str2double(acc_probe(i,y+len+2:y+2*len+1)).'];
    Z=[Z; str2double(acc_probe(i,y+2*len+3:y+3*len+2)).'];
    NormalizedTimestamp=[NormalizedTimestamp; str2double(acc_probe(i,s+1:s+len)).'];
    eTimestamp=[eTimestamp; str2double(acc_probe(i,e+1:e+len)).'];
    %     NormalizedTimestamp=[NormalizedTimestamp; str2double(acc_probe(i,y+3*len+21:y+4*len+20)).'];
end

NormalizedTimestamp=eTimestamp;
% NormalizedTimestamp=(NormalizedTimestamp-min(NormalizedTimestamp))/10^9+min(eTimestamp);

save accTest NormalizedTimestamp X Y Z Data_type

Data_type='Gyro';
y=find(strcmp(gyr_probe(1,:),'Y'));
x=find(strcmp(gyr_probe(1,:),'X'));
s=find(strcmp(gyr_probe(1,:),'SENSOR_TIMESTAMP'));
e=find(strcmp(gyr_probe(1,:),'EVENT_TIMESTAMP'));
len=x-y-1;
X=[];
Y=[];
Z=[];
NormalizedTimestamp=[];
eTimestamp=[];
for i=1:size(gyr_probe,1)
    Y=[Y; str2double(gyr_probe(i,y+1:y+len)).'];
    X=[X; str2double(gyr_probe(i,y+len+2:y+2*len+1)).'];
    Z=[Z; str2double(gyr_probe(i,y+2*len+3:y+3*len+2)).'];
    NormalizedTimestamp=[NormalizedTimestamp; str2double(gyr_probe(i,s+1:s+len)).'];
    eTimestamp=[eTimestamp; str2double(gyr_probe(i,e+1:e+len)).'];
end

NormalizedTimestamp=eTimestamp;
% NormalizedTimestamp=(NormalizedTimestamp-min(NormalizedTimestamp))/10^9+min(eTimestamp);

save gyrTest NormalizedTimestamp X Y Z Data_type

end

if ~isempty(bar_probe)
Data_type='barometer';
s=find(strcmp(bar_probe(1,:),'SENSOR_TIMESTAMP'));
x=find(strcmp(bar_probe(1,:),'PRESSURE'));
a=find(strcmp(bar_probe(1,:),'ALTITUDE'));
e=find(strcmp(bar_probe(1,:),'EVENT_TIMESTAMP'));
len=x-e-1;
Altitude=[];
Pressure=[];
tNormalizedTimestamp=[];
eTimestamp=[];
for i=1:size(bar_probe,1)
    Altitude=[Altitude; str2double(bar_probe(i,a+1:a+len)).'];
    Pressure=[Pressure; str2double(bar_probe(i,x+1:x+len)).'];
    tNormalizedTimestamp=[tNormalizedTimestamp; str2double(bar_probe(i,s+1:s+len)).'];
    eTimestamp=[eTimestamp; str2double(bar_probe(i,e+1:e+len)).'];
end

tNormalizedTimestamp=eTimestamp;
% tNormalizedTimestamp=(tNormalizedTimestamp-min(tNormalizedTimestamp))/10^9+min(eTimestamp);

save barTest tNormalizedTimestamp Altitude Pressure Data_type
end

% separate timestamps and feature values from fall_probe and save in
% labels structure

if ~isempty(fall_probe)
    x=find(strcmp(fall_probe(1,:),'EVALUATION_WINDOW_SIZE'));
    v=find(strcmp(fall_probe(1,:),'FALL_VALUES'));
    %Additional Start and End timestamps for the clip
    xs = find(strcmp(fall_probe(1,:),'MIN_ABSOLUTE_TIMESTAMP'));
    xe = find(strcmp(fall_probe(1,:),'MAX_ABSOLUTE_TIMESTAMP'));
    
    % log timestamp
    log_ts = find(strcmp(fall_probe(1,:),'TIMESTAMP'));
    eval_start_ts = find(strcmp(fall_probe(1,:),'EVALUATION_START'));
    eval_prep_ts = find(strcmp(fall_probe(1,:),'EVALUATION_PREPARED'));
    
    % 1st column - log timestamp; 2nd column - start of evaluation window; 3rd column - evaluation prepared
    labels.timestamp = [str2double(fall_probe(:,log_ts+1)) str2double(fall_probe(:,eval_start_ts+1)) str2double(fall_probe(:,eval_prep_ts+1))];

    labels.winsize=str2double(fall_probe(:,x+1));
    labels.values=str2double(fall_probe(:,v+1:v+419));

    % labels.timestamp=str2double(fall_probe(:,x+1));

    %when absolute time stamps are available for start and end
    labels.timestampSTART_END=[str2double(fall_probe(:,xs+1)) str2double(fall_probe(:,xe+1))];
    %labels.timestampSTART_END=[(str2double(fall_probe(:,x+1))-str2double(fall_probe(:,xe+1))) str2double(fall_probe(:,x+1))]; %TIMESTAMP(end) - EVALUATION_WINDOW_END (duration of window)
    
    
    %parse sample counts for each probe
    ac = find(strcmp(fall_probe(1,:),'ACCELEROMETER_READING_COUNT'));
    gc = find(strcmp(fall_probe(1,:),'GYROSCOPE_READING_COUNT'));
    bc = find(strcmp(fall_probe(1,:),'BAROMETER_READING_COUNT'));
    labels.sensor_counts=[str2double(fall_probe(:,ac+1)) str2double(fall_probe(:,gc+1)) str2double(fall_probe(:,bc+1))];
    
    % parse duration times
    prepd = find(strcmp(fall_probe(1,:),'PREPARE_DURATION'));
    verifd = find(strcmp(fall_probe(1,:),'VERIFY_DURATION')); 
    predd = find(strcmp(fall_probe(1,:),'PREDICT_DURATION'));
    labels.duration=[str2double(fall_probe(:,prepd+1)) str2double(fall_probe(:,verifd+1)) str2double(fall_probe(:,predd+1))];
   
    %{
    %% parse sensors recordings
    bar_ind_s = find(strcmp(fall_probe(1,:),'BAROMETER_SAMPLES')) + 2;
    acc_ind_s = find(strcmp(fall_probe(1,:),'ACCELEROMETER_SAMPLES')) + 2;
    gyr_ind_s = find(strcmp(fall_probe(1,:),'GYROSCOPE_SAMPLES')) + 2;
    
    bar_num_sampl;
    acc_num_sampl;
    gyr_num_sampl;
    %}
    
    
end

%Parses data from failure probe 
if ~isempty(fall_probe_failure_counts)
    x1s=find(strcmp(fall_probe_failure_counts(1,:),'MIN_ABSOLUTE_TIMESTAMP'));
    x1e=find(strcmp(fall_probe_failure_counts(1,:),'MAX_ABSOLUTE_TIMESTAMP'));
    x2=find(strcmp(fall_probe_failure_counts(1,:),'FAILURE_REASON'));
    x3=find(strcmp(fall_probe_failure_counts(1,:),'SAMPLES_COUNT'));
    x4 = find(strcmp(fall_probe_failure_counts(1,:),'PREPARE_DURATION'));
    x5log = find(strcmp(fall_probe_failure_counts(1,:),'TIMESTAMP'));
    x5evalstart = find(strcmp(fall_probe_failure_counts(1,:),'EVALUATION_START'));
    x5prep = find(strcmp(fall_probe_failure_counts(1,:),'EVALUATION_PREPARED'));
    
    labels.failure.timestamp = [str2double(fall_probe_failure_counts(:,x1s+1)) str2double(fall_probe_failure_counts(:,x1e+1))];
    labels.failure.reason = fall_probe_failure_counts(:,x2+1);
    labels.failure.value = str2double(fall_probe_failure_counts(:,x3+1));
    labels.failure.duration = str2double(fall_probe_failure_counts(:,x4+1));
    labels.failure.ts = [str2double(fall_probe_failure_counts(:,x5log+1)) str2double(fall_probe_failure_counts(:,x5evalstart+1)) str2double(fall_probe_failure_counts(:,x5prep+1))];
    
end

if ~isempty(fall_probe_failure_gap)
    x1s=find(strcmp(fall_probe_failure_gap(1,:),'MIN_ABSOLUTE_TIMESTAMP'));
    x1e=find(strcmp(fall_probe_failure_gap(1,:),'MAX_ABSOLUTE_TIMESTAMP'));
    x2=find(strcmp(fall_probe_failure_gap(1,:),'FAILURE_REASON'));
    x3=find(strcmp(fall_probe_failure_gap(1,:),'GAP_SIZE'));
    x4 = find(strcmp(fall_probe_failure_gap(1,:),'PREPARE_DURATION'));
    x5log = find(strcmp(fall_probe_failure_gap(1,:),'TIMESTAMP'));
    x5evalstart = find(strcmp(fall_probe_failure_gap(1,:),'EVALUATION_START'));
    x5prep = find(strcmp(fall_probe_failure_gap(1,:),'EVALUATION_PREPARED'));
    
    labels.failure.timestamp = [labels.failure.timestamp; str2double(fall_probe_failure_gap(:,x1s+1)) str2double(fall_probe_failure_gap(:,x1e+1))];
    labels.failure.reason = [labels.failure.reason; fall_probe_failure_gap(:,x2+1)];
    labels.failure.value = [labels.failure.value; str2double(fall_probe_failure_gap(:,x3+1))];
    labels.failure.duration = [labels.failure.duration; str2double(fall_probe_failure_gap(:,x4+1))];
    labels.failure.ts = [labels.failure.ts; str2double(fall_probe_failure_gap(:,x5log+1)) str2double(fall_probe_failure_gap(:,x5evalstart+1)) str2double(fall_probe_failure_gap(:,x5prep+1))];
    
end

if ~isempty(fall_probe_failure_interp)
    x1s=find(strcmp(fall_probe_failure_interp(1,:),'MIN_ABSOLUTE_TIMESTAMP'));
    x1e=find(strcmp(fall_probe_failure_interp(1,:),'MAX_ABSOLUTE_TIMESTAMP'));
    x2=find(strcmp(fall_probe_failure_interp(1,:),'FAILURE_REASON'));
    x4 = find(strcmp(fall_probe_failure_interp(1,:),'PREPARE_DURATION'));
    x5log = find(strcmp(fall_probe_failure_interp(1,:),'TIMESTAMP'));
    x5evalstart = find(strcmp(fall_probe_failure_interp(1,:),'EVALUATION_START'));
    x5prep = find(strcmp(fall_probe_failure_interp(1,:),'EVALUATION_PREPARED'));
    
    labels.failure.timestamp = [labels.failure.timestamp; str2double(fall_probe_failure_interp(:,x1s+1)) str2double(fall_probe_failure_interp(:,x1e+1))];
    labels.failure.reason = [labels.failure.reason; fall_probe_failure_interp(:,x2+1)];
    labels.failure.value = [labels.failure.value; zeros(size(fall_probe_failure_interp,1),1)];
    labels.failure.duration = [labels.failure.duration; str2double(fall_probe_failure_interp(:,x4+1))];
    labels.failure.ts = [labels.failure.ts; str2double(fall_probe_failure_interp(:,x5log+1)) str2double(fall_probe_failure_interp(:,x5evalstart+1)) str2double(fall_probe_failure_interp(:,x5prep+1))];
end

if ~isempty(fall_probe_failure_other)
    x1s=find(strcmp(fall_probe_failure_other(1,:),'MIN_ABSOLUTE_TIMESTAMP'));
    x1e=find(strcmp(fall_probe_failure_other(1,:),'MAX_ABSOLUTE_TIMESTAMP'));
    x2=find(strcmp(fall_probe_failure_other(1,:),'FAILURE_REASON'));
    x4 = find(strcmp(fall_probe_failure_other(1,:),'PREPARE_DURATION'));
    x5log = find(strcmp(fall_probe_failure_other(1,:),'TIMESTAMP'));
    x5evalstart = find(strcmp(fall_probe_failure_other(1,:),'EVALUATION_START'));
    x5prep = find(strcmp(fall_probe_failure_other(1,:),'EVALUATION_PREPARED'));
    
    labels.failure.timestamp = [labels.failure.timestamp; str2double(fall_probe_failure_other(:,x1s+1)) str2double(fall_probe_failure_other(:,x1e+1))];
    labels.failure.reason = [labels.failure.reason; fall_probe_failure_other(:,x2+1)];
    labels.failure.value = [labels.failure.value; zeros(size(fall_probe_failure_other,1),1)];
    labels.failure.duration = [labels.failure.duration; str2double(fall_probe_failure_other(:,x4+1))];
    labels.failure.ts = [labels.failure.ts; str2double(fall_probe_failure_other(:,x5log+1)) str2double(fall_probe_failure_other(:,x5evalstart+1)) str2double(fall_probe_failure_other(:,x5prep+1))];
end

save FallProbe_TestData labels

%Look at a specific time range
cutoff_start = datetime(2016,9,28,16,30,0,0);  %This is in UTC (+5h from Chicago time)
cutoff_end = datetime(2016,9,28,19,0,0,0);  %This is in UTC (+5h from Chicago time)

indstart = datetime(1970,1,1,0,0,0,labels.timestampSTART_END(:,1)) < cutoff_end & datetime(1970,1,1,0,0,0,labels.timestampSTART_END(:,1)) > cutoff_start;
indend = datetime(1970,1,1,0,0,0,labels.timestampSTART_END(:,2)) < cutoff_end & datetime(1970,1,1,0,0,0,labels.timestampSTART_END(:,2)) > cutoff_start;
labels.timestampSTART_END = labels.timestampSTART_END(indstart & indend,:);
labels.winsize = labels.winsize(indstart & indend,:);
labels.values = labels.values(indstart & indend,:);
labels.sensor_counts = labels.sensor_counts(indstart & indend,:);
labels.duration = labels.duration(indstart & indend,:);
labels.timestamp = labels.timestamp(indstart & indend,:);

indstart = datetime(1970,1,1,0,0,0,labels.failure.timestamp(:,1)) < cutoff_end & datetime(1970,1,1,0,0,0,labels.failure.timestamp(:,1)) > cutoff_start;
indend = datetime(1970,1,1,0,0,0,labels.failure.timestamp(:,2)) < cutoff_end & datetime(1970,1,1,0,0,0,labels.failure.timestamp(:,2)) > cutoff_start;
labels.failure.timestamp = labels.failure.timestamp(indstart & indend,:);
labels.failure.reason = labels.failure.reason(indstart & indend, :);
labels.failure.value = labels.failure.value(indstart & indend,:);
labels.failure.duration = labels.failure.duration(indstart & indend,:);
labels.failure.ts = labels.failure.ts(indstart & indend,:);


save FallProbe_TestData_cut labels

% statistics

%histogram of duration of clips
%{
td = (labels.timestampSTART_END(:,2)-labels.timestampSTART_END(:,1));
figure, histogram(td), xlabel('Clip Duration [s]'), ylabel('Frequency of clips')
figure, hold on, subplot(311), histogram(labels.sensor_counts(:,1)), title('acc'), xlabel('# of samples in clip')
subplot(312), histogram(labels.sensor_counts(:,2)), title('gyr')
subplot(313), histogram(labels.sensor_counts(:,3)), title('bar')
%}

% histogram of time for successful prediction
figure, hold on, subplot(411), histogram(labels.duration(:,1),10), title('Success Preparation'), xlabel('duration')
subplot(412), histogram(labels.duration(:,2),10), title('Success Verification')
subplot(413), histogram(labels.duration(:,3),10), title('Success Prediction')
subplot(414), histogram(labels.failure.duration(:,1),10), title('Failures Preparation')

%histogram of fallnet failures
reasons = labels.failure.reason;
unique_reasons =unique(reasons,'stable');
%reasons_count=cellfun(@(x) sum(ismember(reasons,x)),unique_reasons,'un',0);

B = categorical(reasons,unique_reasons);
figure, histogram(B,'BarWidth',0.5)


%% plot start and end timestamps for failure and success 

figure
plot(labels.timestampSTART_END(:,1)-labels.timestampSTART_END(1,1),zeros(length(labels.timestampSTART_END(:,1)),1),'ob'), hold on
plot(labels.timestampSTART_END(:,2)-labels.timestampSTART_END(1,1),zeros(length(labels.timestampSTART_END(:,2)),1),'xb')

plot(labels.failure.timestamp(:,1)-labels.failure.timestamp(1,1),ones(length(labels.failure.timestamp(:,1)),1),'or'), hold on
plot(labels.failure.timestamp(:,2)-labels.failure.timestamp(1,1),ones(length(labels.failure.timestamp(:,2)),1),'xr')
ylim([-0.5 1.5])

figure, subplot(2,1,1), hold on, title('success')
histogram(labels.timestampSTART_END(:,2)-labels.timestampSTART_END(:,1))
subplot(212), hold on, title('Failed')
histogram(labels.failure.timestamp(:,2)-labels.failure.timestamp(:,1))

%% plots fallnet failures values
%{
for i=1:size(unique_reasons)
    reason = unique_reasons{i,1};
    failures_ind = strcmp(reasons,reason);
    failures_count = labels.failure.value(failures_ind);
    figure,  histogram(failures_count), title(reason)
end
%}