% Loads a user-selected .txt file containing exported data from Fallnet
% and sensor probes and separates probe data fro analysis

clear all

% load file
[FileName,PathName,~] = uigetfile('*.txt');
Payload=readtable([PathName '\' FileName],'Delimiter','\t');
Payload=Payload.Payload;

fall_probe=[];
acc_probe=[];
gyr_probe=[];
bar_probe=[];

for i=1:length(Payload)
    % skip junk data
    if length(Payload{i})<9
        continue
    end
    
    temp=strsplit(Payload{i},{'{' '"' ':' ',' ' ' '[' ']' '}'});
    % for files exported after 07/12/16
    % if strcmp(Payload{i}(3:9),'fallnet')
    
    % for files exported after 07/12/16
    if strcmp(Payload{i}(71:77),'fallnet')
        % only stores fall_probe data of the same size
        % prevents errors if data contains probe data for 
        % multiple sets of model parameters
        if size(temp,2)==size(fall_probe,2) || size(fall_probe,2)==0
            fall_probe=[fall_probe; temp];
        end
    elseif strcmp(Payload{i}(71:79),'Gyroscope')
        gyr_probe=[gyr_probe; temp];
    elseif strcmp(Payload{i}(71:83),'Accelerometer')
        acc_probe=[acc_probe; temp];
    elseif strcmp(Payload{i}(71:78),'Pressure')
        bar_probe=[bar_probe; temp];
    else
        dummy=1;
    end
    
end

% combine data from each sensor into single vectors

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

% separate timestamps and feature values from fall_probe and save in
% labels structure

if ~isempty(fall_probe)
    x=find(strcmp(fall_probe(1,:),'TIMESTAMP'));
    v=find(strcmp(fall_probe(1,:),'FALL_VALUES'));
    %Additional Start and End timestamps for the clip
    xs = find(strcmp(fall_probe(1,:),'EVALUATION_WINDOW_START'));
    xe = find(strcmp(fall_probe(1,:),'EVALUATION_WINDOW_END'));

    labels.timestamp=str2double(fall_probe(:,x+1));
    labels.values=str2double(fall_probe(:,v+1:v+175));
    labels.timestampSTART_END=[str2double(fall_probe(:,xs+1)) str2double(fall_probe(:,xe+1))]; 

    %parse sample counts for each probe
    ac = find(strcmp(fall_probe(1,:),'ACCELEROMETER_READING_COUNT'));
    gc = find(strcmp(fall_probe(1,:),'GYROSCOPE_READING_COUNT'));
    bc = find(strcmp(fall_probe(1,:),'BAROMETER_READING_COUNT'));
    labels.sensor_counts=[str2double(fall_probe(:,ac+1)) str2double(fall_probe(:,gc+1)) str2double(fall_probe(:,bc+1))]; 

    save FallProbe_TestData labels

    %histogram of duration of clips
    td = (labels.timestampSTART_END(:,2)-labels.timestampSTART_END(:,1))/1000;
    figure, histogram(td), xlabel('Clip Duration [s]'), ylabel('Frequency of clips')
    figure, hold on, subplot(311), histogram(labels.sensor_counts(:,1)), title('acc'), xlabel('# of samples in clip')
    subplot(312), histogram(labels.sensor_counts(:,2)), title('gyr')
    subplot(313), histogram(labels.sensor_counts(:,3)), title('bar')
end