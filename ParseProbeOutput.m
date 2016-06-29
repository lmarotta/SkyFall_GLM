clear all

[FileName,PathName,~] = uigetfile('*.txt');
Payload=readtable([PathName '\' FileName],'Delimiter','\t');
Payload=Payload.Payload;

fall_probe=[];
acc_probe=[];
gyr_probe=[];
bar_probe=[];

for i=1:length(Payload)
    if length(Payload{i})<9
        continue
    end
    
    temp=strsplit(Payload{i},{'{' '"' ':' ',' ' ' '[' ']' '}'});
    if strcmp(Payload{i}(3:9),'IS_FALL')
        if size(temp,2)==size(fall_probe,2) || size(fall_probe,2)==0
            fall_probe=[fall_probe; temp];
        end
    end
    if strcmp(Payload{i}(71:79),'Gyroscope')
        gyr_probe=[gyr_probe; temp];
    end
    if strcmp(Payload{i}(71:83),'Accelerometer')
        acc_probe=[acc_probe; temp];
    end
    if strcmp(Payload{i}(71:78),'Pressure')
        bar_probe=[bar_probe; temp];
    end
    
end

Data_type='accelerometer';
y=find(strcmp(acc_probe(1,:),'Y'));
x=find(strcmp(acc_probe(1,:),'X'));
len=x-y-1;
X=[];
Y=[];
Z=[];
NormalizedTimestamp=[];
for i=1:size(acc_probe,1)
    Y=[Y; str2double(acc_probe(i,y+1:y+len)).'];
    X=[X; str2double(acc_probe(i,y+len+2:y+2*len+1)).'];
    Z=[Z; str2double(acc_probe(i,y+2*len+3:y+3*len+2)).'];
    NormalizedTimestamp=[NormalizedTimestamp; str2double(acc_probe(i,y+3*len+21:y+4*len+20)).'];
end

save accTest NormalizedTimestamp X Y Z Data_type

Data_type='Gyro';
y=find(strcmp(gyr_probe(1,:),'Y'));
x=find(strcmp(gyr_probe(1,:),'X'));
e=find(strcmp(gyr_probe(1,:),'EVENT_TIMESTAMP'));
len=x-y-1;
X=[];
Y=[];
Z=[];
NormalizedTimestamp=[];
for i=1:size(gyr_probe,1)
    Y=[Y; str2double(gyr_probe(i,y+1:y+len)).'];
    X=[X; str2double(gyr_probe(i,y+len+2:y+2*len+1)).'];
    Z=[Z; str2double(gyr_probe(i,y+2*len+3:y+3*len+2)).'];
    NormalizedTimestamp=[NormalizedTimestamp; str2double(gyr_probe(i,e+1:e+len)).'];
end

save gyrTest NormalizedTimestamp X Y Z Data_type

Data_type='barometer';
y=find(strcmp(bar_probe(1,:),'EVENT_TIMESTAMP'));
x=find(strcmp(bar_probe(1,:),'PRESSURE'));
a=find(strcmp(bar_probe(1,:),'ALTITUDE'));
len=x-y-1;
Altitude=[];
Pressure=[];
tNormalizedTimestamp=[];
for i=1:size(bar_probe,1)
    Altitude=[Altitude; str2double(bar_probe(i,a+1:a+len)).'];
    Pressure=[Pressure; str2double(bar_probe(i,x+1:x+len)).'];
    tNormalizedTimestamp=[tNormalizedTimestamp; str2double(bar_probe(i,y+1:y+len)).'];
end

save barTest tNormalizedTimestamp Altitude Pressure Data_type

x=find(strcmp(fall_probe(1,:),'TIMESTAMP'));
v=find(strcmp(fall_probe(1,:),'FALL_VALUES'));
labels.timestamp=str2double(fall_probe(:,x+1));
labels.values=str2double(fall_probe(:,v+1:v+43));
save FallProbe_TestData labels