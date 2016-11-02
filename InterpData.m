%% InterpData
function [acce, gyro, baro]=InterpData(labels)

acce=cell(length(labels.acce),1);
gyro=cell(length(labels.acce),1);
baro=cell(length(labels.acce),1);


t=0:.02:4.98;
tbar=0:1/6:29/6;

for i=1:length(labels.acce)
    if size(labels.acce{i},1)>199 && max(diff(labels.acce{i}(:,1)))<.2
        data=labels.acce{i};
        
        % remove duplicates
        data=sortrows(data);
        data(diff(data(:,1))==0,:)=[];
        
        acce{i}=[t' spline(data(:,1)'-data(1,1),data(:,2:end)',t')'];
    else acce{i}=[];
    end
    if size(labels.gyro{i},1)>199 && max(diff(labels.gyro{i}(:,1)))<.2
        data=labels.gyro{i};
        
        % remove duplicates
        data=sortrows(data);
        data(diff(data(:,1))==0,:)=[];
        
        gyro{i}=[t' spline(data(:,1)'-data(1,1),data(:,2:end)',t')'];
    else gyro{i}=[];
    end
    if size(labels.baro{i},1)>9 && range(labels.baro{i}(:,1))>4
        data=labels.baro{i};
        
        % remove duplicates
        data=sortrows(data);
        data(diff(data(:,1))==0,:)=[];
        
        baro{i}=[tbar' spline(data(:,1)'-data(1,1),data(:,2:end)',tbar')'];
    else baro{i}=[];
    end
end