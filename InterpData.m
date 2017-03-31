%% InterpData
function [acce, gyro, baro]=InterpData(labels)

acce=cell(length(labels.acce),1);
gyro=cell(length(labels.acce),1);
baro=cell(length(labels.acce),1);


t=0:.02:4.98;
tbar=0:1/6:29/6;

% [TDFilt(1,:),TDFilt(2,:)]=butter(2,25/25);

for i=1:length(labels.acce)
    if size(labels.acce{i},1)>199 && max(diff(labels.acce{i}(:,1)))<.2
        data=labels.acce{i};
        
        % remove duplicates
        data=sortrows(data);
        data(diff(data(:,1))==0,:)=[];
        
        t_ind=min(floor(max(data(:,1)-data(1,1))*50),250);
        
        acce{i}=[t(1:t_ind)' spline(data(:,1)'-data(1,1),data(:,2:end)',t(1:t_ind)')'];
%         acce{i}=filtfilt(TDFilt(1,:),TDFilt(2,:),acce{i});
    else acce{i}=[];
    end
    if size(labels.gyro{i},1)>199 && max(diff(labels.gyro{i}(:,1)))<.2
        data=labels.gyro{i};
        
        % remove duplicates
        data=sortrows(data);
        data(diff(data(:,1))==0,:)=[];
        
        t_ind=min(floor(max(data(:,1)-data(1,1))*50),250);
        
        gyro{i}=[t(1:t_ind)' spline(data(:,1)'-data(1,1),data(:,2:end)',t(1:t_ind)')'];
%         gyro{i}=filtfilt(TDFilt(1,:),TDFilt(2,:),gyro{i});
    else gyro{i}=[];
    end
    if size(labels.baro{i},1)>9 && range(labels.baro{i}(:,1))>4
        data=labels.baro{i};
        
        % remove duplicates
        data=sortrows(data);
        data(diff(data(:,1))==0,:)=[];
        
        t_ind=min(floor(max(data(:,1)-data(1,1))*6),30);
        
        baro{i}=[tbar(1:t_ind)' spline(data(:,1)'-data(1,1),data(:,2:end)',tbar(1:t_ind)')'];
    else baro{i}=[];
    end
end