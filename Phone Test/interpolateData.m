%% InterpData
function data = interpolateData(labels)

% load FallProbe_TestData

t=0:.02:4.98;
% tbar=0:1/6:29/6;

for i=1:length(labels.acce)
    if size(labels.acce{i},1)>99 && range(labels.acce{i}(:,1))>4
        currData=labels.acce{i};
        
        % remove duplicates
        currData=sortrows(currData);
        currData(diff(currData(:,1))==0,:)=[];
        
        labels.acce{i}=[t' spline(currData(:,1)'-currData(1,1),currData(:,2:end)',t')'];
    else labels.acce{i}=[];
    end
    if size(labels.gyro{i},1)>99 && range(labels.gyro{i}(:,1))>4
        currData=labels.gyro{i};
        
        % remove duplicates
        currData=sortrows(currData);
        currData(diff(currData(:,1))==0,:)=[];
        
        labels.gyro{i}=[t' spline(currData(:,1)'-currData(1,1),currData(:,2:end)',t')'];
    else labels.gyro{i}=[];
    end
    %{
    if size(labels.baro{i},1)>9 && range(labels.baro{i}(:,1))>4
        data=labels.baro{i};
        
        % remove duplicates
        data=sortrows(data);
        data(diff(data(:,1))==0,:)=[];
        
        labels.baro{i}=[tbar' spline(data(:,1)'-data(1,1),data(:,2:end)',tbar')'];
    else labels.baro{i}=[];
    end
    %}
end

data = labels;
% save FallProbe_TestData_interp.mat labels

end