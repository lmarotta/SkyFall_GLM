 
%%%%% Reading Gyroscope and Accelerometer
clear all
load FallProbe_TestData.mat
win=5;

files_list= dir;
fsz= numel(files_list);

clip_sizes = zeros(numel(labels.timestamp),3);


for i=1:fsz
    fn=files_list(i).name;
    flag= regexp(fn,'Test.mat');
    if ~isempty(flag)
        load(fn);
        switch Data_type
            case 'accelerometer'
                for j=1:numel(labels.timestamp)
                    clipStart = labels.timestampSTART_END(j,1);
                    clipEnd = labels.timestampSTART_END(j,2);
                    id= NormalizedTimestamp>clipStart & NormalizedTimestamp<clipEnd;
                    if any(id==true)
                        labels.acce{j}= [NormalizedTimestamp(id,:),X(id,:),Y(id,:),Z(id,:)];
                        clip_sizes(j,1) = size(labels.acce{j},1);
                    else
                        labels.acce{j}=[];
                    end;               
                end;
            case 'Gyro'
                for j=1:numel(labels.timestamp)
                    clipStart = labels.timestampSTART_END(j,1);
                    clipEnd = labels.timestampSTART_END(j,2);
                    id= NormalizedTimestamp>clipStart & NormalizedTimestamp<clipEnd;
                    if any(id==true)
                        labels.gyro{j}= [NormalizedTimestamp(id,:),X(id,:),Y(id,:),Z(id,:)];
                        clip_sizes(j,2) = size(labels.gyro{j},1);
                    else
                        labels.gyro{j}=[];
                    end;
                end;
                
            case 'barometer'
                for j=1:numel(labels.timestamp)
                    clipStart = labels.timestampSTART_END(j,1);
                    clipEnd = labels.timestampSTART_END(j,2);
                    id= tNormalizedTimestamp>clipStart & tNormalizedTimestamp<clipEnd;
                    if any(id==true)
                        labels.baro{j}= [tNormalizedTimestamp(id,:),Altitude(id,:),Pressure(id,:)];
                        clip_sizes(j,3) = size(labels.baro{j},1);
                    else
                        labels.baro{j}=[];
                    end;
                end;
        end;
    end;
end;

labels.clipSizes = clip_sizes;

save PhoneProbe_data  labels

t = 1:size(labels.clipSizes,1);
figure, plot(t,labels.clipSizes(:,1), t,labels.clipSizes(:,2), t,labels.clipSizes(:,3)), legend('Accelerometer','Gyroscope','Barometer')
xlabel('clip'), ylabel('# of samples')
title('# of samples in clips')

%histogram of # of samples in a clip
figure, hold on, subplot(311), histogram(labels.clipSizes(:,1)), title('acc'), xlabel('# of samples in clip')
subplot(312), histogram(labels.clipSizes(:,2)), title('gyr')
subplot(313), histogram(labels.clipSizes(:,3)), title('bar')