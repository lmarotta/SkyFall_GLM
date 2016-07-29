 
%%%%% Reading Gyroscope and Accelerometer
load FallProbe_TestData.mat
win=5;

files_list= dir;
fsz= numel(files_list);


for i=1:fsz
    fn=files_list(i).name;
    flag= regexp(fn,'Test.mat');
    if ~isempty(flag)
        load(fn);
        switch Data_type
            case 'accelerometer'
                for j=1:numel(labels.timestamp)
                    ctstamp=labels.timestamp(j);
                    id= NormalizedTimestamp>ctstamp-win & NormalizedTimestamp<ctstamp;
                    if any(id==true)
                        labels.acce{j}= [NormalizedTimestamp(id,:),X(id,:),Y(id,:),Z(id,:)];
                    else
                        labels.acce{j}=[];
                    end;
                end;
            case 'Gyro'
                for j=1:numel(labels.timestamp)
                    ctstamp=labels.timestamp(j);
                    id= NormalizedTimestamp>ctstamp-win & NormalizedTimestamp<ctstamp;
                    if any(id==true)
                        labels.gyro{j}= [NormalizedTimestamp(id,:),X(id,:),Y(id,:),Z(id,:)];
                    else
                        labels.gyro{j}=[];
                    end;
                end;
                
            case 'barometer'
                for j=1:numel(labels.timestamp)
                    ctstamp=labels.timestamp(j);
                    id= tNormalizedTimestamp>ctstamp-win & tNormalizedTimestamp<ctstamp;
                    if any(id==true)
                        labels.baro{j}= [tNormalizedTimestamp(id,:),Altitude(id,:),Pressure(id,:)];
                    else
                        labels.baro{j}=[];
                    end;
                end;
        end;
    end;
end;

save PhoneProbe_data  labels
c=0;