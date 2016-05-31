
%%%%% Reading Gyroscope and Accelerometer
load labelsData.mat
win=5;

labelsNum= numel(labels.value);

fl=find(labels.value<=8);
lv=repmat(fl',1,3);
lvsz= length(fl);
lind= true(lvsz,3);
files_list= dir;
fsz= numel(files_list);


for i=1:fsz
    fn=files_list(i).name;
    flag= regexp(fn,'Sensor');
    if ~isempty(flag)
        load(fn);
        switch Data_type
            case 'accelerometer'
                cind=lv(lind(:,1),1);
                for j=1:numel(cind)
                    ctstamp=labels.timestamp(cind(j));
                    id= NormalizedTimestamp>ctstamp-win/4& NormalizedTimestamp<ctstamp+3*win/4;
                    if any(id==true)
                        labels.acce{cind(j)}= [NormalizedTimestamp(id,:),X(id,:),Y(id,:),Z(id,:)];
                        indf = find(lv(:,1)==cind(j));
                        lind(indf,1)=false;
                    end;
                end;
            case 'Gyro'
                 cind=lv(lind(:,2),2);
                for j=1:numel(cind)
                    ctstamp=labels.timestamp(cind(j));
                    id= NormalizedTimestamp>ctstamp-win/4& NormalizedTimestamp<ctstamp+3*win/4;
                    if any(id==true)
                        labels.gyro{cind(j)}= [NormalizedTimestamp(id,:),X(id,:),Y(id,:),Z(id,:)];
                        indf = find(lv(:,2)==cind(j));
                        lind(indf,2)=false;
                    end;
                end;
                
            case 'barometer'
                cind=lv(lind(:,3),3);
                for j=1:numel(cind)
                    ctstamp=labels.timestamp(cind(j));
                    id= tNormalizedTimestamp>ctstamp-win/4& tNormalizedTimestamp<ctstamp+3*win/4;
                    if any(id==true)
                        labels.baro{cind(j)}= [tNormalizedTimestamp(id,:),Altitude(id,:),Pressure(id,:)];
                        indf = find(lv(:,3)==cind(j));
                        lind(indf,3)=false;
                    end;
              end;
        end;
    end;
end;

save labels_plus_data  labels
c=0;

%for j=1:labelsNum
%    if ~isempty(labels.gyro{j})& ~isempty(labels.acce{j})&~isempty(labels.baro{j})
%    c=c+1;
%    end;
%end;