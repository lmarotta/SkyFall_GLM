clear all
%%%%% Reading Gyroscope and Accelerometer
load labelsData.mat
% specify clip window size
win=5;

labelsNum= numel(labels.value);

% Only use Activity TimeStamps
fl=find(labels.value>8);
lv=repmat(fl',1,3);
lvsz= length(fl);
lind= true(lvsz,3);

files_list= dir;
fsz= numel(files_list);

endStamps=[];
startStamps=[];
values=[];
subject={};

for i=9:2:17
   
    endInd=find(labels.value==i+1);
    startInd=find(labels.value==i);
    
    endStamps=[endStamps; labels.timestamp(endInd).'];
    startStamps=[startStamps; labels.timestamp(startInd).'];
    values=[values; repmat(i,[length(endInd) 1])];
    subject=[subject labels.subject(endInd)]; 
    
end

labels.value=values;
labels.subject=subject;
startStamps=startStamps([1:231 233:length(startStamps)]);

lind=true(length(startStamps),3);

%%
for i=1:fsz
    fn=files_list(i).name;
    flag= regexp(fn,'Sensor');
    if ~isempty(flag)
        load(fn);
        switch Data_type
            case 'accelerometer'
                curr_ststamp=startStamps(lind(:,1));
                curr_etstamp=endStamps(lind(:,1));                
                for j=1:numel(curr_ststamp)
                    id= NormalizedTimestamp>curr_ststamp(j) & NormalizedTimestamp<curr_etstamp(j);
                    if any(id==true)
                        labels.acce{j}= [NormalizedTimestamp(id,:),X(id,:),Y(id,:),Z(id,:)];
                        lind(j,1)=false;
                    end
                end
            case 'Gyro'
                curr_ststamp=startStamps(lind(:,2));
                curr_etstamp=endStamps(lind(:,2));     
                for j=1:numel(curr_ststamp)
                    id= NormalizedTimestamp>curr_ststamp(j) & NormalizedTimestamp<curr_etstamp(j);
                    if any(id==true)
                        labels.gyro{j}= [NormalizedTimestamp(id,:),X(id,:),Y(id,:),Z(id,:)];
                        lind(j,2)=false;
                    end
                end
                
            case 'barometer'
                curr_ststamp=startStamps(lind(:,3));
                curr_etstamp=endStamps(lind(:,3));     
                for j=1:numel(curr_ststamp)
                    id= tNormalizedTimestamp>curr_ststamp(j) & tNormalizedTimestamp<curr_etstamp(j);
                    if any(id==true)
                        labels.baro{j}= [tNormalizedTimestamp(id,:),Altitude(id,:),Pressure(id,:)];
                        lind(j,3)=false;
                    end
                end
        end
    end
end

%% breaks each set of activity data into 5s clips

newClips=0; % counter on added clips
for i=1:length(labels.acce)
    display(i)
    subject=labels.subject{i};
    acce=labels.acce{i};
    if isempty(acce)
        continue
    end
    baro=labels.baro{i};
    if isempty(baro)
        continue
    end
    gyro=labels.gyro{i};
    if isempty(gyro)
        continue
    end
    tStart=max([min(acce(:,1)),min(baro(:,1)),min(gyro(:,1))]);
    tEnd=min([max(acce(:,1)),max(baro(:,1)),max(gyro(:,1))]);

    duration=min([tEnd-tStart 500]);
    numClips=min(floor(duration/win));
    for j=1:numClips
        a_ind = acce(:,1)>=tStart+win*(j-1) & acce(:,1)<tStart+win*j;
        b_ind = baro(:,1)>=tStart+win*(j-1) & baro(:,1)<tStart+win*j;
        g_ind = gyro(:,1)>=tStart+win*(j-1) & gyro(:,1)<tStart+win*j;
        new_acce{i+newClips}=acce(a_ind,:);
        new_baro{i+newClips}=baro(b_ind,:);
        new_gyro{i+newClips}=gyro(g_ind,:);
        new_subject{i+newClips}=subject;
        if numClips>1
            newClips=newClips+1;
        end
    end
end

labels.acce=new_acce;
labels.baro=new_baro;
labels.gyro=new_gyro;

labels.value=repmat(9,[1,length(labels.acce)]);
labels.subject=new_subject;
labels.text=repmat(labels.text(1),[1,length(labels.acce)]);
labels.timestamp=repmat(labels.timestamp(1),[1,length(labels.acce)]);

save labels_plus_data_ACT  labels
c=0;

%for j=1:labelsNum
%    if ~isempty(labels.gyro{j})& ~isempty(labels.acce{j})&~isempty(labels.baro{j})
%    c=c+1;
%    end;
%end;