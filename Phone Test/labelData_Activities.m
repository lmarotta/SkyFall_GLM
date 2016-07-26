clear all
Subj='Nick_Test';
win=5;
max_clips=100000;
%% breaks each set of activity data into 5s clips

acc=load('accTest.mat');
gyr=load('gyrTest.mat');
bar=load('barTest.mat');

newClips=0; % counter on added clips
acce=[acc.NormalizedTimestamp acc.X acc.Y acc.Z];
if isempty(acce)
    continue
end
baro=[bar.tNormalizedTimestamp bar.Altitude bar.Pressure];
if isempty(baro)
    continue
end
gyro=[gyr.NormalizedTimestamp gyr.X gyr.Y gyr.Z];
if isempty(gyro)
    continue
end
tStart=max([min(acce(:,1)),min(baro(:,1)),min(gyro(:,1))]);
tEnd=min([max(acce(:,1)),max(baro(:,1)),max(gyro(:,1))]);

duration=min([tEnd-tStart 5*max_clips]);
numClips=min(floor(duration/win));
for j=1:numClips
    a_ind = acce(:,1)>=tStart+win*(j-1) & acce(:,1)<tStart+win*j;
    b_ind = baro(:,1)>=tStart+win*(j-1) & baro(:,1)<tStart+win*j;
    g_ind = gyro(:,1)>=tStart+win*(j-1) & gyro(:,1)<tStart+win*j;
    new_acce{1+newClips}=acce(a_ind,:);
    new_baro{1+newClips}=baro(b_ind,:);
    new_gyro{1+newClips}=gyro(g_ind,:);
    if numClips>1
        newClips=newClips+1;
    end
end

labels.acce=new_acce;
labels.baro=new_baro;
labels.gyro=new_gyro;

labels.value=repmat(9,[1,length(labels.acce)]);
labels.subject=repmat(Subj,[1,length(labels.acce)]);
labels.text=repmat('Non_Fall_Act',[1,length(labels.acce)]);
labels.timestamp=repmat(0,[1,length(labels.acce)]);

save labels_plus_data_ACT  labels
c=0;

%for j=1:labelsNum
%    if ~isempty(labels.gyro{j})& ~isempty(labels.acce{j})&~isempty(labels.baro{j})
%    c=c+1;
%    end;
%end;