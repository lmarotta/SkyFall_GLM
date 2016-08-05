len=length(labels.acce);
acc_gap=zeros(len,1);
gyr_gap=zeros(len,1);
bar_gap=zeros(len,1);

acc_fs=[];
gyr_fs=[];
bar_fs=[];

acc_counts=zeros(len,1);
gyr_counts=zeros(len,1);
bar_counts=zeros(len,1);

for i=1:len
    acc=labels.acce{i};
    gyr=labels.gyro{i};
    bar=labels.baro{i};
    
    acc_counts(i)=size(labels.acce{i},1);
    gyr_counts(i)=size(labels.gyro{i},1);
    bar_counts(i)=size(labels.baro{i},1);
    
%     if size(acc,1)>=100 && size(gyr,1)>=100 && size(bar,1)>=10
    if size(acc,1)>=2 && size(gyr,1)>=2 && size(bar,1)>=2
        acc_gap(i)=max(diff(sort(acc(:,1))));
        gyr_gap(i)=max(diff(sort(gyr(:,1))));
        bar_gap(i)=max(diff(sort(bar(:,1))));
        
        acc_fs=[acc_fs; diff(sort(acc(:,1)))];
        gyr_fs=[gyr_fs; diff(sort(gyr(:,1)))];
        bar_fs=[bar_fs; diff(sort(bar(:,1)))];
    end
end

z_inds=find(acc_gap==0 & gyr_gap==0 & bar_gap==0);
acc_gap(z_inds)=[];
gyr_gap(z_inds)=[];
bar_gap(z_inds)=[];