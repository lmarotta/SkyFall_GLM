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
    acce=labels.acce{i};
    gyro=labels.gyro{i};
    baro=labels.baro{i};
    
    acc_counts(i)=size(labels.acce{i},1);
    gyr_counts(i)=size(labels.gyro{i},1);
    bar_counts(i)=size(labels.baro{i},1);
    
%     if size(acc,1)>=100 && size(gyr,1)>=100 && size(bar,1)>=10
    if size(acce,1)>=2 && size(gyro,1)>=2 && size(baro,1)>=2
        acc_gap(i)=max(diff(sort(acce(:,1))));
        gyr_gap(i)=max(diff(sort(gyro(:,1))));
        bar_gap(i)=max(diff(sort(baro(:,1))));
        
        acc_fs=[acc_fs; diff(sort(acce(:,1)))];
        gyr_fs=[gyr_fs; diff(sort(gyro(:,1)))];
        bar_fs=[bar_fs; diff(sort(baro(:,1)))];
    end
end

z_inds=find(acc_gap==0 & gyr_gap==0 & bar_gap==0);
acc_gap(z_inds)=[];
gyr_gap(z_inds)=[];
bar_gap(z_inds)=[];