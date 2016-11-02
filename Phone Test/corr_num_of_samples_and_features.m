 load FallProbe_TestData
 load features_labels
 load class_params_ACT_nobar.mat
 
 labels_features_used=featuresLabels(nz_ind);
 
 acc_ind=strmatch('acce',labels_features_used);
 gyr_ind=strmatch('gyro',labels_features_used);
 bar_ind=strmatch('baro',labels_features_used);
 
 corrs=zeros(length(labels_features_used),1);
 
 X=[];
 Y=[];
 
 Y=labels.sensor_counts(:,2);

 for i=gyr_ind
     X=labels.values(:,i+1);
 
     corrs(i)=corr(X,Y);
 end
 
 Y=labels.sensor_counts(:,1);

 for i=acc_ind
     X=labels.values(:,i+1);
 
     corrs(i)=corr(X,Y);
 end