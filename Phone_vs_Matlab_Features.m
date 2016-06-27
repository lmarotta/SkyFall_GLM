load FallProbe_TestData.mat
load class_params_ACT_old.mat
FN=zeros(1,sum(nz_ind));

for i=1:length(labels.acce)
    
    temp=extract_feature_test_phone(labels.acce{i}, labels.gyro{i}, labels.baro{i}, fvar);
    if ~isempty(temp)
        FN(i,:)=temp(nz_ind);
    end
    
end
