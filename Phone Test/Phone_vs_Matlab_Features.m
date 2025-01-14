load FallProbe_TestData
load class_params_ACT_nobar.mat

FN=zeros(1,sum(nz_ind));

for i=1:length(labels.acce)
    
    temp=extract_feature_test_phone(labels.acce{i}, labels.gyro{i}, labels.baro{i}, fvar);
    if ~isempty(temp)
        FN(i,:)=temp(nz_ind);
    end
    
end

labels.FN=FN;

save PhoneProbe_data_plus_features labels