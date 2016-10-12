load FallProbe_TestData
% load class_params_ACT_nobar.mat

% labels = interpolateData(labels);

javaaddpath('purple-robot-skyfall.jar');

save labels_struct_test labels

com.company.ModelFeatureExtractor.extractFeatures('labels_struct_test.mat', 'class_params_ACT_nobar.mat', 'java_features.csv');

FN = csvread('java_features.csv');
FN=FN(:,1:end-1); %delete last column with all zeros

labels.FN=FN;

save PhoneProbe_data_plus_features labels