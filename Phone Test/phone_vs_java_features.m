load PhoneProbe_data_plus_features;
load features_labels_full;
load class_params_ACT_nobar;

features_labels = featuresLabels(:,fvar.nzstd);
features_labels = features_labels(:,nz_ind);

%javaFeatures = csvread('java_features.csv');
%javaFeatures = javaFeatures(:,1:end-1); %delete last column with all zeros
javaFeatures = labels.FN;
phoneFeatures = labels.features;
phoneFeatures = phoneFeatures(:,2:end); %delete first column with all zeros

%transpose the java matrix to find clips (columns) with all zero indices
ind = any(javaFeatures,2);
javaFeaturesNZ = javaFeatures(ind,:);
phoneFeaturesNZ = phoneFeatures(ind,:);

%Phone features vs Java features
ratios_phone_vs_java = phoneFeaturesNZ.*(javaFeaturesNZ.^-1);
figure, bar(mean(ratios_phone_vs_java))%, ylim([-1 3])
title('Ratios Phone VS Java Features')
xlabel('Features')
ylabel('Feature Value')

badFeaturesRatioInd = mean(ratios_phone_vs_java) > 2 | mean(ratios_phone_vs_java) < -2;
badFeaturesRatio = features_labels(badFeaturesRatioInd);

errors = (phoneFeaturesNZ-javaFeaturesNZ)./phoneFeaturesNZ;
figure, bar(mean(errors))%, ylim([-1 3])
title('Errors Phone VS Java Features')
xlabel('Features')
ylabel('Error Value')

badFeaturesErrorsInd = mean(errors) > 2 | mean(errors) < -2;
badFeaturesErrors = features_labels(badFeaturesErrorsInd);



%%Phone features vs Matlab features
%{
matlabFeaturesNZ = labels.FN(ind,:);
ratios_phone_vs_matlab = phoneFeaturesNZ.*(matlabFeaturesNZ.^-1);
figure, bar(mean(ratios_phone_vs_matlab))%, ylim([-1 3])
title('Ratios Phone VS Matlab Features')
xlabel('Features')
ylabel('Feature Value')

errors = (phoneFeaturesNZ-matlabFeaturesNZ)./phoneFeaturesNZ;
figure, bar(mean(errors))%, ylim([-1 3])
title('Errors Phone VS Matlab Features')
xlabel('Features')
ylabel('Error Value')

%Matlab features vs Java features
Ratios = javaFeaturesNZ.*(matlabFeaturesNZ.^-1);
figure, bar(mean(Ratios))%, ylim([-1 3])
title('Ratios Java VS Matlab Features')
xlabel('Features')
ylabel('Feature Value')

errors = (javaFeaturesNZ-matlabFeaturesNZ)./javaFeaturesNZ;
figure, bar(mean(errors))%, ylim([-1 3])
title('Errors Java VS Matlab Features')
xlabel('Features')
ylabel('Error Value')
%}
