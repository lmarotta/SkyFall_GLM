load PhoneProbe_data_plus_features;

javaFeatures = csvread('java_features.csv');
javaFeatures = javaFeatures(:,1:end-1); %delete last column with all zeros
phoneFeatures = labels.values;
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

errors = (phoneFeaturesNZ-javaFeaturesNZ)./phoneFeaturesNZ;
figure, bar(mean(errors)), ylim([-1 3])
title('Errors Phone VS Java Features')
xlabel('Features')
ylabel('Error Value')

%%Phone features vs Matlab features 
matlabFeaturesNZ = labels.FN(ind,:);
ratios_phone_vs_matlab = phoneFeaturesNZ.*(matlabFeaturesNZ.^-1);
figure, bar(mean(ratios_phone_vs_matlab))%, ylim([-1 3])
title('Ratios Phone VS Matlab Features')
xlabel('Features')
ylabel('Feature Value')

errors = (phoneFeaturesNZ-matlabFeaturesNZ)./phoneFeaturesNZ;
figure, bar(mean(errors)), ylim([-1 3])
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
figure, bar(mean(errors)), ylim([-1 3])
title('Errors Java VS Matlab Features')
xlabel('Features')
ylabel('Error Value')
