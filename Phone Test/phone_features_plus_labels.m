load PhoneProbe_data_plus_features
load features_labels_full
load class_params_ACT_nobar

phone_features = labels.features;
phone_features = phone_features(:,2:end);

features_labels = featuresLabels(:,fvar.nzstd);
features_labels = features_labels(:,nz_ind);

features_bad = mean(phone_features,1) > 25;
features_bad_labels = features_labels(:,features_bad);