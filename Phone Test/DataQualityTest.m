%Available Datasets
%labels_plus_data - original momo's datasets (fall and fall-like, )
%labels_plus_data_ACT - activities dataset (only contains activities from 3
%amputees)
%labels_full - cleaned dataset with Falls and Activities

%% Analyze Original dataset
load labels_plus_data.mat %original momo's dataset
% load labels_plus_data_ACT.mat %original momo's dataset

%remove empty clips
ind1 = cellfun(@isempty, labels.acce);
ind2 = cellfun(@isempty, labels.gyro);
ind3 = cellfun(@isempty, labels.baro);
inds = ind1+ind2+ind3;

%remove empty clips
inds = logical(inds);
labels.timestamp(inds) =[];
labels.value(inds) = [];
labels.subject(inds) = [];
labels.acce(inds) = [];
labels.baro(inds) = [];
labels.gyro(inds) = [];
labels = rmfield(labels,'text');
labels

%check size of clips
acce_clipsize = cellfun(@length,labels.acce);
gyro_clipsize = cellfun(@length,labels.gyro);
baro_clipsize = cellfun(@length,labels.baro);

ind1 = acce_clipsize < 100;
ind2 = gyro_clipsize < 100;
ind3 = baro_clipsize < 10;
inds = logical(ind1+ind2+ind3);
labels.timestamp(inds) =[];
labels.value(inds) = [];
labels.subject(inds) = [];
labels.acce(inds) = [];
labels.baro(inds) = [];
labels.gyro(inds) = [];
labels
unique(labels.subject)
unique(labels.value)

%% Cleaned Dataset
load labels_full.mat

%discard clips with more than 1000 samples and less than 100
acce_clipsize = cellfun(@length,labels.acce);
gyro_clipsize = cellfun(@length,labels.gyro);
baro_clipsize = cellfun(@length,labels.baro);

figure,
subplot(311), histogram(acce_clipsize,40)
subplot(312), histogram(gyro_clipsize)
subplot(313), histogram(baro_clipsize)

%resample clips 
