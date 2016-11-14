%Available Datasets
%labels_plus_data - original momo's datasets (fall and fall-like, )
%labels_plus_data_ACT - activities dataset (only contains activities from 3
%amputees)
%labels_full - cleaned dataset with Falls and Activities

%% Analyze Original dataset
clear all
% load labels_plus_data.mat %original momo's dataset
load labels_plus_data_ACT.mat %original momo's dataset

%remove empty clips
ind1 = cellfun(@isempty, labels.acce);
ind2 = cellfun(@isempty, labels.gyro);
ind3 = cellfun(@isempty, labels.baro);
inds = ind1+ind2+ind3;
inds = logical(inds);
labels.timestamp(inds) =[];
labels.value(inds) = [];
labels.subject(inds) = [];
labels.acce(inds) = [];
labels.baro(inds) = [];
labels.gyro(inds) = [];
labels = rmfield(labels,'text');
labels

%% check size of clips
acce_clipsize = cellfun(@length,labels.acce);
gyro_clipsize = cellfun(@length,labels.gyro);
baro_clipsize = cellfun(@length,labels.baro);
%distribution of clip sizes
figure,
subplot(311), histogram(acce_clipsize)
subplot(312), histogram(gyro_clipsize)
subplot(313), histogram(baro_clipsize)

%% remove clips with less than 100 (and more than 1000 samples)
ind1 = acce_clipsize < 100; %| acce_clipsize > 1000;
ind2 = gyro_clipsize < 100; %| gyro_clipsize > 1000;
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

%% show distribution of data across subjects
us = unique(labels.subject);
labels.subjcode = zeros(length(labels.subject),1);
for s = 1:length(us)
    us(s);
    ind = find(strcmp(labels.subject,us(s)));
    labels.subjcode(ind) = s;
end
figure, histogram(labels.subjcode)

%% check size of clips and fix aliasing issue
acce_clipsize = cellfun(@length,labels.acce);
gyro_clipsize = cellfun(@length,labels.gyro);
baro_clipsize = cellfun(@length,labels.baro);

%accelerometer clips
ind = acce_clipsize > 600; %clips with repeating data (logical)
N = length(find(ind));
inds = find(ind); %clip indices
%--plot clips 
find(ind)
acce_clipsize(ind)
acc = labels.acce(ind);
for clip = 1:length(find(ind))
    figure, plot(acc{clip}(:,2:end))
end
%--find periodicity
for clip = 1:N
    [autocor,lags]=xcorr(labels.acce{inds(clip)}(:,2:end),'unbiased');
    autocorXYZ = [autocor(:,1) autocor(:,5) autocor(:,end)];
    autocorXYZ = autocorXYZ(lags > 0,:);
    lags = lags(lags>0);
    autocor = mean(autocorXYZ,2);
    %locate peaks after smoothing
    autocor_s = smooth(autocor);
    autocor_s(autocor_s < 0) = 0;
%     findpeaks(autocor_s); %to plot peaks
    [pks,locs]=findpeaks(autocor_s);
    %find peak between 200 and 300 and take the max
    ind = (locs<300 & locs>200);
    locs = locs(ind); pks = pks(ind);
    [~,indpk] = max(pks);
    endclip = locs(indpk); %the end clip index
    labels.acce{inds(clip)} = labels.acce{inds(clip)}(1:endclip,:);
end

%% gyro clips 
ind = gyro_clipsize > 1000; %clips with repeating data (logical)
N = length(find(ind));
inds = find(ind); %clip indices
%--plot clips 
for clip = 1:length(inds)
    figure, plot(labels.gyro{inds(clip)}(:,2:end))
end
%--find periodicity
for clip = 1:N
    [autocor,lags]=xcorr(labels.gyro{inds(clip)}(:,2:end),'unbiased');
    autocorXYZ = [autocor(:,1) autocor(:,5) autocor(:,end)];
    autocorXYZ = autocorXYZ(lags > 0,:);
    lags = lags(lags>0);
    autocor = mean(autocorXYZ,2);
    %locate peaks after smoothing
    autocor_s = smooth(autocor);
    autocor_s(autocor_s < 0) = 0;
%     findpeaks(autocor_s); %to plot peaks
    [pks,locs]=findpeaks(autocor_s);
    %find peak between 200 and 300 and take the max
    ind = (locs<300 & locs>200);
    locs = locs(ind); pks = pks(ind);
    [~,indpk] = max(pks);
    endclip = locs(indpk); %the end clip index
    labels.gyro{inds(clip)} = labels.gyro{inds(clip)}(1:endclip,:);
end
for clip = 1:length(inds)
    figure, plot(labels.gyro{inds(clip)}(:,2:end))
end


%% check
%distribution of clip sizes
acce_clipsize = cellfun(@length,labels.acce);
gyro_clipsize = cellfun(@length,labels.gyro);
baro_clipsize = cellfun(@length,labels.baro);

figure,
subplot(311), histogram(acce_clipsize)
subplot(312), histogram(gyro_clipsize)
subplot(313), histogram(baro_clipsize)

%% the distribution of fall and fall-like
figure,
histogram(labels.value)

%% Save the data 
save labels_plus_data_Cleaned labels
%% Cleaned Dataset w Falls and Activities
load labels_full.mat

%discard clips with more than 1000 samples and less than 100
acce_clipsize = cellfun(@length,labels.acce);
gyro_clipsize = cellfun(@length,labels.gyro);
baro_clipsize = cellfun(@length,labels.baro);

figure,
subplot(311), histogram(acce_clipsize)
subplot(312), histogram(gyro_clipsize)
subplot(313), histogram(baro_clipsize)

%resample clips 
