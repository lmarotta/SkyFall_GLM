%Test trained GLMnet model with jittered datapoints
function Test_Jittered

%% Load the classifier parameters  
% load class_params_dec.mat
%  List of variables in class_params:
%   b:  coefficients of logistic regression classifier
%   fvar: a structure containing the necessary vaiables for feature
%   nz_ind: indices of selected features extraction

%the list of jittered clips
files = dir('labels_plus_data_jittered_*');

for f = 1:length(files)
    clear id conf
    filename = files(f).name;
    load(filename), disp(filename)
    
    %sensors data, each k corresponds to 1 data point
    for k = 1:length(labels.acce)
        accel_data = labels.acce{k};
        gyro_data = labels.gyro{k};
        baro_data = labels.baro{k};
        [id(k), conf(k)]= fall_model_eval(accel_data, gyro_data,  baro_data);
    end
    
    isfall = labels.value < 5;
    fall_err(f) = sum(~id(isfall))/sum(isfall);
    nfall_err(f) = sum(id(~isfall))/sum(~isfall);
    err(f) = (sum(id ~= isfall))/length(id); %error rate
end

err = err([1 2 4 3 6 5 7]); %last corresponds to uniform distribution
stdwin = [.01 .1 1 1.5  2 2.5];

%plot acc with unif jittered clips

figure, plot(stdwin,err(1:6),'o-')

figure, plot([stdwin(1) 5], err([1 7]),'o-')

%plot falls acc with unif jittered clips

figure, plot(stdwin,fall_err(1:6),'o-')

figure, plot([stdwin(1) 5], fall_err([1 7]),'o-')

%plot ~falls acc with unif jittered clips

figure, plot(stdwin,nfall_err(1:6),'o-')

figure, plot([stdwin(1) 5], nfall_err([1 7]),'o-')