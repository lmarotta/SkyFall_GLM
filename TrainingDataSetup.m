%% TrainingDataSetup
% Inputs:
% locations     Cell array of desired data locations {'pouch', 'pocket', 'hand'}
% subjs         Cell array of desired subject IDs
% condition     0: healthy, 1: Amputee, 2: healthy outdoor 
% leave inputs empty ([]) to use all

function F = TrainingDataSetup(locations, subjs, n, condition, savedfilename)

if nargin<5
    savedfilename='Training_Data';
end
if nargin < 4
    condition = 0;
end
if nargin<2
    subjs=[]; n =1;
    savedfilename='Training_Data';
end
if nargin<1
    locations=[]; n=1;
    savedfilename='Training_Data';

end    

switch condition
    case 0 
        filepath = './TrainingData/Healthy/';
    case 1
        filepath = './TrainingData/';
    case 2
        filepath = './TrainingData/Outdoors/';
        f = dir('./TrainingData/Healthy/*.mat');
        S=[];
        for i=1:length(f)
            tempData=load(['./TrainingData/Healthy/' f(i).name]);
            S=[S; unique(tempData.data.subject)];
            S=unique(S);
        end
end

filenames=dir([filepath,'falls_data_10sec_*.mat']);

for indF=1:length(filenames)
    d=load([filepath filenames(indF).name]);
    
    if indF==1, data=d.data;
    else
    
        data.acce=[data.acce; d.data.acce];
        data.gyro=[data.gyro; d.data.gyro];
        data.baro=[data.baro; d.data.baro];

        data.value=[data.value; d.data.value];
        data.subject=[data.subject; d.data.subject];
        data.location=[data.location; d.data.location];
        data.type_str=[data.type_str; d.data.type_str];
    end
    
end

%% Remove undesired subjects and locations
if ~isempty(locations)
    loc_inds=cellfun(@(x) any(strcmp(x,locations)),data.location);
else
    loc_inds=ones(length(data.location),1);
end
if ~isempty(subjs)
    subj_inds=cellfun(@(x) any(strcmp(x,subjs)),data.subject);
else
    subj_inds=ones(length(data.subject),1);
end

empty_inds=cellfun(@isempty, data.acce) | cellfun(@isempty, data.gyro) | cellfun(@isempty, data.baro);

empty_inds=empty_inds | ~loc_inds | ~subj_inds;

data.acce(empty_inds)=[];
data.gyro(empty_inds)=[];
data.baro(empty_inds)=[];
data.value(empty_inds)=[];
data.subject(empty_inds)=[];
data.location(empty_inds)=[];
data.type_str(empty_inds)=[];

%% separate falls and non-falls
fall_inds=data.value~=9;
falls_data.acce=data.acce(fall_inds);
falls_data.gyro=data.gyro(fall_inds);
falls_data.baro=data.baro(fall_inds);
falls_data.value=data.value(fall_inds);
falls_data.subject=data.subject(fall_inds);
falls_data.location=data.location(fall_inds);
falls_data.type_str=data.type_str(fall_inds);

%distribution of max acc
acc = cellfun(@(x) x(:,2:end)./9.81, falls_data.acce,'UniformOutput',false);
maxacc = sqrt(cellfun(@(x) max(sum(x.^2,2)),acc)); %in [g]
figure, hold on, subplot(121)
boxplot(maxacc), ylim([0 4])


act_data.acce=data.acce(~fall_inds);
act_data.gyro=data.gyro(~fall_inds);
act_data.baro=data.baro(~fall_inds);
act_data.value=data.value(~fall_inds);
act_data.subject=data.subject(~fall_inds);
act_data.location=data.location(~fall_inds);
act_data.type_str=data.type_str(~fall_inds);

%distribution of max acc
acc = cellfun(@(x) x(:,2:end)./9.81, act_data.acce,'UniformOutput',false);
maxacc = sqrt(cellfun(@(x) max(sum(x.^2,2)),acc)); %in [g]
subplot(122)
boxplot(maxacc), ylim([0 4])


falls_data=JitterClips(falls_data,n);

%distribution of max acc following jitter
acc = cellfun(@(x) x(:,2:end)./9.81, falls_data.acce,'UniformOutput',false);
maxacc = sqrt(cellfun(@(x) max(sum(x.^2,2)),acc)); %in [g]
figure, hold on, subplot(121)
boxplot(maxacc), ylim([0 4])

%jittering is taking data which are not falls - TO FIX!
indm = find(maxacc < 2);
if ~isempty(indm)
    figure, plot(acc{indm(1)})
else
    display('No mistaken falls')
end

%% remove Falls/Activities under 2g and combine data

fallInds=cellfun(@(x) max(sum(x(:,2:end).^2,2))>(2*9.81)^2,falls_data.acce);
actInds=cellfun(@(x) max(sum(x(:,2:end).^2,2))>(2*9.81)^2,act_data.acce);

data.acce=[falls_data.acce(fallInds); act_data.acce(actInds)];
data.gyro=[falls_data.gyro(fallInds); act_data.gyro(actInds)];
data.baro=[falls_data.baro(fallInds); act_data.baro(actInds)];
data.value=[falls_data.value(fallInds); act_data.value(actInds)];
data.subject=[falls_data.subject(fallInds); act_data.subject(actInds)];
data.location=[falls_data.location(fallInds); act_data.location(actInds)];
data.type_str=[falls_data.type_str(fallInds); act_data.type_str(actInds)];

%distribution of max acc
figure; hold on

acc = cellfun(@(x) x(:,2:end)./9.81, falls_data.acce,'UniformOutput',false);
maxacc = sqrt(cellfun(@(x) max(sum(x.^2,2)),acc)); %in [g]
subplot(121)
boxplot(maxacc), ylim([0 4])

acc = cellfun(@(x) x(:,2:end)./9.81, act_data.acce(actInds),'UniformOutput',false);
maxacc = sqrt(cellfun(@(x) max(sum(x.^2,2)),acc)); %in [g]
subplot(122)
boxplot(maxacc), ylim([0 4])

%% Interpolate Data

[acce, gyro, baro]=InterpData(data);

data.acce=acce;
data.gyro=gyro;
data.baro=baro;

empty_inds=cellfun(@isempty, data.acce) | cellfun(@isempty, data.gyro) | cellfun(@isempty, data.baro);
data.acce(empty_inds)=[];
data.gyro(empty_inds)=[];
data.baro(empty_inds)=[];
data.value(empty_inds)=[];
data.subject(empty_inds)=[];
data.location(empty_inds)=[];
data.type_str(empty_inds)=[];

labels=data;

%% Display histogram of max acceleration norm in each clip 
fall_inds=data.value~=9;
acc = cellfun(@(x) x(:,2:end)./9.81, data.acce,'UniformOutput',false);
maxacc = sqrt(cellfun(@(x) max(sum(x.^2,2)),acc)); %in [g]
figure, subplot(121)
histogram(maxacc)

%% Feature calculation 

F = extract_feature_matlab(labels);

if condition==2
    subjid=cellfun(@(x) iff(length(x)>5,find(strcmp([x(1:3) x(5:end)],S)),find(strcmp(x,S))), labels.subject);
    F(:,1)=subjid;
end

save([savedfilename '.mat'],'F')
save([savedfilename '_labels.mat'], 'labels')
% save Training_Data F
% save Training_Data_labels labels

end
function result=iff(condition, trueResult, falseResult)
    if condition
        result=trueResult;
    else
        result=falseResult;
    end
end
