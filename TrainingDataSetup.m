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

act_data.acce=data.acce(~fall_inds);
act_data.gyro=data.gyro(~fall_inds);
act_data.baro=data.baro(~fall_inds);
act_data.value=data.value(~fall_inds);
act_data.subject=data.subject(~fall_inds);
act_data.location=data.location(~fall_inds);
act_data.type_str=data.type_str(~fall_inds);

falls_data=JitterClips(falls_data,n);

data.acce=[falls_data.acce; act_data.acce];
data.gyro=[falls_data.gyro; act_data.gyro];
data.baro=[falls_data.baro; act_data.baro];
data.value=[falls_data.value; act_data.value];
data.subject=[falls_data.subject; act_data.subject];
data.location=[falls_data.location; act_data.location];
data.type_str=[falls_data.type_str; act_data.type_str];

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

%% Feature calculation 

F = extract_feature_matlab(labels);
save([savedfilename '.mat'],'F')
save([savedfilename '_labels.mat'], 'labels')
% save Training_Data F
% save Training_Data_labels labels
