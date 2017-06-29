function F = extract_feature_matlab(labels)

% labelsNum= numel(labels.acce);
% 
% %extract feature dimension
% i = 1;
% data{1}= labels.gyro{i}(:, 2:end);
% data{2}= labels.acce{i}(:, 2:end);
% data{3}= labels.baro{i}(:, 2:end);
% stamp{1}= labels.gyro{i}(:, 1);
% stamp{2}= labels.acce{i}(:, 1);
% stamp{3}= labels.baro{i}(:, 1);
%     
% F1 = CalcFeatures(data,stamp);
% Nfeat = size(F1,2);
% F=zeros(labelsNum,Nfeat); %hard code feature size
% F(1,:) = F1;
% 
% acce = labels.acce;
% gyro = labels.gyro;
% baro = labels.baro;
% %loop through the clips (matlab)
% parfor i=2:labelsNum
%             
%     data1= gyro{i}(:, 2:end);
%     data2= acce{i}(:, 2:end);
%     data3= baro{i}(:, 2:end);
%     stamp1= gyro{i}(:, 1);
%     stamp2= acce{i}(:, 1);
%     stamp3= baro{i}(:, 1);
%     data = {data1 data2 data3};
%     stamp = {stamp1 stamp2 stamp3};
%     
%     F(i,:)=CalcFeatures(data,stamp);
%    
% end

javaaddpath('FeatureGeneration/purple-robot-skyfall.jar');

chunkSize=4000;

data= labels;

F=[];

for indChunk=1:ceil(length(labels.acce)/chunkSize)
    
    startind=(indChunk-1)*chunkSize+1;
    endind=min(indChunk*chunkSize,length(data.acce));
    
    labels.acce=data.acce(startind:endind);
    labels.gyro=data.gyro(startind:endind);
    labels.baro=data.baro(startind:endind);
    
    save labels_temp.mat labels
    
    com.company.TrainFeatureExtractor.extractFeatures('labels_temp.mat', 'chunk_features.csv');
    
    F=[F; csvread('chunk_features.csv')];
    
end

%assign subject code and location to feature matrix
%subject codes are increasing
cvSubj=unique(labels.subject);
subj_id=zeros(length(labels.subject), 1);    
% Assign subject and location folds
for indCVSubj=1:length(cvSubj)
    inds=strcmp(labels.subject,cvSubj{indCVSubj});
    subj_id(inds)=indCVSubj;
end

%Location (0 = NA, 1=Pouch, 2=Pocket, 3=Hand)
locs = unique(labels.location);
location = zeros(size(labels.location,1),1);
for i =1:length(locs)
    inds = strcmp(labels.location,locs(i));
    switch locs{i}
        case 'NA'
            location(inds) = 0;
        case 'pouch'
            location(inds) = 1;
        case 'pocket'
            location(inds) = 2;
        case 'hand'
            location(inds) = 3;
        case 'hand '
            location(inds) = 3;
    end
end    
          

%convert subjectcode to numeric
%1 = CF (healthy) 2 = AF (amputee)
if ~isempty(labels.value)
    subjcode = cellfun(@(x) strcmp(x(1:2),'CF'),labels.subject);
    % subjtype = cellfun(@(x) x(1:2),labels.subject,'UniformOutput',false)
end
    F = [subj_id location subjcode labels.value F];
end