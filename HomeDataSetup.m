%generate features from phone unlabeled (home) data
%INPUT: data structure with raw and interpolated phone data
function F = HomeDataSetup(data)

%% Interpolated Data

datai.acce = data.acce_inerp;
datai.gyro = data.gyro_inerp;
datai.baro = data.baro;

javaaddpath('FeatureGeneration/purple-robot-skyfall.jar');
chunkSize=4000;
F=[];
for indChunk=1:ceil(length(datai.acce)/chunkSize)
    
    startind=(indChunk-1)*chunkSize+1;
    endind=min(indChunk*chunkSize,length(data.acce));
    
    labels.acce=datai.acce(startind:endind);
    labels.gyro=datai.gyro(startind:endind);
    labels.baro=datai.baro(startind:endind);
    
    save labels_temp.mat labels
    
    com.company.TrainFeatureExtractor.extractFeatures('labels_temp.mat', 'chunk_features.csv');
    
    F=[F; xlsread('chunk_features.csv')];
    
end
