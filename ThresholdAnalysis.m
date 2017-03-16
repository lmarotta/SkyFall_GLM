%% Threshold Analysis

if ~exist('Test_Data_Amputees.mat','file')
    X_Amp = TrainingDataSetup([], [], 10, 1,'Test_Data_Amputees');
else
    X_Amp = load('Test_Data_Amputees_labels');
    X_Amp = X_Amp.labels;
end

if ~exist('HealthyData.mat','file')
    X = TrainingDataSetup([], [], 10, 0,'HealthyData');
else
    X = load('HealthyData_labels');
    X = X.labels;
end

UFT = min(cellfun(@(x) sqrt(max(sum(x(:,2:end).^2,2))),X.acce(X.value<9)));
LFT = max(cellfun(@(x) sqrt(min(sum(x(:,2:end).^2,2))),X.acce(X.value<9)));

pred = ThresholdDetection(UFT,LFT,X_Amp.acce);
FNR_HA = sum(~pred(X_Amp.value<9))/sum(X_Amp.value<9)
FPR_HA = sum(pred(X_Amp.value==9))/sum(X_Amp.value==9)

filespath = 'Z:/Amputee Phones-R01/Home Data Collection/Amputees/';
f = dir([filespath '/Raw/*.mat']);
for i=1:length(f)
    l = load([filespath '/Raw/' f(i).name]);
    HomeData = l.labels;
    HighAccInds = cellfun(@(x) sqrt(max(sum(x(:,2:end).^2,2)))/9.81>2,HomeData.acce);
    pred = ThresholdDetection(UFT,LFT,HomeData.acce(HighAccInds));
    HomeAcc(i) = sum(~pred)/length(pred);
    Size(i) = sum(HighAccInds);
end

HomeAcc = nansum(HomeAcc.*Size)/sum(Size)

%% Healthy to Healthy CV

subs=unique(X.subject);
for indSub=1:length(subs)
    TrainSet=X.acce(X.value<9 & strcmp(X.subject,subs(indSub)));
    TestSet=X.acce(~strcmp(X.subject,subs(indSub)));
    TestVals=X.value(~strcmp(X.subject,subs(indSub)));
    
    UFT = min(cellfun(@(x) sqrt(max(sum(x(:,2:end).^2,2))),TrainSet));
    LFT = max(cellfun(@(x) sqrt(min(sum(x(:,2:end).^2,2))),TrainSet));

    pred = ThresholdDetection(UFT,LFT,TestSet);
    FNR_HH(indSub) = sum(~pred(TestVals<9))/sum(TestVals<9);
    FPR_HH(indSub) = sum(pred(TestVals==9))/sum(TestVals==9);
    
end