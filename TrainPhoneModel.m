%Current Labels
%slip=1
%trip=2;
%rightfall=3;
%leftfall=4;
%Activities=9; (everything else)

function TrainPhoneModel(saveName,SensCutoff)

close all

features_used = zeros(18,1); features_used(8) = 1; %only magnitude features
% features_used([ 1 2 5 8 9]) = 1; %expanded set

featureset=getFeatureInds(features_used);
% rng(10001)
%% Loading the data

if ~exist('Test_Data_Amputees.mat','file')
    X_Amp = TrainingDataSetup([], [], 10, 1,'Test_Data_Amputees');
else
    X_Amp = load('Test_Data_Amputees');
    X_Amp = X_Amp.F;
end

if ~exist('HealthyData.mat','file')
    X = TrainingDataSetup([], [], 10, 0,'HealthyData');
else
    X = load('HealthyData');
    X = X.F;
end

if ~exist('OutdoorFalls.mat','file')
    O = TrainingDataSetup([], [], 10, 2, 'OutdoorFalls');
else
    O = load('OutdoorFalls.mat');
    O = O.F;
end

F=[X;O;X_Amp]; % dataset with all subjects (Healthy Indoor/Outdoor and Amputees)
%Feature matrix assembled as follows: F = [subj_id location subjcode labels Features];

%convert labels to binary (1-4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = F(:,4); %labels for classification
LB=(L<9);%+1;  %binary labeling
L=LB; %binary labels for fall detection

subjid = F(:,1:3);  %[subjid, location subjcode]
F = F(:,5:end); %the feature matrix

F = F(:,featureset);
subj=unique(subjid(:,1));

%% LOSO CV
%This CV is used to find an optimal threshold to achieve desired
%Sensitivity

%initialize fvar structure
f = cell(length(subj),1);
fvar= struct('std',f,'mu',f,'eps',f,'nzstd',f);

for indCV=1:length(subj)

    test_subj=subj(indCV);
    indtrain = subjid(:,1)~=test_subj;
    indtest = ~indtrain;
    
    % Train GLMnet
    %default values - no grid search over params
    alpha = 0.6;
    lambda = 0.015;

    [fvar_si,b,nz_ind]=Modeltrain(F(indtrain,:),L(indtrain),alpha,lambda);
    fvar(indCV) = fvar_si;
    b_all{indCV}=b;
    nz_ind_all{indCV}=nz_ind;
 
    % Testing the model on Test Data (left out subject)
    [pred,conf,confmat] = Modeleval(F(indtest,:),L(indtest),fvar(indCV),nz_ind,b,.5,0);
    conf_all{indCV}=conf;
    confmat_all(:,:,indCV)=confmat;
    isfall = logical(L(indtest));
    isfall_all{indCV}=isfall;

    % Save Optimal threshold for each subject
    if length(unique(isfall)) >= 2
        [~, TPR, Thresh]=perfcurve(isfall, conf, true);
        OptThres(indCV)=Thresh(find(TPR>=SensCutoff,1));
    else
        OptThres(indCV)=nan;
    end    
    
end

%% Plot ROC curve for all data
[X, Y, T]=perfcurve(isfall_all, conf_all, true,'TVals',[0:0.025:1]);
figure; errorbar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));

isfall=cell2mat(isfall_all');
conf=cell2mat(conf_all');
[X, Y, T, AUC]=perfcurve(isfall, conf, true);
figure; plot(X,Y)
AUC

%% Train on all data and save

[fvar,b,nz_ind]=Modeltrain(F,L,alpha,lambda);
fvar.nzstd = featureset';
Thres=nanmean(OptThres);

save(['.\PhoneModels\' saveName],'fvar','b','nz_ind','Thres')