%Generate results for the paper

% LOSO CV on Healthy - Train/Test on Indoor falls
% Input: Features set (0:reduced, 1:full)
%Number of locations (1:waist, 2:waist+pocket, 3:all)
% Output: AUC, Sens, Spec
function results = GeneratePaperResults
close all



nData=160; %number of data points for training on 1 location only

% features_used = ones(18,1); %full feature set
features_used = zeros(18,1); features_used(8) = 1; %only magnitude features

featureset=getFeatureInds(features_used);

%% Store Amputees data (all 3 locations) for testing
%Feature matrix assembled as follows
% F = [subj_id location subjcode labels Features];
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

% [AUC{2},Sens{2},Spec{2}] = LOSOCV(X,X_Amp,1,0,0); %RF
% [AUC{4},Sens{4},Spec{4}] = LOSOCV(X,X_Amp,2,0,0);
% [AUC{6},Sens{6},Spec{6}] = LOSOCV(X,X_Amp,3,0,0);
% [AUC{1},Sens{1},Spec{1}] = LOSOCV(X,X_Amp,1,1,0); %GLMnet
% [AUC{3},Sens{3},Spec{3}] = LOSOCV(X,X_Amp,2,1,0);


%Train and Test on all 3 locations
cvtype = [1 2 3]; %all cv
% cvtype = 2; %H-A only
[AUC,Sens,Spec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,1:3,1:3,100,1,featureset,cvtype,0);
results.AUC = AUC;
results.AUCErr = AUCErr;
results.Sens = Sens;
results.Spec = Spec;
results.SpecCI = SpecCI;
results.mAUC = cellfun(@nanmean,AUC,'UniformOutput',false);
results.sAUC = cellfun(@nanstd,AUC,'UniformOutput',false);
results.mSens = cellfun(@nanmean,Sens,'UniformOutput',false);
results.sSens = cellfun(@nanstd,Sens,'UniformOutput',false);
results.mSpec = cellfun(@nanmean,Spec,'UniformOutput',false);
results.sSpec = cellfun(@nanstd,Spec,'UniformOutput',false);

%% Train on 1 location and test on 3
cvtype = 2; %H-A only

% Train on waist
[wAUC,wSens,wSpec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,1,1:3,nData,1,featureset,cvtype,0);
results.waist.AUC = wAUC;
results.waist.AUCErr = AUCErr;
results.waist.Sens = wSens;
results.waist.Spec = wSpec;
results.waist.SpecCI = SpecCI;
results.waist.mAUC = cellfun(@nanmean,wAUC,'UniformOutput',false);
results.waist.sAUC = cellfun(@nanstd,wAUC,'UniformOutput',false);
results.waist.mSens = cellfun(@nanmean,wSens,'UniformOutput',false);
results.waist.sSens = cellfun(@nanstd,wSens,'UniformOutput',false);
results.waist.mSpec = cellfun(@nanmean,wSpec,'UniformOutput',false);
results.waist.sSpec = cellfun(@nanstd,wSpec,'UniformOutput',false);

% Train on pocket
[pAUC,pSens,pSpec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,2,1:3,nData,1,featureset,cvtype,0);
results.pock.AUC = pAUC;
results.pock.AUCErr = AUCErr;
results.pock.Sens = pSens;
results.pock.Spec = pSpec;
results.pock.SpecCI = SpecCI;
results.pock.mAUC = cellfun(@nanmean,pAUC,'UniformOutput',false);
results.pock.sAUC = cellfun(@nanstd,pAUC,'UniformOutput',false);
results.pock.mSens = cellfun(@nanmean,pSens,'UniformOutput',false);
results.pock.sSens = cellfun(@nanstd,pSens,'UniformOutput',false);
results.pock.mSpec = cellfun(@nanmean,pSpec,'UniformOutput',false);
results.pock.sSpec = cellfun(@nanstd,pSpec,'UniformOutput',false);

% Train on hand
[hAUC,hSens,hSpec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,3,1:3,nData,1,featureset,cvtype,0);
results.hand.AUC = hAUC;
results.hand.AUCErr = AUCErr;
results.hand.Sens = hSens;
results.hand.Spec = hSpec;
results.hand.SpecCI = SpecCI;
results.hand.mAUC = cellfun(@nanmean,hAUC,'UniformOutput',false);
results.hand.sAUC = cellfun(@nanstd,hAUC,'UniformOutput',false);
results.hand.mSens = cellfun(@nanmean,hSens,'UniformOutput',false);
results.hand.sSens = cellfun(@nanstd,hSens,'UniformOutput',false);
results.hand.mSpec = cellfun(@nanmean,hSpec,'UniformOutput',false);
results.hand.sSpec = cellfun(@nanstd,hSpec,'UniformOutput',false);

%% Plot location results (Healthy-Amputee)
%need to add error bars
figure, hold on
bar(1:4,[results.waist.mAUC{cvtype} results.pock.mAUC{cvtype} results.hand.mAUC{cvtype} results.mAUC{cvtype}])
figauc = errorbar(1:4,[results.waist.mAUC{cvtype} results.pock.mAUC{cvtype} results.hand.mAUC{cvtype} results.mAUC{cvtype}],...
    [results.waist.mAUC{cvtype}-results.waist.AUCErr{cvtype}(1) results.pock.mAUC{cvtype}-results.pock.AUCErr{cvtype}(1) results.hand.mAUC{cvtype}-results.hand.AUCErr{cvtype}(1) results.mAUC{cvtype}-results.AUCErr{cvtype}(1)],...
    [results.waist.AUCErr{cvtype}(2)-results.waist.mAUC{cvtype} results.pock.AUCErr{cvtype}(2)-results.pock.mAUC{cvtype} results.hand.AUCErr{cvtype}(2)-results.hand.mAUC{cvtype} results.AUCErr{cvtype}(2)-results.mAUC{cvtype}],...
    'linewidth',1.5,'linestyle','none','color','k');
h = gca;
h.YLim = [0.4 1];
title('mean AUC')
%plot Sens-Spec
figure, hold on
figSS = bar([results.waist.mSens{cvtype} results.pock.mSens{cvtype} results.hand.mSens{cvtype} results.mSens{cvtype}; ...
    results.waist.mSpec{cvtype} results.pock.mSpec{cvtype} results.hand.mSpec{cvtype} results.mSpec{cvtype}]');
h = gca;
h.YLim = [0.6 1];
title('mean Sens and Spec')

%plot Spec at 0.9 Sens
figure, hold on
bar(1:4,[results.waist.mSpec{cvtype} results.pock.mSpec{cvtype} results.hand.mSpec{cvtype} results.mSpec{cvtype}])
figSSCI = errorbar(1:4,[results.waist.mSpec{cvtype} results.pock.mSpec{cvtype} results.hand.mSpec{cvtype} results.mSpec{cvtype}],...
    [results.waist.mSpec{cvtype}-results.waist.SpecCI{cvtype}(1) results.pock.mSpec{cvtype}-results.pock.SpecCI{cvtype}(1) results.hand.mSpec{cvtype}-results.hand.SpecCI{cvtype}(1) results.mSpec{cvtype}-results.SpecCI{cvtype}(1)],...
    [results.waist.SpecCI{cvtype}(2)-results.waist.mSpec{cvtype} results.pock.SpecCI{cvtype}(2)-results.pock.mSpec{cvtype} results.hand.SpecCI{cvtype}(2)-results.hand.mSpec{cvtype} results.SpecCI{cvtype}(2)-results.mSpec{cvtype}],...
    'linewidth',1.5,'linestyle','none','color','k');
h = gca;
h.YLim = [0.4 1];
title('mean Spec at 90% Sens')



%% HOME DATA ANALYSIS 
cvtype = 2;

filespath = 'Z:/Amputee Phones-R01/Home Data Collection/Amputees/';
l = load([filespath 'HomeDataAmp.mat']);
F = l.F;
sprintf('Data length = %.2f h',size(F,1)*5/60/60)

inds= X_Amp(:,1)==5 | X_Amp(:,1)==6;

X_Amp = [X_Amp(inds,:);F];

% Train on waist
[wAUC,wSens,wSpec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,1,0:3,nData,1,featureset,cvtype,0);
results.waist.AUC = wAUC;
results.waist.AUCErr = AUCErr;
results.waist.Sens = wSens;
results.waist.Spec = wSpec;
results.waist.Spec = SpecCI;
results.waist.mAUC = cellfun(@nanmean,wAUC,'UniformOutput',false);
results.waist.sAUC = cellfun(@nanstd,wAUC,'UniformOutput',false);
results.waist.mSens = cellfun(@nanmean,wSens,'UniformOutput',false);
results.waist.sSens = cellfun(@nanstd,wSens,'UniformOutput',false);
results.waist.mSpec = cellfun(@nanmean,wSpec,'UniformOutput',false);
results.waist.sSpec = cellfun(@nanstd,wSpec,'UniformOutput',false);

% Train on pocket
[pAUC,pSens,pSpec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,2,1:3,nData,1,featureset,cvtype,0);
results.pock.AUC = pAUC;
results.pock.AUCErr = AUCErr;
results.pock.Sens = pSens;
results.pock.Spec = pSpec;
results.pock.Spec = SpecCI;
results.pock.mAUC = cellfun(@nanmean,pAUC,'UniformOutput',false);
results.pock.sAUC = cellfun(@nanstd,pAUC,'UniformOutput',false);
results.pock.mSens = cellfun(@nanmean,pSens,'UniformOutput',false);
results.pock.sSens = cellfun(@nanstd,pSens,'UniformOutput',false);
results.pock.mSpec = cellfun(@nanmean,pSpec,'UniformOutput',false);
results.pock.sSpec = cellfun(@nanstd,pSpec,'UniformOutput',false);

% Train on hand
[hAUC,hSens,hSpec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,3,1:3,nData,1,featureset,cvtype,0);
results.hand.AUC = hAUC;
results.hand.AUCErr = AUCErr;
results.hand.Sens = hSens;
results.hand.Spec = hSpec;
results.hand.Spec = SpecCI;
results.hand.mAUC = cellfun(@nanmean,hAUC,'UniformOutput',false);
results.hand.sAUC = cellfun(@nanstd,hAUC,'UniformOutput',false);
results.hand.mSens = cellfun(@nanmean,hSens,'UniformOutput',false);
results.hand.sSens = cellfun(@nanstd,hSens,'UniformOutput',false);
results.hand.mSpec = cellfun(@nanmean,hSpec,'UniformOutput',false);
results.hand.sSpec = cellfun(@nanstd,hSpec,'UniformOutput',false);

% 3 Locations
[AUC,Sens,Spec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,1:3,1:3,nData,1,featureset,cvtype,0);
results.AUC = AUC;
results.AUCErr = AUCErr;
results.Sens = Sens;
results.Spec = Spec;
results.Spec = SpecCI;
results.mAUC = cellfun(@nanmean,AUC,'UniformOutput',false);
results.sAUC = cellfun(@nanstd,AUC,'UniformOutput',false);
results.mSens = cellfun(@nanmean,Sens,'UniformOutput',false);
results.sSens = cellfun(@nanstd,Sens,'UniformOutput',false);
results.mSpec = cellfun(@nanmean,Spec,'UniformOutput',false);
results.sSpec = cellfun(@nanstd,Spec,'UniformOutput',false);

%% Plot location results (Healthy-Amputee)
%need to add error bars
figure, hold on
bar(1:4,[results.waist.mAUC{cvtype} results.pock.mAUC{cvtype} results.hand.mAUC{cvtype} results.mAUC{cvtype}])
figauc = errorbar(1:4,[results.waist.mAUC{cvtype} results.pock.mAUC{cvtype} results.hand.mAUC{cvtype} results.mAUC{cvtype}],...
    [results.waist.mAUC{cvtype}-results.waist.AUCErr{cvtype}(1) results.pock.mAUC{cvtype}-results.pock.AUCErr{cvtype}(1) results.hand.mAUC{cvtype}-results.hand.AUCErr{cvtype}(1) results.mAUC{cvtype}-results.AUCErr{cvtype}(1)],...
    [results.waist.AUCErr{cvtype}(2)-results.waist.mAUC{cvtype} results.pock.AUCErr{cvtype}(2)-results.pock.mAUC{cvtype} results.hand.AUCErr{cvtype}(2)-results.hand.mAUC{cvtype} results.AUCErr{cvtype}(2)-results.mAUC{cvtype}],...
    'linewidth',1.5,'linestyle','none','color','k');h = gca;
h.YLim = [0.8 1];
title('mean AUC')
%plot Sens-Spec
figure, hold on
figSS = bar([results.waist.mSens{cvtype} results.pock.mSens{cvtype} results.hand.mSens{cvtype} results.mSens{cvtype}; ...
    results.waist.mSpec{cvtype} results.pock.mSpec{cvtype} results.hand.mSpec{cvtype} results.mSpec{cvtype}]');
h = gca;
h.YLim = [0.6 1];
title('mean Sens and Spec')
%plot Spec at 0.9 Sens
figure, hold on
bar(1:4,[results.waist.mSpec{cvtype} results.pock.mSpec{cvtype} results.hand.mSpec{cvtype} results.mSpec{cvtype}])
figSSCI = errorbar(1:4,[results.waist.mSpec{cvtype} results.pock.mSpec{cvtype} results.hand.mSpec{cvtype} results.mSpec{cvtype}],...
    [results.waist.mSpec{cvtype}-results.waist.SpecCI{cvtype}(1) results.pock.mSpec{cvtype}-results.pock.SpecCI{cvtype}(1) results.hand.mSpec{cvtype}-results.hand.SpecCI{cvtype}(1) results.mSpec{cvtype}-results.SpecCI{cvtype}(1)],...
    [results.waist.SpecCI{cvtype}(2)-results.waist.mSpec{cvtype} results.pock.SpecCI{cvtype}(2)-results.pock.mSpec{cvtype} results.hand.SpecCI{cvtype}(2)-results.hand.mSpec{cvtype} results.SpecCI{cvtype}(2)-results.mSpec{cvtype}],...
    'linewidth',1.5,'linestyle','none','color','k');
h = gca;
h.YLim = [0.4 1];
title('mean Spec at 90% Sens')


%% Home Retraining

% Train on waist
[wAUC,wSens,wSpec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,1,0:3,nData,1,featureset,cvtype,1);
results.waist.AUC = wAUC;
results.wiast.AUCErr = AUCErr;
results.waist.Sens = wSens;
results.waist.Spec = wSpec;
results.waist.Spec = SpecCI;
results.waist.mAUC = cellfun(@nanmean,wAUC,'UniformOutput',false);
results.waist.sAUC = cellfun(@nanstd,wAUC,'UniformOutput',false);
results.waist.mSens = cellfun(@nanmean,wSens,'UniformOutput',false);
results.waist.sSens = cellfun(@nanstd,wSens,'UniformOutput',false);
results.waist.mSpec = cellfun(@nanmean,wSpec,'UniformOutput',false);
results.waist.sSpec = cellfun(@nanstd,wSpec,'UniformOutput',false);

% Train on pocket
[pAUC,pSens,pSpec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,2,1:3,nData,1,featureset,cvtype,1);
results.pock.AUC = pAUC;
results.pock.AUCErr = AUCErr;
results.pock.Sens = pSens;
results.pock.Spec = pSpec;
results.pock.Spec = SpecCI;
results.pock.mAUC = cellfun(@nanmean,pAUC,'UniformOutput',false);
results.pock.sAUC = cellfun(@nanstd,pAUC,'UniformOutput',false);
results.pock.mSens = cellfun(@nanmean,pSens,'UniformOutput',false);
results.pock.sSens = cellfun(@nanstd,pSens,'UniformOutput',false);
results.pock.mSpec = cellfun(@nanmean,pSpec,'UniformOutput',false);
results.pock.sSpec = cellfun(@nanstd,pSpec,'UniformOutput',false);

% Train on hand
[hAUC,hSens,hSpec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,3,1:3,nData,1,featureset,cvtype,1);
results.hand.AUC = hAUC;
results.hand.AUCErr = AUCErr;
results.hand.Sens = hSens;
results.hand.Spec = hSpec;
results.hand.Spec = SpecCI;
results.hand.mAUC = cellfun(@nanmean,hAUC,'UniformOutput',false);
results.hand.sAUC = cellfun(@nanstd,hAUC,'UniformOutput',false);
results.hand.mSens = cellfun(@nanmean,hSens,'UniformOutput',false);
results.hand.sSens = cellfun(@nanstd,hSens,'UniformOutput',false);
results.hand.mSpec = cellfun(@nanmean,hSpec,'UniformOutput',false);
results.hand.sSpec = cellfun(@nanstd,hSpec,'UniformOutput',false);

% 3 Locations
[AUC,Sens,Spec,AUCErr,SpecCI] = LOSOCV(X,X_Amp,1:3,1:3,nData,1,featureset,cvtype,1);
results.AUC = AUC;
results.AUCErr = AUCErr;
results.Sens = Sens;
results.Spec = Spec;
results.Spec = SpecCI;
results.mAUC = cellfun(@nanmean,AUC,'UniformOutput',false);
results.sAUC = cellfun(@nanstd,AUC,'UniformOutput',false);
results.mSens = cellfun(@nanmean,Sens,'UniformOutput',false);
results.sSens = cellfun(@nanstd,Sens,'UniformOutput',false);
results.mSpec = cellfun(@nanmean,Spec,'UniformOutput',false);
results.sSpec = cellfun(@nanstd,Spec,'UniformOutput',false);

%% Plot location results (Healthy-Amputee)
%need to add error bars
figure, hold on
bar(1:4,[results.waist.mAUC{cvtype} results.pock.mAUC{cvtype} results.hand.mAUC{cvtype} results.mAUC{cvtype}])
figauc = errorbar(1:4,[results.waist.mAUC{cvtype} results.pock.mAUC{cvtype} results.hand.mAUC{cvtype} results.mAUC{cvtype}],...
    [results.waist.mAUC{cvtype}-results.waist.AUCErr{cvtype}(1) results.pock.mAUC{cvtype}-results.pock.AUCErr{cvtype}(1) results.hand.mAUC{cvtype}-results.hand.AUCErr{cvtype}(1) results.mAUC{cvtype}-results.AUCErr{cvtype}(1)],...
    [results.waist.AUCErr{cvtype}(2)-results.waist.mAUC{cvtype} results.pock.AUCErr{cvtype}(2)-results.pock.mAUC{cvtype} results.hand.AUCErr{cvtype}(2)-results.hand.mAUC{cvtype} results.AUCErr{cvtype}(2)-results.mAUC{cvtype}],...
    'linewidth',1.5,'linestyle','none','color','k');h = gca;
h.YLim = [0.8 1];
title('mean AUC')
%plot Sens-Spec
figure, hold on
figSS = bar([results.waist.mSens{cvtype} results.pock.mSens{cvtype} results.hand.mSens{cvtype} results.mSens{cvtype}; ...
    results.waist.mSpec{cvtype} results.pock.mSpec{cvtype} results.hand.mSpec{cvtype} results.mSpec{cvtype}]');
h = gca;
h.YLim = [0.6 1];
title('mean Sens and Spec')
%plot Spec at 0.9 Sens
figure, hold on
bar(1:4,[results.waist.mSpec{cvtype} results.pock.mSpec{cvtype} results.hand.mSpec{cvtype} results.mSpec{cvtype}])
figSSCI = errorbar(1:4,[results.waist.mSpec{cvtype} results.pock.mSpec{cvtype} results.hand.mSpec{cvtype} results.mSpec{cvtype}],...
    [results.waist.mSpec{cvtype}-results.waist.SpecCI{cvtype}(1) results.pock.mSpec{cvtype}-results.pock.SpecCI{cvtype}(1) results.hand.mSpec{cvtype}-results.hand.SpecCI{cvtype}(1) results.mSpec{cvtype}-results.SpecCI{cvtype}(1)],...
    [results.waist.SpecCI{cvtype}(2)-results.waist.mSpec{cvtype} results.pock.SpecCI{cvtype}(2)-results.pock.mSpec{cvtype} results.hand.SpecCI{cvtype}(2)-results.hand.mSpec{cvtype} results.SpecCI{cvtype}(2)-results.mSpec{cvtype}],...
    'linewidth',1.5,'linestyle','none','color','k');
h = gca;
h.YLim = [0.4 1];
title('mean Spec at 90% Sens')

end



%cvtype is a vector with the cases
%[1 2 3] = H-H, H-A, A-A
function [AUC,Sens,Spec,AUCErr,SpecCI] = LOSOCV(X,X_test,locations_train,locations_test,nData,model,featureset,cvtype,HomeFP_retrain)

rng(400)

if HomeFP_retrain
    l=load('Z:/Amputee Phones-R01/Home Data Collection/Healthy/HomeDataHealthy.mat');
    HomeData=l.F;
    HomeData = HomeData(:,featureset);
    HomeL = zeros(size(HomeData,1),1); %all labels are non-fall
end

%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = X(:,4); %labels for classification
L=(L<9);%+1;  %binary labeling
subjid = X(:,1:3);  %[subjid, location subjcode]
subj=unique(subjid(:,1));
F=X(:,5:end); %remove metadata
F = F(:,featureset);

L_w = L(any(bsxfun(@eq,subjid(:,2),locations_train),2));
ratio = sum(L_w)/sum(~L_w);   %falls/non_falls
sprintf('Nlocs = %d, falls/non_falls = %f',length(locations_train),ratio)

conf_all=cell(1,length(subj));
isfall_all=cell(1,length(subj));
confmat_all=zeros(2,2,length(subj));
Nruns = 5;  %average multiple runs for each subject
alpha = 0.6;
lambda = 0.015;

%% LOSO CV HEALTHY
AUC_HH=zeros(Nruns,length(subj));
Sens_HH=zeros(Nruns,length(subj));
Spec_HH=zeros(Nruns,length(subj));

if find(cvtype==1)
    
    for indCV=1:length(subj)
        
        
        test_subj=subj(indCV);
        indtrain = subjid(:,1)~=test_subj & any(bsxfun(@eq,subjid(:,2),locations_train),2);
        indtest = subjid(:,1) == test_subj & any(bsxfun(@eq,subjid(:,2),locations_test),2);
        
        falls_ind = find(indtrain & L);
        nfalls_ind = find(indtrain & ~L);
        
        %balance the dataset - randomly sample nData instances from each class
        %repeat it 5 times and average results
        for run = 1:Nruns
            indtrain = [falls_ind(randperm(length(falls_ind),nData)); nfalls_ind(randperm(length(nfalls_ind),nData))];
            
            if model
                % Train GLMnet
                %default values - no grid search over params
                alpha = 0.6;
                lambda = 0.015;
                
                [fvar_si,b,nz_ind]=Modeltrain(F(indtrain,:),L(indtrain),alpha,lambda,0);
                % Testing the model on Test Data (left out subject)
                [pred,conf,confmat] = Modeleval(F(indtest,:),L(indtest),fvar_si,nz_ind,b,0.5,0);
                conf_all{indCV}=conf;
                confmat_all(:,:,indCV)=confmat;
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                [~, ~, ~, AUC_HH(run,indCV)]=perfcurve(isfall, conf, true);
                Sens_HH(run,indCV) = sum(pred & isfall)/sum(isfall);
                Spec_HH(run,indCV) = sum(~pred & ~isfall)/sum(~isfall);
            else
                RFModel=TreeBagger(100,F(indtrain,:),L(indtrain));
                [pred,conf]=predict(RFModel,F(indtest,:));
                conf = conf(:,2); %posterior prob of a fall
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                [~, ~, ~, AUC_HH(run,indCV)]=perfcurve(isfall, conf, true);
                
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                pred=cellfun(@(x) logical(str2double(x)),pred);
                
                confmat=confusionmat(isfall,pred);
                confmat_all(:,:,indCV)=confmat;
                
                Sens_HH(run,indCV) = sum(pred & isfall)/sum(isfall);
                Spec_HH(run,indCV) = sum(~pred & ~isfall)/sum(~isfall);
            end
        end
    end
    %average over runs
    AUC_HH = nanmean(AUC_HH,1);
    Sens_HH = nanmean(Sens_HH,1);
    Spec_HH = nanmean(Spec_HH,1);
    
    confmat=sum(confmat_all,3);
    plotConfmat(confmat,'Healthy-Healthy');
    
    %plot ROC curves w confidence bounds
    figroc = figure;
    [X, Y, T, AUC]=perfcurve(isfall_all, conf_all, true,'XVals',[0:0.05:1]); %conf bounds with CV
%     [X, Y, T, AUC]=perfcurve(cell2mat(isfall_all'), cell2mat(conf_all'), true,'Nboot',1000,'XVals',[0:0.05:1]); %cb with Bootstrap
    e = errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
    e.LineWidth = 2; e.Marker = 'o';
    xlabel('False positive rate')
    ylabel('True positive rate')
    
end
%% Train on Healthy, Test on Amputees
AUC_HA=NaN(Nruns,length(unique(X_test(:,1))),'double');
Sens_HA=NaN(Nruns,length(unique(X_test(:,1))),'double');
Spec_HA=NaN(Nruns,length(unique(X_test(:,1))),'double');
if any(cvtype == 2)
    
    isfall_all = {};
    conf_all = {};
    Ft = X_test(:,5:end);
    Ft = Ft(:,featureset);
    
    indtrain = any(bsxfun(@eq,subjid(:,2),locations_train),2);
    
    falls_ind = find(indtrain & L);
    nfalls_ind = find(indtrain & ~L);
    
    for run = 1:Nruns
        indtrain = [falls_ind(randperm(length(falls_ind),nData)); nfalls_ind(randperm(length(nfalls_ind),nData))];
        
        if model
            [fvar_si,b,nz_ind]=Modeltrain(F(indtrain,:),L(indtrain),alpha,lambda,0);
            
            if HomeFP_retrain
                [pred,~,~] = Modeleval(HomeData,HomeL,fvar_si,nz_ind,b,.5,0);
                sprintf('Spec = %.2f%', length(pred)/(length(pred)+sum(pred)))
                indmisc=find(pred); %misclassified clips
                sprintf('# of misclassified clips = %d/%d',length(indmisc),length(F))
                Fmisc = HomeData(indmisc,:); Lmisc = false(length(indmisc),1);
                FPHome = [F(indtrain,:);Fmisc]; LPHome = [L(indtrain);Lmisc];
                disp('Training model on Healthy lab + misclassified home data')
                [fvar_si,b,nz_ind]=Modeltrain(FPHome,LPHome,alpha,lambda,0);
            end
            
            % Testing the model on each amputee subj (external set)
            S=unique(X_test(:,1));
            for i = 1:length(unique(X_test(:,1)))
                subj=S(i);
                rowid = (X_test(:,1)==subj) & any(bsxfun(@eq,X_test(:,2),locations_test),2);
                [pred,conf,confmat] = Modeleval(Ft(rowid,:),X_test(rowid,4)<5,fvar_si,nz_ind,b,0.5,0);
                conf_all{i}=conf;
                confmat_all(:,:,subj)=confmat;
                isfall = X_test(rowid,4)<5;
                isfall_all{i}=isfall;
                if length(unique(isfall)) <2 %subject AF06 with missing activities - TO FIX
                    isfall_all{i} = [];  conf_all{i}=[];
                else
                    [~, ~, ~, AUC_HA(run,i)]=perfcurve(isfall, conf, true);
                    Sens_HA(run,i) = sum(pred & isfall)/sum(isfall);
                    Spec_HA(run,i) = sum(~pred & ~isfall)/sum(~isfall);
                end
            end
        else
            RFModel=TreeBagger(100,F(indtrain,:),L(indtrain));
            for subj = 1:length(unique(X_test(:,1)))
                rowid = (X_test(:,1)==subj) & any(bsxfun(@eq,X_test(:,2),locations_test),2);
                [pred,conf]=predict(RFModel,Ft(rowid,:));
                conf = conf(:,2); %posterior prob of a fall
                conf_all{subj}=conf;
                isfall = X_test(rowid,4)<5;
                isfall_all{subj}=isfall;
                [~, ~, ~, AUC_HA(run,subj)]=perfcurve(isfall, conf, true);
                Sens_HA(run,subj) = sum(pred & isfall)/sum(isfall);
                Spec_HA(run,subj) = sum(~pred & ~isfall)/sum(~isfall);
                confmat=confusionmat(X_test(:,4)<9,cellfun(@(x) logical(str2double(x)),pred));
                confmat_all(:,:,subj)=confmat;
                
            end
        end
    end
    AUC_HA = nanmean(AUC_HA,1);
    Sens_HA = nanmean(Sens_HA,1);
    Spec_HA = nanmean(Spec_HA,1);
    
    confmat=sum(confmat_all,3);
    plotConfmat(confmat,'Healthy-Amputee');
    if ~exist('figroc','var')
        figroc = figure;
    end
    figure(figroc), hold on
    [X, Y, T, AUC]=perfcurve(isfall_all(~cellfun(@isempty,isfall_all)), conf_all(~cellfun(@isempty,isfall_all)), true,'XVals',[0:0.05:1]);
%     [X, Y, T, AUC]=perfcurve(cell2mat(isfall_all'), cell2mat(conf_all'), true,'Nboot',1000,'XVals',[0:0.05:1]); %cb with Bootstrap
    e = errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1)); %plot ROC curve

    %95% CI on specificity for fixed sensitivity
    sens = 0.9;
    TL = [cell2mat(conf_all') cell2mat(cellfun(@double,isfall_all','UniformOutput',false))];
    bootstat = bootstrp(1000,@specAUC,TL,sens);
    ci = bootci(1000,@specAUC,TL,sens);
    ci_spec = ci(:,1); ci_AUC = ci(:,2);
    
    Spec_HA = mean(bootstat(:,1)); %the spec at 90% sensitivity
    AUC_HA = mean(bootstat(:,2));  %mean AUC
    
    e.LineWidth = 2; e.Marker = 'o';
    AUCerr_HA=ci_AUC;
    SpecCI_HA = ci_spec;

    

end

%% LOSO on Amputees

AUC_AA=NaN(Nruns,length(subj),'double');
Sens_AA=NaN(Nruns,length(subj),'double');
Spec_AA=NaN(Nruns,length(subj),'double');
if any(cvtype == 3)
    
    
    isfall_all = {};
    conf_all = {};
    
    F = X_test(:,5:end);
    F = F(:,featureset);
    L = X_test(:,4); %labels for classification
    L=(L<9);%+1;  %binary labeling
    subjid = X_test(:,1:3);  %[subjid, location subjcode]
    subj=unique(subjid(:,1));
    
    for indCV=1:length(subj)
        
        test_subj=subj(indCV);
        indtrain = subjid(:,1)~=test_subj & any(bsxfun(@eq,subjid(:,2),locations_train),2);
        indtest = subjid(:,1) == test_subj & any(bsxfun(@eq,subjid(:,2),locations_test),2);
        
        falls_ind = find(indtrain & L);
        nfalls_ind = find(indtrain & ~L);
        
        for run = 1:Nruns
            %balance the dataset
            indtrain = [falls_ind(randperm(length(falls_ind),nData)); nfalls_ind(randperm(length(nfalls_ind),nData))];
            
            if model
                % Train GLMnet
                %default values - no grid search over params
                alpha = 0.6;
                lambda = 0.015;
                
                [fvar_si,b,nz_ind]=Modeltrain(F(indtrain,:),L(indtrain),alpha,lambda,0);
                % Testing the model on Test Data (left out subject)
                [pred,conf,confmat] = Modeleval(F(indtest,:),L(indtest),fvar_si,nz_ind,b,0.5,0);
                conf_all{indCV}=conf;
                confmat_all(:,:,indCV)=confmat;
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                if length(unique(isfall)) <2 %subject AF06 with missing activities - TO FIX
                    isfall_all{indCV} = [];  conf_all{indCV}=[];
                else
                    [~, ~, ~, AUC_AA(run,indCV)]=perfcurve(isfall, conf, true);
                    Sens_AA(run,indCV) = sum(pred & isfall)/sum(isfall);
                    Spec_AA(run,indCV) = sum(~pred & ~isfall)/sum(~isfall);
                end
            else
                RFModel=TreeBagger(100,F(indtrain,:),L(indtrain));
                [pred,conf]=predict(RFModel,F(indtest,:));
                conf = conf(:,2); %posterior prob of a fall
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                [~, ~, ~, AUC_AA(run,indCV)]=perfcurve(isfall, conf, true);
                
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                pred=cellfun(@(x) logical(str2double(x)),pred);
                
                confmat=confusionmat(isfall,pred);
                confmat_all(:,:,indCV)=confmat;
                
                Sens_AA(run,indCV) = sum(pred & isfall)/sum(isfall);
                Spec_AA(run,indCV) = sum(~pred & ~isfall)/sum(~isfall);
            end
        end
    end
    AUC_AA = nanmean(AUC_AA,1);
    Sens_AA = nanmean(Sens_AA,1);
    Spec_AA = nanmean(Spec_AA,1);
    
    confmat=sum(confmat_all,3);
    plotConfmat(confmat,'Amputee-Amputee');
    if ~exist('figroc','var')
        figroc = figure;
    end
    figure(figroc), hold on
    [X, Y, T, AUC]=perfcurve(isfall_all(~cellfun(@isempty,isfall_all)), conf_all(~cellfun(@isempty,isfall_all)), true,'XVals',[0:0.05:1]);
%     [X, Y, T, AUC]=perfcurve(cell2mat(isfall_all'),cell2mat(conf_all'),true,'Nboot',1000,'XVals',[0:0.05:1]);
    e = errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
    e.LineWidth = 2; e.Marker = 'o';
    legend('Healthy-Healthy','Healthy-Amputee','Amputee-Amputee')
end


AUC = {AUC_HH;AUC_HA;AUC_AA};
Sens = {Sens_HH;Sens_HA;Sens_AA};
Spec = {Spec_HH;Spec_HA;Spec_AA};
AUCErr = {[];AUCerr_HA;[]};
SpecCI = {[];SpecCI_HA;[]};


end

function plotConfmat(confmat,figtitle)

activities = {'Non-Fall','Fall'};
figure; imagesc(confmat./repmat(sum(confmat,2),[1 2]));
confmat = confmat./repmat(sum(confmat,2),[1 2]);
% [cmin,cmax] = caxis;
caxis([0,1]) %set colormap to 0 1
ax = gca;
ax.XTick = 1:size(activities,2);
ax.YTick = 1:size(activities,2);
set(gca,'XTickLabel',activities,'FontSize',14)
set(gca,'YTickLabel',activities,'FontSize',14)
ax.XTickLabelRotation = 45;
axis square
if isempty(figtitle)
    title('Fall Detection')
else
    title(figtitle);
end
%add text
thres = 0.6;
for i = 1:length(activities)
    for j = 1:length(activities)
        if confmat(i,j) < thres
            col = [1 1 1];
        else
            col = [0 0 0];
        end
        text(i-0.2,j,sprintf('%.3f',confmat(j,i)),'FontSize',12,'FontWeight','bold','Color',col);
    end
end
end
