%Current Labels
%slip=1
%trip=2;
%rightfall=3;
%leftfall=4;
%Activities=9; (everything else)

% function [confmat_all, conf_all, isfall_all, b_all, fvar, nz_ind_all, FP]=Main_GLM_ACT_Matlab_ThrTuning()

FP=[]; %to store FP and train model stage 2
amputee_test = 1;   %use amputee data only for test
epsF = 1e-6; %threshold on std dev for standardizing features
class=0; % flag for fall classification (rather than detection only)
no_baro=0; % 0 - use barometer
nbaro_features = 106;
% rng(10001)
%% Loading the data

load Training_Data
%Feature matrix assembled as follows
% F = [subj_id location subjcode labels Features];

%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = F(:,4); %labels for classification
LB=(L<9);%+1;  %binary labeling
L=LB; %binary labels for fall detection

if ~amputee_test
    subjid = F(:,1:3);  %[subjid, location subjcode]
    F = F(:,5:end); %the feature matrix
else
    ind = F(:,3) == 0;      %index of amputee subjects
    FAmputee = F(ind,5:end);    %features from amputee only (test)
    subjid = F(~ind,1:3);  %[subjid, location subjcode]
    F = F(~ind,5:end);
    LAmputee = L(ind);      %labels for amputees
    L = L(~ind);            %labels for all others
end

subj=unique(subjid(:,1));
folds_nr = length(subj)-1; %for glmnet

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

    [fvar_si,b,nz_ind]=Modeltrain(F(indtrain,:),L(indtrain),alpha,lambda,no_baro);
    fvar(indCV) = fvar_si;
    b_all{indCV}=b;
    nz_ind_all{indCV}=nz_ind;

 
    % Testing the model on Test Data (left out subject)
    [pred,conf,confmat] = Modeleval(F(indtest,:),L(indtest),fvar(indCV),nz_ind,b);
    conf_all{indCV}=conf;
    confmat_all(:,:,indCV)=confmat;
    isfall = logical(L(indtest));
    isfall_all{indCV}=isfall;

    % Save Optimal threshold for each subject
    if length(unique(isfall)) >= 2
        [~, TPR, Thresh]=perfcurve(isfall, conf, true);
        OptThres(indCV)=Thresh(find(TPR>=0.99,1));
    else
        OptThres(indCV)=nan;
    end    
    
end

fvar_all=fvar;
%% Plot ROC curve for all data
[X, Y, T]=perfcurve(isfall_all(2:end), conf_all(2:end), true,'XVals',[0:0.05:1]); %skip subject 1 (CK) data
figure; errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));

isfall=cell2mat(isfall_all');
conf=cell2mat(conf_all');
[X, Y, T, AUC]=perfcurve(isfall, conf, true);
figure; plot(X,Y)
AUC

%% Train a model on all Healthy lab data
if amputee_test
    Thres = nanmean(OptThres);  %the optimal threshold
    [fvar,b,nz_ind]=Modeltrain(F,L,alpha,lambda,no_baro);
    disp('Training model on all healthy data')
    
    %Evaluate model on AMPUTEES
    [pred,conf,confmat] = Modeleval(FAmputee,LAmputee,fvar,nz_ind,b,Thres,1);
    %Save model trained on healthy and threshold
    save('HealthyModel.mat', 'fvar', 'b', 'nz_ind', 'Thres')
    isfall = LAmputee;
    [X, Y, T, AUC]=perfcurve(isfall, conf, true);
    figure; plot(X,Y,'LineWidth',3)
    xlabel('False positive rate')
    ylabel('True positive rate')
    title('ROC Fall detection in Amputee')

end