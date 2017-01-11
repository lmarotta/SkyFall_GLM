%Generate results for the paper

% LOSO CV on Healthy - Train/Test on Indoor falls
% Input: Features set (0:reduced, 1:full)
         %Number of locations (1:waist, 2:waist+pocket, 3:all)
% Output: AUC, Sens, Spec
function GeneratePaperResults

%% Store Amputees data (all 3 locations) for testing
%Feature matrix assembled as follows
% F = [subj_id location subjcode labels Features];
X = TrainingDataSetup([], [], 10, 1,'Test_Data_Amputees');

[F_amp,L_amp,subjid_amp] = dataextraction(X); %AMPUTEE DATA

%LOSO CV on Healthy 1 location (waist)
X = TrainingDataSetup([], [], 10, 0,'HealthyData');
[AUC,Sens,Spec] = LOSOCV(X,3,1)



end 

function [F,L,subjid] = dataextraction(X)

%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = X(:,4); %labels for classification
L=(L<9);%+1;  %binary labeling
subjid = X(:,1:3);  %[subjid, location subjcode]
F = X(:,5:end); %the feature matrix

end

function [AUC,Sens,Spec] = LOSOCV(X,n_locations,featureset)

%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = X(:,4); %labels for classification
L=(L<9);%+1;  %binary labeling
subjid = X(:,1:3);  %[subjid, location subjcode]
F = X(:,5:end); %the feature matrix
subj=unique(subjid(:,1));

L_w = L(subjid(:,2)<=n_locations);
ratio = sum(L_w)/sum(~L_w);   %falls/non_falls
sprintf('Nlocs = %d, falls/non_falls = %f',n_locations,ratio)

for indCV=1:length(subj)
    
    
    test_subj=subj(indCV);
    indtrain = subjid(:,1)~=test_subj & subjid(:,2)<=n_locations;
    indtest = subjid(:,1) == test_subj;
    
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
    [~, ~, ~, AUC(indCV)]=perfcurve(isfall, conf, true);
    Sens(indCV) = sum(pred & isfall)/sum(isfall);
    Spec(indCV) = sum(~pred & ~isfall)/sum(~isfall);
end

[X, Y, T]=perfcurve(isfall_all(2:end), conf_all(2:end), true,'XVals',[0:0.05:1]); %skip subject 1 (CK) data


end
