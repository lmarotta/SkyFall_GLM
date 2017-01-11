%Generate results for the paper

% LOSO CV on Healthy - Train/Test on Indoor falls
% Input: Features set (0:reduced, 1:full)
         %Number of locations (1:waist, 2:waist+pocket, 3:all)
% Output: AUC, Sens, Spec
function GeneratePaperResults

close all
%% Store Amputees data (all 3 locations) for testing
%Feature matrix assembled as follows
% F = [subj_id location subjcode labels Features];
if ~exist('Test_Data_Amputees.mat','file')
    X_Amp = TrainingDataSetup([], [], 10, 1,'Test_Data_Amputees');
else
    X_Amp = load('Test_Data_Amputees');
    X_Amp = X_Amp.F;
end

%LOSO CV on Healthy 1 location (waist)
if ~exist('HealthyData.mat','file')
    X = TrainingDataSetup([], [], 10, 0,'HealthyData');
else
    X = load('HealthyData');
    X = X.F;
end

[AUC{2},Sens{2},Spec{2}] = LOSOCV(X,X_Amp,1,1);
[AUC{4},Sens{4},Spec{4}] = LOSOCV(X,X_Amp,2,1);
[AUC{6},Sens{6},Spec{6}] = LOSOCV(X,X_Amp,3,1);
[AUC{1},Sens{1},Spec{1}] = LOSOCV(X,X_Amp,1,0);
[AUC{3},Sens{3},Spec{3}] = LOSOCV(X,X_Amp,2,0);
[AUC{5},Sens{5},Spec{5}] = LOSOCV(X,X_Amp,3,0);

mAUC=cellfun(@mean,AUC);
sAUC=cellfun(@std,AUC);

figure; hold on; bar(mAUC); for i=1:6; errorbar(i,mAUC(i), 0, sAUC(i),'Color','k'); end

mSens=cellfun(@mean,Sens);
sSens=cellfun(@std,Sens);

figure; hold on; bar(mSens); for i=1:6; errorbar(i,mSens(i), 0, sSens(i),'Color','k'); end

mSpec=cellfun(@mean,Spec);
sSpec=cellfun(@std,Spec);

figure; hold on; bar(mSpec); for i=1:6; errorbar(i,mSpec(i), 0, sSpec(i),'Color','k'); end

end 

function [F,L,subjid] = dataextraction(X)

%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = X(:,4); %labels for classification
L=(L<9);%+1;  %binary labeling
subjid = X(:,1:3);  %[subjid, location subjcode]
F = X(:,5:end); %the feature matrix

end

function [AUC,Sens,Spec] = LOSOCV(X,X_test,n_locations,featureset)

%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = X(:,4); %labels for classification
L=(L<9);%+1;  %binary labeling
subjid = X(:,1:3);  %[subjid, location subjcode]
if featureset
    F = X(:,5:end); %the feature matrix
else
    F = X(:,938:956); % only magnitude features
end
subj=unique(subjid(:,1));

L_w = L(subjid(:,2)<=n_locations);
ratio = sum(L_w)/sum(~L_w);   %falls/non_falls
sprintf('Nlocs = %d, falls/non_falls = %f',n_locations,ratio)

AUC=zeros(1,length(subj));
Sens=zeros(1,length(subj));
Spec=zeros(1,length(subj));

conf_all=cell(1,length(subj));
isfall_all=cell(1,length(subj));
confmat_all=zeros(2,2,length(subj));

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

confmat=sum(confmat_all,3);
plotConfmat(confmat);

if ~isempty(X_test)
    
    if featureset
        Ft = X_test(:,5:end); %the feature matrix
    else
        Ft = X_test(:,938:956); % only magnitude features
    end
    
    indtrain = subjid(:,2)<=n_locations;
    [fvar_si,b,nz_ind]=Modeltrain(F(indtrain,:),L(indtrain),alpha,lambda,0);
    % Testing the model on Test Data (external set)
    [~,conf,confmat] = Modeleval(Ft,X_test(:,4)<5,fvar_si,nz_ind,b,0.5,0);
    
    plotConfmat(confmat)
    [X, Y, ~]=perfcurve(X_test(:,4)<9, conf, true,'XVals',[0:0.05:1]); %skip subject 1 (CK) data
    figure; plot(X,Y);
    
end


[X, Y, ~]=perfcurve(isfall_all, conf_all, true,'XVals',[0:0.05:1]); %skip subject 1 (CK) data
figure; errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));


end

function plotConfmat(confmat)
    
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
    title('Fall Detection')
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
