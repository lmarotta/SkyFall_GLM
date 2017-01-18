%Generate results for the paper

% LOSO CV on Healthy - Train/Test on Indoor falls
% Input: Features set (0:reduced, 1:full)
%Number of locations (1:waist, 2:waist+pocket, 3:all)
% Output: AUC, Sens, Spec
function results = GeneratePaperResults

nData=80;

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

% [AUC{2},Sens{2},Spec{2}] = LOSOCV(X,X_Amp,1,0,0); %RF
% [AUC{4},Sens{4},Spec{4}] = LOSOCV(X,X_Amp,2,0,0);
% [AUC{6},Sens{6},Spec{6}] = LOSOCV(X,X_Amp,3,0,0);
% [AUC{1},Sens{1},Spec{1}] = LOSOCV(X,X_Amp,1,1,0); %GLMnet
% [AUC{3},Sens{3},Spec{3}] = LOSOCV(X,X_Amp,2,1,0);


[AUC,Sens,Spec] = LOSOCV(X,X_Amp,1:3,1:3,400,1,0);

% Train on waist
[wAUC{1},wSens{1},wSpec{1}] = LOSOCV(X,X_Amp,1,1,nData,1,0);
[wAUC{2},wSens{2},wSpec{2}] = LOSOCV(X,X_Amp,1,2,nData,1,0);
[wAUC{3},wSens{3},wSpec{3}] = LOSOCV(X,X_Amp,1,3,nData,1,0);

% Train on pocket
[pAUC{1},pSens{1},pSpec{1}] = LOSOCV(X,X_Amp,2,1,nData,1,0);
[pAUC{2},pSens{2},pSpec{2}] = LOSOCV(X,X_Amp,2,2,nData,1,0);
[pAUC{3},pSens{3},pSpec{3}] = LOSOCV(X,X_Amp,2,3,nData,1,0);

% Train on hand
[hAUC{1},hSens{1},hSpec{1}] = LOSOCV(X,X_Amp,3,1,nData,1,0);
[hAUC{2},hSens{2},hSpec{2}] = LOSOCV(X,X_Amp,3,2,nData,1,0);
[hAUC{3},hSens{3},hSpec{3}] = LOSOCV(X,X_Amp,3,3,nData,1,0);

AUC=[cellfun(@(x) cellfun(@nanmean, x),wAUC,'UniformOutput',false); ...
    cellfun(@(x) cellfun(@nanmean, x),pAUC,'UniformOutput',false); ...
    cellfun(@(x) cellfun(@nanmean, x),hAUC,'UniformOutput',false)];

HH_AUC=cellfun(@(x) x(1),AUC);
HA_AUC=cellfun(@(x) x(2),AUC);
AA_AUC=cellfun(@(x) x(3),AUC);

figure; imagesc(HH_AUC); caxis([0 1]);
figure; imagesc(HA_AUC); caxis([0 1]);
figure; imagesc(AA_AUC); caxis([0 1]);

mAUC=cellfun(@mean,AUC);
sAUC=cellfun(@std,AUC);

% figure; hold on; bar(mAUC); for i=1:6; errorbar(i,mAUC(i), 0, sAUC(i),'Color','k'); end
%
% mSens=cellfun(@mean,Sens);
% sSens=cellfun(@std,Sens);
%
% figure; hold on; bar(mSens); for i=1:6; errorbar(i,mSens(i), 0, sSens(i),'Color','k'); end
%
% mSpec=cellfun(@mean,Spec);
% sSpec=cellfun(@std,Spec);
%
% figure; hold on; bar(mSpec); for i=1:6; errorbar(i,mSpec(i), 0, sSpec(i),'Color','k'); end

results.AUC = AUC;
results.Sens = Sens;
results.Spec = Spec;

results.waist.AUC = wAUC;
results.waist.Sens = wSens;
results.waist.Spec = wSpec;

results.pock.AUC = pAUC;
results.pock.Sens = pSens;
results.pock.Spec = pSpec;

results.hand.AUC = hAUC;
results.hand.Sens = hSens;
results.hand.Spec = hSpec;

end

function [F,L,subjid] = dataextraction(X)

%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = X(:,4); %labels for classification
L=(L<9);%+1;  %binary labeling
subjid = X(:,1:3);  %[subjid, location subjcode]
F = X(:,5:end); %the feature matrix

end

function [AUC,Sens,Spec] = LOSOCV(X,X_test,locations_train,locations_test,nData,model,featureset)

%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = X(:,4); %labels for classification
L=(L<9);%+1;  %binary labeling
subjid = X(:,1:3);  %[subjid, location subjcode]
if featureset
    F=X(:,5:end);
else
    F = X(:,943:962); % only magnitude features
end

subj=unique(subjid(:,1));

L_w = L(any(bsxfun(@eq,subjid(:,2),locations_train),2));
ratio = sum(L_w)/sum(~L_w);   %falls/non_falls
sprintf('Nlocs = %d, falls/non_falls = %f',length(locations_train),ratio)


conf_all=cell(1,length(subj));
isfall_all=cell(1,length(subj));
confmat_all=zeros(2,2,length(subj));

%% LOSO CV HEALTHY
Nruns = 1;  %average multiple runs for each subject
AUC_HH=zeros(Nruns,length(subj));
Sens_HH=zeros(Nruns,length(subj));
Spec_HH=zeros(Nruns,length(subj));

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
% [X, Y, T, AUC]=perfcurve(cell2mat(isfall_all'), cell2mat(conf_all'), true,'Nboot',0,'XVals',[0:0.05:1]); %cb with Bootstrap
e = errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
e.LineWidth = 2; e.Marker = 'o';
xlabel('False positive rate')
ylabel('True positive rate')

%% Train on Healthy, Test on Amputees
if ~isempty(X_test)
    AUC_HA=NaN(Nruns,length(unique(X_test(:,1))),'double');
    Sens_HA=NaN(Nruns,length(unique(X_test(:,1))),'double');
    Spec_HA=NaN(Nruns,length(unique(X_test(:,1))),'double');
    isfall_all = {};
    conf_all = {};
    
    if featureset
        Ft = X_test(:,5:end);
    else
        Ft = X_test(:,943:962); % only magnitude features
    end
    
    indtrain = any(bsxfun(@eq,subjid(:,2),locations_train),2);
    
    falls_ind = find(indtrain & L);
    nfalls_ind = find(indtrain & ~L);
    
    for run = 1:Nruns
        indtrain = [falls_ind(randperm(length(falls_ind),nData)); nfalls_ind(randperm(length(nfalls_ind),nData))];
        
        if model
            [fvar_si,b,nz_ind]=Modeltrain(F(indtrain,:),L(indtrain),alpha,lambda,0);
            % Testing the model on each amputee subj (external set)
            for subj = 1:length(unique(X_test(:,1)))
                rowid = (X_test(:,1)==subj) & any(bsxfun(@eq,X_test(:,2),locations_test),2);
                [pred,conf,confmat] = Modeleval(Ft(rowid,:),X_test(rowid,4)<5,fvar_si,nz_ind,b,0.5,0);
                conf_all{subj}=conf;
                confmat_all(:,:,subj)=confmat;
                isfall = X_test(rowid,4)<5;
                isfall_all{subj}=isfall;
                if length(unique(isfall)) <2 %subject AF06 with missing activities - TO FIX
                    ;
                else
                    [~, ~, ~, AUC_HA(run,subj)]=perfcurve(isfall, conf, true);
                    Sens_HA(run,subj) = sum(pred & isfall)/sum(isfall);
                    Spec_HA(run,subj) = sum(~pred & ~isfall)/sum(~isfall);
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
    figure(figroc), hold on
    [X, Y, T, AUC]=perfcurve(isfall_all(2:end), conf_all(2:end), true,'XVals',[0:0.05:1]);
%     [X, Y, T, AUC]=perfcurve(cell2mat(isfall_all'), cell2mat(conf_all'), true,'Nboot',0,'XVals',[0:0.05:1]); %cb with Bootstrap
    e = errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
    e.LineWidth = 2; e.Marker = 'o';

    
    %% LOSO on Amputees
    AUC_AA=NaN(Nruns,length(subj),'double');
    Sens_AA=NaN(Nruns,length(subj),'double');
    Spec_AA=NaN(Nruns,length(subj),'double');
    isfall_all = {};
    conf_all = {};

    if featureset
        F = X_test(:,5:end);
    else
        F = X_test(:,943:962); % only magnitude features
    end
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
                    ;
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
    figure(figroc), hold on
    [X, Y, T, AUC]=perfcurve(isfall_all(2:end), conf_all(2:end), true,'XVals',[0:0.05:1]);
%     [X, Y, T, AUC]=perfcurve(cell2mat(isfall_all'),cell2mat(conf_all'),true,'Nboot',0,'XVals',[0:0.05:1]);
    e = errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
    e.LineWidth = 2; e.Marker = 'o';
    legend('Healthy-Healthy','Healthy-Amputee','Amputee-Amputee')

end

AUC = {AUC_HH;AUC_AA;AUC_HA};
Sens = {Sens_HH;Sens_AA;Sens_HA};
Spec = {Spec_HH;Spec_AA;Spec_HA};


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
