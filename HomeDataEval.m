clear all

features_used = ones(18,1); %features_used([4 7 10:18]) = 0; %full feature set
% features_used = zeros(18,1); features_used([8]) = 1; %only magnitude features
featureInds=getFeatureInds(features_used); 

%default values - no grid search over params
alpha = 0.6;
lambda = 0.015;
no_baro=0; % 0 - use barometer

if ~exist('HealthyData.mat','file')
    TrainingDataSetup([], [], 10, 0, 'HealthyData');
else
    X = load('HealthyData.mat');
    X = X.F;
end
%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = X(:,4); %labels for classification
L=(L<9);%+1;  %binary labeling
subjid = X(:,1:3);  %[subjid, location subjcode]
subj=unique(subjid(:,1));

F = X(:,5:end); %remove meta-data
F = F(:,featureInds);

disp('Training model on all healthy data')
[fvar,b,nz_ind]=Modeltrain(F,L,alpha,lambda,no_baro);
Flab = F; 
Llab = L;

%% Test on amputee lab
Thres = 0.5; %model threshold
disp('Test model on amputees - lab data')
display = 1;
X_Amp = load('Test_Data_Amputees');
X_Amp = X_Amp.F;
L = X_Amp(:,4); L = L<9;
X_Amp(:,1:4) = []; %remove meta-data
F = X_Amp(:,featureInds);
[pred,conf,confmat] = Modeleval(F,L,fvar,nz_ind,b,Thres,display);
conf_lab = conf;
L_lab = L;


%% Load Home data and test - Amputee
display = 0;
F = [];
lf = [];
filespath = 'Z:/Amputee Phones-R01/Home Data Collection/Amputees/';
if ~exist([filespath 'HomeDataAmp.mat'],'file')
    f = dir([filespath,'/Raw/*.mat']);
    for i = 1:length(f)
        sprintf('Filename %s',f(i).name)
        l = load([filespath '/Raw/' f(i).name]); labels = l.labels;
        %extract features
        Fnew = HomeDataSetup(labels,1.5); %1.5g threshold for acceleration clips
        F = [F;Fnew];
        lf = [lf;size(Fnew,1)*5/60/60]; %store duration of each day
        sprintf('Data length = %.2f h',size(F,1)*5/60/60) 
        disp(lf)
    end
    save([filespath 'HomeDataAmp.mat'],'F')

else
    l = load([filespath 'HomeDataAmp.mat']);
    F = l.F;
    sprintf('Data length = %.2f h',size(F,1)*5/60/60)
end

%feature selection
F = F(:,featureInds);
L = false(size(F,1),1); %all labels are non-fall
[pred,conf,confmat] = Modeleval(F,L,fvar,nz_ind,b,Thres,display);
sprintf('Spec = %.2f%', length(pred)/(length(pred)+sum(pred)))
figure, histogram(conf)    
figure, histogram(pred), hold on, title(sprintf('Spec = %.2f%', length(pred)/(length(pred)+sum(pred))))

roc = figure, hold on
[X, Y, T, AUC]=perfcurve([L;L_lab], [conf;conf_lab], true,'XVals',[0:0.05:1]); %conf bounds with CV
% [X, Y, T, AUC]=perfcurve(cell2mat(isfall_all'), cell2mat(conf_all'), true,'Nboot',0,'XVals',[0:0.05:1]); %cb with Bootstrap
e = plot(X,Y);
e.LineWidth = 2; e.Marker = 'o';

%% Train on Healthy Lab + home data misclassified clips
display = 0;
F = [];
lf = [];
filespath = 'Z:/Amputee Phones-R01/Home Data Collection/Healthy/';
if ~exist([filespath 'HomeDataHealthy.mat'],'file')
    f = dir([filespath,'/Raw/*.mat']);
    for i = 1:length(f)
        sprintf('Filename %s',f(i).name)
        l = load([filespath '/Raw/' f(i).name]); labels = l.labels;
        %extract features
        Fnew = HomeDataSetup(labels,1.5); %1.5g threshold for acceleration clips
        F = [F;Fnew];
        lf = [lf;size(Fnew,1)*5/60/60]; %store duration of each day
        sprintf('Data length = %.2f h',size(F,1)*5/60/60) 
        disp(lf)
    end
    save([filespath 'HomeDataHealthy.mat'],'F')

else
    sprintf('Load %s',[filespath 'HomeDataHealthy.mat'])
    l = load([filespath 'HomeDataHealthy.mat']);
    F = l.F;
    sprintf('Data length = %.2f h',size(F,1)*5/60/60)
end

%feature selection
F = F(:,featureInds);
L = zeros(size(F,1),1); %all labels are non-fall
[pred,conf,confmat] = Modeleval(F,L,fvar,nz_ind,b,Thres,display);
sprintf('Spec = %.2f%', length(pred)/(length(pred)+sum(pred)))
indmisc=find(pred); %misclassified clips
sprintf('# of misclassified clips = %d/%d',length(indmisc),length(F))
Fmisc = F(indmisc,:); Lmisc = zeros(length(indmisc),1);
F = [Flab;Fmisc]; L = [Llab;Lmisc];
disp('Training model on Healthy lab + misclassified home data')
[fvar,b,nz_ind]=Modeltrain(F,L,alpha,lambda,no_baro);

%test model on amputee lab dataset
Thres = 0.5; %model threshold
disp('Test model on amputees - lab data')
display = 1;
X_Amp = load('Test_Data_Amputees');
X_Amp = X_Amp.F;
L = X_Amp(:,4); L = L<9;
X_Amp(:,1:4) = []; %remove meta-data
F = X_Amp(:,featureInds);
[pred,conf,confmat] = Modeleval(F,L,fvar,nz_ind,b,Thres,display);
conf_lab = conf;
L_lab = L;

%test on amputee home data
display = 0;
filespath = 'Z:/Amputee Phones-R01/Home Data Collection/Amputees/';
disp('Test model on amputees - home data')
l = load([filespath 'HomeDataAmp.mat']);
F = l.F;
F = F(:,featureInds);
L = false(size(F,1),1);
[pred,conf,confmat] = Modeleval(F,L,fvar,nz_ind,b,Thres,display);
sprintf('Spec = %.2f%', length(pred)/(length(pred)+sum(pred)))
figure, histogram(conf)    
figure, histogram(pred), title(sprintf('Spec = %.2f%', length(pred)/(length(pred)+sum(pred))))

%test on home + lab data
figure(roc)
[X, Y, T, AUC]=perfcurve([L;L_lab], [conf;conf_lab], true,'XVals',[0:0.05:1]); %conf bounds with CV
% [X, Y, T, AUC]=perfcurve(cell2mat(isfall_all'), cell2mat(conf_all'), true,'Nboot',0,'XVals',[0:0.05:1]); %cb with Bootstrap
e = plot(X,Y);
e.LineWidth = 2; e.Marker = 'o';
xlabel('False positive rate')
ylabel('True positive rate')

%% Train on Healthy Lab + home data
load ../SkyFall_HomeData/Nick_Luca_01132017/NickHomeData.mat
%extract features 
F = [];
F = HomeDataSetup(labelsNick,1.5); %1.5g threshold for acceleration clips
if full_featureset
    F=F(:,1:end);
else
    F = F(:,939:958); % only magnitude features
end
L = zeros(size(F,1),1);

F = [Flab;F]; L = [Llab;L];
disp('Training model on lab + home data')
[fvar,b,nz_ind]=Modeltrain(F,L,alpha,lambda,no_baro);

%% Load Home data and test
Thres = 0.5;
display = 0;
%extract features 
F = [];
F = HomeDataSetup(labelsLuca,1.5); %1.5g threshold for acceleration clips
if full_featureset
    F=F(:,1:end);
else
    F = F(:,939:958); % only magnitude features
end
L = zeros(size(F,1),1);
[pred,conf,confmat] = Modeleval(F,L,fvar,nz_ind,b,Thres,display);
sprintf('Spec = %.2f%', length(pred)/(length(pred)+sum(pred)))

figure, histogram(conf)    
figure, histogram(pred)

%% Test on amputee lab

display = 1;
X_Amp = load('Test_Data_Amputees');
X_Amp = X_Amp.F;
L = X_Amp(:,4); L = L<9;
X_Amp(:,1:4) = [];
if full_featureset
    F=X_Amp(:,1:end);
else
    F = X_Amp(:,939:958); % only magnitude features
end
[pred,conf,confmat] = Modeleval(F,L,fvar,nz_ind,b,Thres,display);

%% Test phone model w Gyro features
load ./PhoneModels/MagFeat.mat
F = HomeDataSetup(labels,1.5); %1.5g threshold for acceleration clips
F=F(:,1:end);

L=false(size(F,1),1);

FNZ = F(:,fvar.nzstd);
FN = (FNZ - repmat(fvar.mu,[size(FNZ,1),1])) ./ repmat(fvar.std,[size(FNZ,1),1]); %features normalized
FN = FN(:,nz_ind);
conf= glmval(b, FN, 'logit');
pred= ceil(conf-Thres);

%results
isfall = logical(L);
confmat(:,:)=confusionmat(isfall,pred==1,'order',[false true]);

figure, histogram(conf)    
figure, histogram(pred)
    



