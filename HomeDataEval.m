clear all

full_featureset = 0;
%default values - no grid search over params
alpha = 0.6;
lambda = 0.015;
no_baro=0; % 0 - use barometer

X = load('HealthyData.mat');
X = X.F;
%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = X(:,4); %labels for classification
L=(L<9);%+1;  %binary labeling
subjid = X(:,1:3);  %[subjid, location subjcode]
if full_featureset
    F=X(:,5:end);
else
    F = X(:,943:962); % only magnitude features
end
subj=unique(subjid(:,1));

disp('Training model on all healthy data')
[fvar,b,nz_ind]=Modeltrain(F,L,alpha,lambda,no_baro);
Flab = F; 
Llab = L;

%% Load Home data and test
Thres = 0.5;
display = 0;
load ../SkyFall_HomeData/Nick_Luca_01132017/LucaHomeData.mat
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

figure, histogram(conf)    
figure, histogram(pred)

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
F = HomeDataSetup(labelsLuca,1.5); %1.5g threshold for acceleration clips
F=F(:,1:end);

[pred,conf,confmat] = Modeleval(F,L,fvar,nz_ind,b,Thres,display);
figure, histogram(conf)    
figure, histogram(pred)
    


