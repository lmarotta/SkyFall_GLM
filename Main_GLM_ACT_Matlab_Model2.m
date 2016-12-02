load class_params_Stage1Model.mat
load AmputeeHome.mat %this is the home datafile to train on
epsF = 1e-6; %threshold on std dev for standardizing features

%evaluate model on home data to generate FP
[pred,conf,confmat] = Modeleval(F,L,fvar,nz_ind,b,Thres,0);
figure, plot(pred,'o-r')
FP = length(find(pred)) %number of FP on home data
indFP = find(pred);

%% Train Model 2 on Home FP + Lab Falls, Evaluate on Amputee
load Training_Data.mat
Fhome = F(indFP,:);
Lhome = zeros(size(Fhome,1),1);
L = F(:,4); %labels for lab data
ind = F(:,3) == 0;
FAmputee = F(ind,5:end);    %features from amputee only (test)
LAmputee = L(ind);      %labels for amputees
L = L(~ind);            %labels for all others
F = F(~ind,:);
%use only falls data from lab
ind = L<9;  %only falls data
L = ones(sum(ind),1); 
F = F(ind,5:end);

F = [Fhome;F];
L = [Lhome;L];

alpha = 0.6;
lambda = 0.015;
[fvar,b,nz_ind]=Modeltrain(F,L,alpha,lambda,no_baro); %Model 2
[pred,conf,confmat] = Modeleval(FAmputee,LAmputee,fvar,nz_ind,b);
isfall = LAmputee<9;
[X, Y, T, AUC]=perfcurve(isfall, conf, true);
figure; plot(X,Y)

%%
[pred1,conf1] = Modeleval(FAmputee);
[X, Y, T, AUC]=perfcurve(isfall, conf1, true);
figure; plot(X,Y)


conf_combined=conf(pred1==1);
[X, Y, T, AUC]=perfcurve(isfall(pred1==1), conf_combined, true);
figure; plot(X,Y)


