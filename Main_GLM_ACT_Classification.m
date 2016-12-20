%Current Labels
%slip=1
%trip=2;
%rightfall=3;
%leftfall=4;
%Activities=9; (everything else)

%feature matrix structure
% F = [subj_id location subjcode labels.value Features];
%Location: (0 = NA, 1=Pouch, 2=Pocket, 3=Hand)
%Subjcode: 1 = CF (healthy); 0 = AF (amputee)

FP=[]; %to store FP and train model stage 2
amputee_CV = 0;   %use only amputee data for CV (train and test)
epsF = 1e-6; %threshold on std dev for standardizing features
class=0; % flag for fall classification (rather than detection only)
no_baro=0; % 0 - use barometer
nbaro_features = 6;
RF=1;


alpha = 0.6;
lambda = 0.015;

% rng(10001)
%% Loading the data

load Training_Data
%Feature matrix assembled as follows
% F = [subj_id location subjcode labels Features];

%convert labels to binary (1-4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = F(:,4); %labels for classification
LB=(L<9);%+1;  %binary labeling

% Remove non-falls for classification
F(~LB,:)=[];
L(~LB)=[];
%combine left and right together
% L(L==4)=3;
% L(L==2)=1;
% L(L==3)=2;
L(L==9)=5; %non-falls

if ~amputee_CV %train on healthy and test on amputee
    ind = F(:,3) == 0;      %index of amputee subjects (test)
    FAmputee= F(ind,5:end);
    subjid = F(~ind,1:3);  %[subjid, location subjcode]
    F = F(~ind,5:end);     %features for training
    LAmputee = L(ind);      %labels for amputees
    L = L(~ind);            %labels for all others
else %CV on amputee only
    ind = F(:,3) == 0;      %index of amputee subjects
    subjid = F(ind,1:3);  %[subjid, location subjcode]
    F = F(ind,5:end);
    L = L(ind);      %labels for amputees
    
end

subj=unique(subjid(:,1));

%% LOSO CV
%This CV is used to find an optimal threshold to achieve desired
%Sensitivity

%initialize fvar structure
f = cell(length(subj),1);
fvar= struct('std',f,'mu',f,'eps',f,'nzstd',f);

if amputee_CV %LOSO CV on amputees
    
    for indCV=1:length(subj)
        
        test_subj=subj(indCV);
        indtrain = subjid(:,1)~=test_subj;
        indtest = ~indtrain;
        
        disp(['Training model with subject' num2str(indCV) 'out'])
        disp(['Training Data size ' num2str(size(F(indtrain,:),1)) ' instances - # Features ' num2str(size(F(indtrain,:),2))])
        
        if RF
            RFModel=TreeBagger(50,F(indtrain,:),L(indtrain),'OOBPred','on');
            figure; plot(oobError(RFModel))
            %         t = templateTree('MinLeafSize',5);
            %         RFModel=fitensemble(F(indtrain,:),L(indtrain),'RUSBoost',20,t,'LearnRate',0.1);
        else
            L1=L; L2=L; L3=L; L4=L;
            L1(L~=1)=2;
            L2(L==2)=1; L2(L~=2)=2;
            L3(L==3)=1; L3(L~=3)=2;
            L4(L==4)=1; L4(L~=4)=2;
            
            Model1=lassoglm(F(indtrain,:),L1(indtrain),'binomial','Alpha',alpha,'Lambda',lambda);
            Model2=lassoglm
        end
        
        
        %Evaluate model on left out AMPUTEE
        P=predict(RFModel,F(indtest,:));
        confmat=confusionmat(L(indtest),str2double(P),'order',1:length(unique(L)));
        figure, imagesc(confmat./repmat(sum(confmat,2),[1 length(unique(L))]));

    end
    
else %Train on all healthy and test on amputees
    disp('Training model on all healthy data')
    disp(['Training Data size ' num2str(size(F,1)) ' instances - # Features ' num2str(size(F,2))])

    RFModel=TreeBagger(30,F,L,'OOBPred','on');
    figure; plot(oobError(RFModel))
    %     t = templateTree('MinLeafSize',5);
    %     RFModel=fitensemble(F,L,'RUSBoost',20,t,'LearnRate',0.1);
    
    %Evaluate model on AMPUTEES
    P=predict(RFModel,FAmputee);
    P = str2double(P);
    confmat=confusionmat(LAmputee,P,'order',1:length(unique(L)));
    figure, imagesc(confmat./repmat(sum(confmat,2),[1 length(unique(L))]));
    acc = sum(P == LAmputee)/length(LAmputee)
    title(['Accuracy = ' num2str(acc)])
    
end