%% Testing the model on Test Data 
clear all

load TestData_Home_unlab.mat
load class_params_ACT_nobar_NoFalllike.mat
%include labels if not there
N = length(labels.winsize); %number of clips
labels.value = 9*ones(N,1);   %only non-falls (labels >=9)

javaaddpath('../purple-robot-skyfall.jar');
save labels_test labels
com.company.TrainFeatureExtractor.extractFeatures('labels_test.mat', 'test_features.csv');
F=csvread('test_features.csv');rple-robot-skyfall.j
F=F(:,1:1781);

% find zeroed rows in F (feature cannot be calculated)
z_inds=max(abs(F).')==0;

F(z_inds,:)=[];
F=F(:,fvar.nzstd);
dsz= size(F,1);
v1= ones(dsz,1);
FN=(F-v1*fvar.mu)./(v1*fvar.std);
FN_nz=FN(:,nz_ind);

conf= glmval(b, FN_nz, 'logit');
id= ceil(conf-0.5); %0=non fall, 1=fall

%plot the data
figure
plot(((1:length(id))*5/60),id), ylim([-0.1 1.1])
xlabel('Time [min]'), ylabel('IsFall')

isfall = labels.value < 9;
isfall = isfall(~z_inds);
%fall_err = sum(~id(isfall))/sum(isfall);
%nfall_err = sum(id(~isfall))/sum(~isfall);
% err = (sum(id.' ~= isfall))/length(id); %error rate
% confmat(:,:)=confusionmat(isfall,id==1,'order',[1 2]);
% figure; imagesc(confmat./repmat(sum(confmat,2),[1 2])); colorbar; caxis([0 1])


