
function [FN, L, fold, mu, stdF, eps, nzstd]= extract_features(labels, fold_nr)

%----------------------------------------------------------------------------------------------------
% extracts features from  5 seconds interval of  sensor recordings (accelerometer, gyro and barometer)
%----------------------------------------------------------------------------------------------------
% INPUTS
%         data structure
%
% OUTPUT
%         FN: N x p feature matrix where N is number of datapoints and p is the number of features
%         L: numeric labels, N x 1 vector  
%         fold: the fold each feature vector (row) belongs (N x 1)
%         mu,stdF: mean and standard deviation of features (1 x p)
%         eps: the zero std threshold 
%         nzstd: features with non-zero standard deviation in training set (1 x p)
eps=1e-6;

sub_ind_vect=[];

javaaddpath('FeatureGeneration/purple-robot-skyfall.jar');

save labels_struct labels

com.company.TrainFeatureExtractor.extractFeatures('labels_struct.mat', 'full_features.csv');

F=csvread('full_features.csv');
F=F(:,1:1781);

L=labels.value;

% find zeroed rows in F (cannot compute features, missing data)
z_inds=max(abs(F).')==0;

F(z_inds,:)=[];
L(z_inds)=[];
subjects=labels.subject;
subjects(z_inds)=[];

for i=1:length(subjects)

    Amp_flag =regexp(subjects{i},'AF','ignorecase');

    if isempty(Amp_flag)
        sub_ind=str2num(subjects{i}(3:4))+7; %Control ID is higher than 7
    else
        sub_ind=str2num(subjects{i}(3:4)); %Amputee ID is 1-7
    end;

    sub_ind_vect = [sub_ind_vect; sub_ind];

end


% F=F(1:end-50,:);
% sub_ind_vect=sub_ind_vect(1:end-50);
% L=L(1:end-50);
dsz= size(F,1);
v1= ones(dsz,1);
mu= mean(F);
stdF= std(F);


%%%%Standardize the result

nzstd=stdF>eps*mean(stdF);
% remove barometer features
nzstd(1701:1781)=zeros(1,81);
% remove time mean
nzstd([1 851])=[0 0];
mu=mu(nzstd);
F=F(:,nzstd);
stdF=stdF(nzstd);

FN=(F-v1*mu)./(v1*stdF);

fold=zeros(dsz,1);


[~,id]=sort(sub_ind_vect,'ascend');
FN= FN(id, :);
L= L(id);

%fold_ind= linspace(1,dsz+1,folds_nr);
%fold_sz= floor(dsz/fold_nr);

fold_grid=round(linspace(1,dsz,fold_nr+1));

% Arranging data from all subjects in 10 folds
for i=1:fold_nr
 fold(fold_grid(i):fold_grid(i+1))=i;
end;


return;