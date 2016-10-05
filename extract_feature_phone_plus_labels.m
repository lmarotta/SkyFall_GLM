function [FN, fl, L, fold, mu, stdF, eps, nzstd]= extract_feature_phone_plus_labels(labels, fold_nr)

eps=1e-6;

sub_ind_vect=[];
%----------------------------------------------------------------------------------------------------
% extracts features from  5 seconds interval of  sensor recordings (accelerometer, gyro and barometer)
%----------------------------------------------------------------------------------------------------
% INPUTS
%         gyro_data : Tg * 4 matrix, Tg number of time stamps, 1st column: time stamps, 2nd, 3rd and 4th are the the sensor readings (normalized time stamps)
%         accel_data: Ta * 4 matrix, Ta number of time stamps, 1st column: time stamps, 2nd, 3rd and 4th are the the sensor readings (normalized time stamps)
%         baro_data : Tb * 3 matrix, Tb number of time stamps, 1st column: time stamps, 2nd and 3rd are  the sensor readings (normalized time stamps)
%
%         fvar: contains the required vaiables for buliding  features  with
%         the following fields
%                fvar.eps: the resolution threshold
%                fvar.mean: mean of the training set data
%                fvar.std:  standard deviation of training set data
%                fvar.nzstd: columns with non-zero standard deviation in training set
%
% OUTPUT
%         FN: 1*p vector containing  Features, where p is the size of feature set is returned empty

javaaddpath('purple-robot-skyfall.jar');

save labels_struct labels

com.company.TrainFeatureExtractor.extractFeatures('labels_struct.mat', 'full_features.csv');

F=csvread('full_features.csv');
F=F(:,1:1781);

L=labels.value;

% find zeroed rows in F
z_inds=max(abs(F).')==0;

F(z_inds,:)=[];
L(z_inds)=[];
subjects=labels.subject;
subjects(z_inds)=[];

for i=1:length(subjects)

    Amp_flag =regexp(subjects{i},'AF','ignorecase');

    if isempty(Amp_flag)
        sub_ind=str2num(subjects{i}(3:4))+7;
    else
        sub_ind=str2num(subjects{i}(3:4));
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
fl = [];
stdF=stdF(nzstd);

FN=(F-v1*mu)./(v1*stdF);

fold=zeros(dsz,1);


[~,id]=sort(sub_ind_vect,'ascend');
FN= FN(id, :);
L= L(id);

%fold_ind= linspace(1,dsz+1,folds_nr);
%fold_sz= floor(dsz/fold_nr);

fold_grid=round(linspace(1,dsz,fold_nr+1));



for i=1:fold_nr
 fold(fold_grid(i):fold_grid(i+1))=i;
end;


return;