%%
%load LabeledDataNov9
%slip=1
%trip=2;
%rightfall=3;
%leftfall=4;
%sliplike=5;
%triplike=6;
%rightfalllike=7;
%leftfalllike=8;
%walk_start=9;
%walk_end=10;
%stairsup_start=11;
%stairsup_end=12;
%stairsdown_start=13;
%stairsdown_end=14;
%sit_start=15
%sit_end=16
%stand_start=17
%stand_end=18;

clear all

X=70; % set percentage to save for training if split==1
split=1; % flag to split data into test and train sets (25-75) and create cofnusion matrix

for indSet=1:2
skip_like=0; % flag to not use fall-like data
ACTnumber=1000;
rmv_ACT=2-indSet; % flag to remove activities data
class=2-indSet;

addpath(genpath('./glmnet_matlab/'))
%% Initialization

Twin = 2000;
folds_nr=10;
overlap=500;
trNum=100;
clNum=5;
nalpha=10;
nlambda=200;
maxPrincRatio=1e-4;


%%% Loading the data

% load full data set
% any filtering of data past this point should be mirrored on the features
% loaded by extract_feature_phone_plus_labels
load labels_full.mat

% Can use either standard or jittered falls data
% load labels_plus_data.mat
% load labels_plus_data_jittered_Unif
% 
% % load Activities data - used as non falls
% ACT=load('labels_plus_data_ACT.mat');
% ACTnumber=min([ACTnumber length(ACT.labels.subject)]);
% 
% %% Merge Falls and Activities Data
% 
% labels.value=[labels.value ACT.labels.value(1:ACTnumber)];
% labels.subject=[labels.subject ACT.labels.subject(1:ACTnumber)];
% labels.timestamp=[labels.timestamp ACT.labels.timestamp(1:ACTnumber)];
% labels.text=[labels.text ACT.labels.text(1:ACTnumber)];
% labels.acce=[labels.acce ACT.labels.acce(1:ACTnumber)];
% labels.baro=[labels.baro ACT.labels.baro(1:ACTnumber)];
% labels.gyro=[labels.gyro ACT.labels.gyro(1:ACTnumber)];
% if skip_like
%     nlike_inds=find(labels.value<5 | labels.value>8);
% 
%     labels.value=labels.value(nlike_inds);
%     labels.subject=labels.subject(nlike_inds);
%     labels.timestamp=labels.timestamp(nlike_inds);
%     labels.text=labels.text(nlike_inds);
%     labels.acce=labels.acce(nlike_inds);
%     labels.baro=labels.baro(nlike_inds);
%     labels.gyro=labels.gyro(nlike_inds);
% end
% 
if split
    %Take 75% of original data for training 
    inds = randperm(length(labels.timestamp));
    indtrain = inds(1:round(X/100*length(inds)));
    indtest = inds(round(X/100*length(inds))+1:end);
else
    indtrain = 1:length(labels.timestamp);
    indtest=[];
end

if rmv_ACT
    train_ACT=labels.value(indtrain)==9;
    indtrain(train_ACT)=[];
    
    test_ACT=labels.value(indtest)==9;
    indtest(test_ACT)=[]; 
end

% 
labels_test=labels;
% 
%TEST DATA
labels_test.timestamp = labels.timestamp(indtest);
labels_test.value = labels.value(indtest);
labels_test.subject = labels.subject(indtest);
labels_test.text = labels.text(indtest);
labels_test.acce = labels.acce(indtest);
labels_test.gyro = labels.gyro(indtest);
labels_test.baro = labels.baro(indtest);
% 
inds=cellfun(@(x,y,z) length(x)<100 | length(y)<100 | length(z)<10,labels_test.acce,labels_test.gyro,labels_test.baro);
inds=~inds;
labels_test.timestamp = labels_test.timestamp(inds);
labels_test.value = labels_test.value(inds);
labels_test.subject = labels_test.subject(inds);
labels_test.text = labels_test.text(inds);
labels_test.acce = labels_test.acce(inds);
labels_test.gyro = labels_test.gyro(inds);
labels_test.baro = labels_test.baro(inds);
% 
% % save indstest_act.mat indstest labels
% 
%TRAIN DATA
labels.timestamp = labels.timestamp(indtrain);
labels.value = labels.value(indtrain);
labels.subject = labels.subject(indtrain);
labels.text = labels.text(indtrain);
labels.acce = labels.acce(indtrain);
labels.gyro = labels.gyro(indtrain);
labels.baro = labels.baro(indtrain);

%% Feature extraction

[FN, fl, L, fold_id, muF, stdF, epsF, nzstd]=extract_feature_phone_plus_labels(labels, folds_nr, indtrain);
fvar.std= stdF;
fvar.mu= muF;
fvar.eps= epsF;
fvar.nzstd= nzstd;
fvar.fl = fl;

FSz= size(FN,2);
DSz= size(FN,1);

% Assign fall categories as 2 (falls and fall-like) or 1 (non-fall)
LB=(L<9)+1;
LF=L;
LF(L>4 & L<9)=5;
%uncomment for fall-like/activities separate 
% LF(L>4)=5;
%uncomment for fall detection
if ~class
    LF=LB;
end


alpha=linspace(0.5, 1, nalpha);
alpha=alpha(end:-1:1);
d=cell(nalpha, 1);
opts=glmnetSet;
opts.nlambda= nlambda;


min_err_iter= zeros(nalpha, 1);
lam_vect= zeros(nalpha, 1);
lam_ind =zeros(nalpha, 1);


for i=1: nalpha
    opts.alpha= alpha(i);
    
    if max(LF)>2
        d{i}= cvglmnet(FN,LF,'multinomial',opts,'class',folds_nr,fold_id);
    else
        d{i}= cvglmnet(FN,LF,'binomial',opts,'class',folds_nr,fold_id);
    end;
    
    %cvmMat{i=d{i}.cvm;
    err_now= min(d{i}.cvm)
    alp_now= alpha(i)
    min_err_iter(i)= err_now;
    lam_vect(i)= d{i}.lambda_min;
    [~, lam_ind(i)]= min(d{i}.lambda);
end;

[min_err, alp_ind]= min(min_err_iter);
%[min_err,alp_ind]=  min(m_col);


accuracy=(1-min_err)*100
std= d{alp_ind}.cvsd(lam_ind(alp_ind))*100
alp_opt= alpha(alp_ind)
lam_opt= d{alp_ind}.lambda(lam_ind(alp_ind))
sparsity= 100*(FSz-d{alp_ind}.nzero(lam_ind(alp_ind)))/FSz

pred= cvglmnetPredict(d{alp_ind},[],lam_opt,'nonzero');

if max(LF)>2
pred_Mat=cell2mat(pred);
else
pred_Mat= pred;
end;

pred_sum= sum(pred_Mat,2);

nz_ind= pred_sum>0;
FN_nz= FN(:,nz_ind);

min_err_iter_nz= zeros(nalpha, 1);
lam_vect_nz= zeros(nalpha, 1);
lam_ind_nz= zeros(nalpha, 1);
nalpha=40
alpha=linspace(0, 0.5, nalpha);
alpha=alpha(end:-1:1);
for i=1: nalpha
    opts.alpha= alpha(i);
    if max(LF)>2
        d_nz{i}= cvglmnet(FN_nz, LF, 'multinomial', opts, 'class', folds_nr, fold_id);
    else
        d_nz{i}= cvglmnet(FN_nz, LF, 'binomial', opts, 'class', folds_nr, fold_id);
    end;
    %cvmMat{i=d{i}.cvm;
    err_now= min(d_nz{i}.cvm)
    alp_now= alpha(i)
    min_err_iter_nz(i)= err_now;
    lam_vect_nz(i)= d_nz{i}.lambda_min;
    [~, lam_ind_nz(i)]= min(d_nz{i}.lambda);
end;
[min_err_nz, alp_ind_nz]= min(min_err_iter_nz);
%[min_err,alp_ind]=  min(m_col);


accuracy=(1-min_err_nz)*100
std= d_nz{alp_ind_nz}.cvsd(lam_ind_nz(alp_ind_nz))*100
alp_opt_nz= alpha(alp_ind_nz)
lam_opt_nz= d_nz{alp_ind_nz}.lambda(lam_ind_nz(alp_ind_nz))
%sparsity= 100*(FSz-d{alp_ind}.nzero(lam_ind(alp_ind)))/FSz

pred= cvglmnetPredict(d_nz{alp_ind_nz},[],lam_opt_nz,'nonzero');
b= glmnetPredict(d_nz{alp_ind_nz}.glmnet_fit,[],lam_opt_nz,'coefficients');
    
%end;
%% Genrating Random n folds
%DSz=size(FN,1);
%fldsz=round(DSz/folds);
if split
    F=csvread('full_feature_set.csv');
    F=F(indtest,1:1781);

    % find zeroed rows in F
    z_inds=max(abs(F).')==0;

    F(z_inds,:)=[];
    indtest(z_inds)=[];
    F=F(:,fvar.nzstd);
    dsz= size(F,1);
    v1= ones(dsz,1);
    FN=(F-v1*fvar.mu)./(v1*fvar.std);
    FN_nz=FN(:,nz_ind);
    
    if class
        conf(:,1)=glmval(b{1},FN_nz,'logit');
        conf(:,2)=glmval(b{2},FN_nz,'logit');
        conf(:,3)=glmval(b{3},FN_nz,'logit');
        conf(:,4)=glmval(b{4},FN_nz,'logit');
        conf(:,5)=glmval(b{5},FN_nz,'logit');
        [~, id]=max(conf,[],2);
        
        value=labels_test.value;
        value(value>5)=5;
        confmat=confusionmat(id,value);
        figure; imagesc(confmat./repmat(sum(confmat,2),[1 5])); colorbar; caxis([0 1])
    else
        conf= glmval(b, FN_nz, 'logit');
        id= ceil(conf-0.5);
        
        isfall = labels_test.value < 9;
        fall_err = sum(~id(isfall))/sum(isfall);
        nfall_err = sum(id(~isfall))/sum(~isfall);
        err = (sum(id.' ~= isfall))/length(id); %error rate
        confmat=confusionmat(isfall,id==1);
        figure; imagesc(confmat./repmat(sum(confmat,2),[1 2])); colorbar; caxis([0 1])
    end
    
end

all_fvar{indSet}=fvar;
all_b{indSet}=b;
all_nz_ind{indSet}=nz_ind;
all_conf{indSet}=confmat;
end

labels_2_train.value=labels_test.value(isfall);
labels_2_train.text=labels_test.text(isfall);
labels_2_train.timestamp=labels_test.timestamp(isfall);
labels_2_train.subject=labels_test.subject(isfall);
labels_2_train.acce=labels_test.acce(isfall);
labels_2_train.gyro=labels_test.gyro(isfall);
labels_2_train.baro=labels_test.baro(isfall);

indtest2=indtest;

labels=labels_2_train;

Twin = 2000;
folds_nr=10;
overlap=500;
trNum=100;
clNum=5;
nalpha=10;
nlambda=200;
maxPrincRatio=1e-4;

if split
    %Take X% of original data for training 
    inds = randperm(length(labels.timestamp));
    indtrain = inds(1:round(X/100*length(inds)));
    indtest = inds(round(X/100*length(inds))+1:end);
else
    indtrain = 1:length(labels.timestamp);
    indtest=[];
end

if rmv_ACT
    train_ACT=labels.value(indtrain)==9;
    indtrain(train_ACT)=[];
    
    test_ACT=labels.value(indtest)==9;
    indtest(test_ACT)=[]; 
end

% 
labels_test=labels;
% 
%TEST DATA
labels_test.timestamp = labels.timestamp(indtest);
labels_test.value = labels.value(indtest);
labels_test.subject = labels.subject(indtest);
labels_test.text = labels.text(indtest);
labels_test.acce = labels.acce(indtest);
labels_test.gyro = labels.gyro(indtest);
labels_test.baro = labels.baro(indtest);
% 
inds=cellfun(@(x,y,z) length(x)<100 | length(y)<100 | length(z)<10,labels_test.acce,labels_test.gyro,labels_test.baro);
inds=~inds;
labels_test.timestamp = labels_test.timestamp(inds);
labels_test.value = labels_test.value(inds);
labels_test.subject = labels_test.subject(inds);
labels_test.text = labels_test.text(inds);
labels_test.acce = labels_test.acce(inds);
labels_test.gyro = labels_test.gyro(inds);
labels_test.baro = labels_test.baro(inds);
% 
% % save indstest_act.mat indstest labels
% 
%TRAIN DATA
labels.timestamp = labels.timestamp(indtrain);
labels.value = labels.value(indtrain);
labels.subject = labels.subject(indtrain);
labels.text = labels.text(indtrain);
labels.acce = labels.acce(indtrain);
labels.gyro = labels.gyro(indtrain);
labels.baro = labels.baro(indtrain);

%% Feature extraction

vals=indtest2(isfall);
[FN, fl, L, fold_id, muF, stdF, epsF, nzstd]=extract_feature_phone_plus_labels(labels, folds_nr, vals(indtrain));
fvar.std= stdF;
fvar.mu= muF;
fvar.eps= epsF;
fvar.nzstd= nzstd;
fvar.fl = fl;

FSz= size(FN,2);
DSz= size(FN,1);

% Assign fall categories as 2 (falls and fall-like) or 1 (non-fall)
LB=(L<9)+1;
LF=L;
LF(L>4 & L<9)=5;
%uncomment for fall-like/activities together 
LF(L>4)=5;

alpha=linspace(0.5, 1, nalpha);
alpha=alpha(end:-1:1);
d=cell(nalpha, 1);
opts=glmnetSet;
opts.nlambda= nlambda;


min_err_iter= zeros(nalpha, 1);
lam_vect= zeros(nalpha, 1);
lam_ind =zeros(nalpha, 1);


for i=1: nalpha
    opts.alpha= alpha(i);
    
    if max(LF)>2
        d{i}= cvglmnet(FN,LF,'multinomial',opts,'class',folds_nr,fold_id);
    else
        d{i}= cvglmnet(FN,LF,'binomial',opts,'class',folds_nr,fold_id);
    end;
    
    %cvmMat{i=d{i}.cvm;
    err_now= min(d{i}.cvm)
    alp_now= alpha(i)
    min_err_iter(i)= err_now;
    lam_vect(i)= d{i}.lambda_min;
    [~, lam_ind(i)]= min(d{i}.lambda);
end;

[min_err, alp_ind]= min(min_err_iter);
%[min_err,alp_ind]=  min(m_col);


accuracy=(1-min_err)*100
std= d{alp_ind}.cvsd(lam_ind(alp_ind))*100
alp_opt= alpha(alp_ind)
lam_opt= d{alp_ind}.lambda(lam_ind(alp_ind))
sparsity= 100*(FSz-d{alp_ind}.nzero(lam_ind(alp_ind)))/FSz

pred= cvglmnetPredict(d{alp_ind},[],lam_opt,'nonzero');

if max(LF)>2
pred_Mat=cell2mat(pred);
else
pred_Mat= pred;
end;

pred_sum= sum(pred_Mat,2);

nz_ind= pred_sum>0;
FN_nz= FN(:,nz_ind);

min_err_iter_nz= zeros(nalpha, 1);
lam_vect_nz= zeros(nalpha, 1);
lam_ind_nz= zeros(nalpha, 1);
nalpha=40
alpha=linspace(0, 0.5, nalpha);
alpha=alpha(end:-1:1);
for i=1: nalpha
    opts.alpha= alpha(i);
    if max(LF)>2
        d_nz{i}= cvglmnet(FN_nz, LF, 'multinomial', opts, 'class', folds_nr, fold_id);
    else
        d_nz{i}= cvglmnet(FN_nz, LF, 'binomial', opts, 'class', folds_nr, fold_id);
    end;
    %cvmMat{i=d{i}.cvm;
    err_now= min(d_nz{i}.cvm)
    alp_now= alpha(i)
    min_err_iter_nz(i)= err_now;
    lam_vect_nz(i)= d_nz{i}.lambda_min;
    [~, lam_ind_nz(i)]= min(d_nz{i}.lambda);
end;
[min_err_nz, alp_ind_nz]= min(min_err_iter_nz);
%[min_err,alp_ind]=  min(m_col);


accuracy=(1-min_err_nz)*100
std= d_nz{alp_ind_nz}.cvsd(lam_ind_nz(alp_ind_nz))*100
alp_opt_nz= alpha(alp_ind_nz)
lam_opt_nz= d_nz{alp_ind_nz}.lambda(lam_ind_nz(alp_ind_nz))
%sparsity= 100*(FSz-d{alp_ind}.nzero(lam_ind(alp_ind)))/FSz

pred= cvglmnetPredict(d_nz{alp_ind_nz},[],lam_opt_nz,'nonzero');
b= glmnetPredict(d_nz{alp_ind_nz}.glmnet_fit,[],lam_opt_nz,'coefficients');
    
%end;
%%
%DSz=size(FN,1);
%fldsz=round(DSz/folds);
conf=[];
if split
    F=csvread('full_feature_set.csv');
    F=F(vals(indtest),1:1781);

    % find zeroed rows in F
    z_inds=max(abs(F).')==0;

    F(z_inds,:)=[];
    F=F(:,fvar.nzstd);
    dsz= size(F,1);
    v1= ones(dsz,1);
    FN=(F-v1*fvar.mu)./(v1*fvar.std);
    FN_nz=FN(:,nz_ind);
    
    conf(:,1)=glmval(b{1},FN_nz,'logit');
    conf(:,2)=glmval(b{2},FN_nz,'logit');
    conf(:,3)=glmval(b{3},FN_nz,'logit');
    conf(:,4)=glmval(b{4},FN_nz,'logit');
    conf(:,5)=glmval(b{5},FN_nz,'logit');
    [~, id]=max(conf,[],2);

    value=labels_test.value;
    value(value>5)=5;
    id=sortrows([id value.']
    confmat=confusionmat(id,value);
    figure; imagesc(confmat./repmat(sum(confmat,2),[1 5])); colorbar; caxis([0 1])
    
end

all_fvar{3}=fvar;
all_b{3}=b;
all_nz_ind{3}=nz_ind;
all_conf{3}=confmat;