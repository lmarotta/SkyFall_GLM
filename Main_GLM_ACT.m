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

skip_like=0; % flag to not use fall-like data
ACTnumber=500;
split=0; % flag to split data into test and train sets (25-75) and create cofnusion matrix

addpath(genpath('./glmnet_matlab/'))
% rng(10001)
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
% Can use either standard or jittered falls data
% load labels_plus_data.mat
load labels_plus_data_jittered_Unif

% load Activities data - used as non falls
ACT=load('labels_plus_data_ACT.mat');
ACTnumber=min([ACTnumber length(ACT.labels.subject)]);

%% Merge Falls and Activities Data

labels.value=[labels.value ACT.labels.value(1:ACTnumber)];
labels.subject=[labels.subject ACT.labels.subject(1:ACTnumber)];
labels.timestamp=[labels.timestamp ACT.labels.timestamp(1:ACTnumber)];
labels.text=[labels.text ACT.labels.text(1:ACTnumber)];
labels.acce=[labels.acce ACT.labels.acce(1:ACTnumber)];
labels.baro=[labels.baro ACT.labels.baro(1:ACTnumber)];
labels.gyro=[labels.gyro ACT.labels.gyro(1:ACTnumber)];
if skip_like
    nlike_inds=find(labels.value<5 | labels.value>8);

    labels.value=labels.value(nlike_inds);
    labels.subject=labels.subject(nlike_inds);
    labels.timestamp=labels.timestamp(nlike_inds);
    labels.text=labels.text(nlike_inds);
    labels.acce=labels.acce(nlike_inds);
    labels.baro=labels.baro(nlike_inds);
    labels.gyro=labels.gyro(nlike_inds);
end

if split
    %Take 75% of original data for training 
    inds = randperm(length(labels.timestamp));
    indtrain = inds(1:round(0.75*length(inds)));
    indtest = inds(round(0.75*length(inds))+1:end);
else
    indtrain = 1:length(labels.timestamp);
    indtest=[];
end

labels_test=labels;

%TEST DATA
labels_test.timestamp = labels.timestamp(indtest);
labels_test.value = labels.value(indtest);
labels_test.subject = labels.subject(indtest);
labels_test.text = labels.text(indtest);
labels_test.acce = labels.acce(indtest);
labels_test.gyro = labels.gyro(indtest);
labels_test.baro = labels.baro(indtest);

inds=cellfun(@(x,y,z) length(x)<100 | length(y)<100 | length(z)<10,labels_test.acce,labels_test.gyro,labels_test.baro);
inds=~inds;
labels_test.timestamp = labels_test.timestamp(inds);
labels_test.value = labels_test.value(inds);
labels_test.subject = labels_test.subject(inds);
labels_test.text = labels_test.text(inds);
labels_test.acce = labels_test.acce(inds);
labels_test.gyro = labels_test.gyro(inds);
labels_test.baro = labels_test.baro(inds);

% save indstest_act.mat indstest labels

%TRAIN DATA
labels.timestamp = labels.timestamp(indtrain);
labels.value = labels.value(indtrain);
labels.subject = labels.subject(indtrain);
labels.text = labels.text(indtrain);
labels.acce = labels.acce(indtrain);
labels.gyro = labels.gyro(indtrain);
labels.baro = labels.baro(indtrain);

%% Feature extraction

[FN, L, fold_id, muF, stdF, epsF, nzstd]=extract_feature_phone(labels, folds_nr);
fvar.std= stdF;
fvar.mu= muF;
fvar.eps= epsF;
fvar.nzstd= nzstd;

FSz= size(FN,2);
DSz= size(FN,1);

LB=(L<9)+1;
LF=L;
LF(L>4)=5;
%uncomment for fall detection
LF=LB;


alpha=linspace(0.5, 1, nalpha);
alpha=alpha(end:-1:1);
d=cell(nalpha, 1);
opts=glmnetSet;
opts.nlambda= nlambda;


min_err_iter= zeros(nalpha, 1);
lam_vect= zeros(nalpha, 1);
lam_ind =zeros(nalpha, 1);


for i=1: nalpha
    %rng(200);
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
    %rng(200);
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
d= glmval(b, FN_nz, 'logit');
bind=ceil(d-0.5);
L0=(LF-1);
crate=(numel(L0)-sum(abs(bind-L0)))/numel(L0);

Feature_Labels={'DCM',  'DCMed', 'DCFFT_re', 'DCFFT_im',  'DCFFT_abs', 'DCfit121', 'DCfit131',  'DCfit141', 'DCfit122', 'DCfit132',  'DCfit142', 'DCfit123', 'DCfit133', 'DCfit143', 'DCCorrV', 'DCstd', 'DCS', 'DCK', 'DDCM', 'DDCMed', 'DDCCorrV', 'DDCstd', 'DDCS', 'DDCK', 'DDCFFT_re', 'DDCFFT_im',  'DDCFFT_abs', 'DDCfit121', 'DDCfit131',  'DDCfit141', 'DDCfit122', 'DDCfit132',  'DDCfit142', 'DDCfit123', 'DDCfit133',  'DDCfit143' ...
    'DCM',  'DCMed', 'DCFFT_re', 'DCFFT_im',  'DCFFT_abs', 'DCfit121', 'DCfit131',  'DCfit141', 'DCfit122', 'DCfit132',  'DCfit142', 'DCfit123', 'DCfit133',  'DCfit143',  'DCCorrV', 'DCstd', 'DCS', 'DCK', 'DDCM', 'DDCMed', 'DDCCorrV', 'DDCstd', 'DDCS', 'DDCK', 'DDCFFT_re', 'DDCFFT_im',  'DDCFFT_abs', 'DDCfit121', 'DDCfit131',  'DDCfit141', 'DDCfit122', 'DDCfit132',  'DDCfit142', 'DDCfit123', 'DDCfit133',  'DDCfit143' ...
    'DCM', 'DCMed', 'DCCorrV', 'DCfit121',  'DCfit131', 'DCfit122',  'DCfit132',  'DCfit123', 'DCfit133',  'DCfit124',  'DCfit134', 'DCstd', 'DCskew', 'DCkurt', 'DDCM', 'DDCMed', 'DDCstd', 'DCFFT_re', 'DCFFT_im', 'DCFFT_abs'};

save class_params_ACT fvar b  nz_ind Feature_Labels

%end;
%% Genrating Randomn folds
%DSz=size(FN,1);
%fldsz=round(DSz/folds);
if split
    for k = 1:length(labels_test.acce)
        accel_data = labels_test.acce{k};
        gyro_data = labels_test.gyro{k};
        baro_data = labels_test.baro{k};
        [id(k), conf(k)]= fall_model_eval_ACT(accel_data, gyro_data,  baro_data);
    end

    isfall = labels_test.value < 9;
    fall_err = sum(~id(isfall))/sum(isfall);
    nfall_err = sum(id(~isfall))/sum(~isfall);
    err = (sum(id ~= isfall))/length(id); %error rate
    confmat(:,:)=confusionmat(isfall,id==1);
    figure; imagesc(confmat./repmat(sum(confmat,2),[1 2])); colorbar; caxis([0 1])
end