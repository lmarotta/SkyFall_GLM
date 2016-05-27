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

addpath ./glmnet_matlab/
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
load labels_plus_data.mat

%Take 75% of original data for training 
inds = randperm(length(labels.timestamp));
indtrain = inds(1:round(0.75*length(inds)));
indsjit = inds(round(0.75*length(inds))+1:end);
save indsjit.mat indsjit

%TRAIN DATA
labels.timestamp = labels.timestamp(indtrain);
labels.value = labels.value(indtrain);
labels.subject = labels.subject(indtrain);
labels.text = labels.text(indtrain);
labels.acce = labels.acce(indtrain);
labels.gyro = labels.gyro(indtrain);
labels.baro = labels.baro(indtrain);


%% Feature extraction

[FN, L, fold_id, muF, stdF, epsF, nzstd]=extract_feature_train_may(labels, folds_nr);
fvar.std= stdF;
fvar.mu= muF;
fvar.eps= epsF;
fvar.nzstd= nzstd;

FSz= size(FN,2);
DSz= size(FN,1);

LB=(L>4)+1;
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

save class_params fvar b  nz_ind

%end;
%% Testing the trained model on jittered test clips w different std dev

JitterClips(0) %create jittered version of the test clips 
%(0 = only use test clips, 1 = includes train clips and jitter them)
Test_Jittered  %test the current model with the jittered clips






