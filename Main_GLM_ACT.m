%%
%Current Labels
%slip=1
%trip=2;
%rightfall=3;
%leftfall=4;
%sliplike=5;
%triplike=6;
%rightfalllike=7;
%leftfalllike=8;
%Activities=9; (everything else)

%Original (Momo's) labels for activities
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

split=1; %flag to split data into test and train sets (25-75) and create cofnusion matrix
class=0; % flag for fall classification (rather than detection only)
fall_like=9; % 5 to set as non-fall, 9 to set as fall
remove_falllike = 0;
remove_activities = 0;
java_feat=0;

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
nruns = 5;  %how many runs the classification model is trained/tested
% confmat_all = []; %save all confusion matrices when running multiple
% classification iters

%%% Loading the data

% load full data set
% any filtering of data past this point should be mirrored on the features
% loaded by extract_feature_phone_plus_labels


% for r = 1:nruns

load new_labels
% load labels_full.mat
% load labels_full_whome.mat


if remove_falllike
    inds = labels.value > 4 & labels.value < 9;
    labels.timestamp = labels.timestamp(~inds);
    labels.value = labels.value(~inds);
    labels.subject = labels.subject(~inds);
    labels.acce = labels.acce(~inds);
    labels.gyro = labels.gyro(~inds);
    labels.baro = labels.baro(~inds);
end

if remove_activities
    inds = labels.value == 9;
    labels.timestamp = labels.timestamp(~inds);
    labels.value = labels.value(~inds);
    labels.subject = labels.subject(~inds);
    labels.acce = labels.acce(~inds);
    labels.gyro = labels.gyro(~inds);
    labels.baro = labels.baro(~inds);
end


if split
    %Take 75% of original data for training
    inds = randperm(length(labels.acce));
    indtrain = inds(1:round(0.75*length(inds)));
    indtest = inds(round(0.75*length(inds))+1:end);
else
    indtrain = 1:length(labels.acce);
    indtest=[];
end


%
labels_test=labels;

%
%TEST DATA
% labels_test.timestamp = labels.timestamp(indtest);
labels_test.value = labels.value(indtest);
labels_test.subject = labels.subject(indtest);
labels_test.acce = labels.acce(indtest);
labels_test.gyro = labels.gyro(indtest);
labels_test.baro = labels.baro(indtest);

%TRAIN DATA
% labels.timestamp = labels.timestamp(indtrain);
labels.value = labels.value(indtrain);
labels.subject = labels.subject(indtrain);
labels.acce = labels.acce(indtrain);
labels.gyro = labels.gyro(indtrain);
labels.baro = labels.baro(indtrain);

%% Feature extraction

[FN, L, fold_id, muF, stdF, epsF, nzstd]=extract_feature_matlab(labels, folds_nr);
fvar.std= stdF;
fvar.mu= muF;
fvar.eps= epsF;
fvar.nzstd= nzstd;

FSz= size(FN,2); %feature size
DSz= size(FN,1); %dataset size

% Assign fall categories as 2 (falls and fall-like) or 1 (non-fall)
LB=(L<fall_like)+1;  %binary labeling

%fall classification labels
LF=L;
LF(L>4 & L<9)=5;

if ~class
    LF=LB;
end


alpha=linspace(0.5, 1, nalpha); %elastic net hyperparam
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
        d{i}= cvglmnet(FN,LF,'multinomial',opts,'class',folds_nr,fold_id); %fall classification
    else
        d{i}= cvglmnet(FN,LF,'binomial',opts,'class',folds_nr,fold_id); %fall detection model
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
    err_now= min(d_nz{i}.cvm) %mean cross validation error
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
%This part might be removed
% d= glmval(b, FN_nz, 'logit');
% bind=ceil(d-0.5)';
% L0=(LF-1);
% crate=(numel(L0)-sum(abs(bind-L0)))/numel(L0);


%if ~split
save class_params_ACT_nobar fvar b  nz_ind
%end

%end;
%% Testing the model on Test Data (if split is set to 1)
if split
    if java_feat
        javaaddpath('purple-robot-skyfall.jar');
        labels = labels_test;   %because java code expects labels structure
        save labels_struct_test labels
        com.company.TrainFeatureExtractor.extractFeatures('labels_struct_test.mat', 'test_features.csv');
        F=csvread('test_features.csv');
        F=F(:,1:1781);
        
        % find zeroed rows in F (feature cannot be calculated)
        z_inds=max(abs(F).')==0;

        F(z_inds,:)=[];
        F=F(:,fvar.nzstd);
        dsz= size(F,1);
        v1= ones(dsz,1);
        FN=(F-v1*fvar.mu)./(v1*fvar.std);
        FN_nz=FN(:,nz_ind);
    else
        F=[];
        for i=1:length(labels_test.acce)
            F=[F; extract_feature_test_phone(labels_test.acce{i}, labels_test.gyro{i}, labels_test.baro{i},fvar)];
        end
        FN_nz=F(:,nz_ind);
        z_inds=zeros(1,length(labels_test.value));
    end
    
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
        
        %to save confusion matrices from multiple runs
        %             if nruns > 1
        %                 confmat_all(:,:,r) = confmat;
        %             end
        figure; imagesc(confmat./repmat(sum(confmat,2),[1 5])); colorbar; caxis([0 1])
    else
        conf= glmval(b, FN_nz, 'logit');
        id= ceil(conf-0.5); 
        
        isfall = labels_test.value < 9;
        isfall = isfall(~z_inds);
        fall_err = sum(~id(isfall))/sum(isfall);
        nfall_err = sum(id(~isfall))/sum(~isfall);
        err = (sum(id ~= isfall))/length(id); %error rate
        confmat(:,:)=confusionmat(isfall,id==1);
        %plot confusion matrix
        activities = {'Non-Fall','Fall'};
        figure; imagesc(confmat./repmat(sum(confmat,2),[1 2]));
        confmat = confmat./repmat(sum(confmat,2),[1 2]);
        [cmin,cmax] = caxis;
        caxis([0,1]) %set colormap to 0 1
        ax = gca;
        ax.XTick = 1:size(activities,2);
        ax.YTick = 1:size(activities,2);
        set(gca,'XTickLabel',activities,'FontSize',14)
        set(gca,'YTickLabel',activities,'FontSize',14)
        ax.XTickLabelRotation = 45;
        axis square
        title('Fall Detection')
        %add text
        thres = 0.6;
        for i = 1:length(activities)
            for j = 1:length(activities)
                if confmat(i,j) < thres
                    col = [1 1 1];
                else
                    col = [0 0 0];
                end
                text(i-0.2,j,sprintf('%.3f',confmat(j,i)),'FontSize',12,'FontWeight','bold','Color',col);
            end
        end
        saveas(gcf,'./Figs/FallDetec.fig')
        saveas(gcf,'./Figs/FallDetec.jpg')

    end
    
end

%show # features used for accelerometer and gyro
load features_labels_full.mat
featuresLabels = featuresLabels(fvar.nzstd); %remove barometer features
featuresLabels=featuresLabels(nz_ind); %the features selected
Nacc = length(find(cellfun(@(x) strcmp(x(1:4),'acce'),featuresLabels)));
Ngyr = length(find(cellfun(@(x) strcmp(x(1:4),'gyro'),featuresLabels)));
figure, histogram([Nacc Ngyr]),        
set(gca,'XTickLabel',{'Accelerometer','Gyroscope'},'FontSize',14)
saveas(gcf,'./Figs/FeaturesUsed.jpg')
% end