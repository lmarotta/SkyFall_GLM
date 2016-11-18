%Current Labels
%slip=1
%trip=2;
%rightfall=3;
%leftfall=4;
%Activities=9; (everything else)

function [confmat_all, conf_all, isfall_all, b_all, fvar_all, nz_ind_all]=Main_GLM_ACT()

epsF = 1e-6; %threshold on std dev for standardizing features
split=1; %flag to split data into test and train sets (25-75) and create confusion matrix
class=0; % flag for fall classification (rather than detection only)
no_baro=0; % 0 - use barometer
nbaro_features = 80;

% fall_like=9; % 5 to set as non-fall, 9 to set as fall
% remove_falllike = 0;
% remove_activities = 0;

addpath(genpath('./glmnet_matlab/'))
% rng(10001)
%% Initialization

nalpha=10;
nlambda=200;

%%% Loading the data

load Training_Data
%Feature matrix assembled as follows
% F = [subj_id location subjcode labels Features];


% if remove_falllike
%     inds = labels.value > 4 & labels.value < 9;
%     labels.timestamp = labels.timestamp(~inds);
%     labels.value = labels.value(~inds);
%     labels.subject = labels.subject(~inds);
%     labels.acce = labels.acce(~inds);
%     labels.gyro = labels.gyro(~inds);
%     labels.baro = labels.baro(~inds);
% end

% if remove_activities
%     inds = labels.value == 9;
%     labels.timestamp = labels.timestamp(~inds);
%     labels.value = labels.value(~inds);
%     labels.subject = labels.subject(~inds);
%     labels.acce = labels.acce(~inds);
%     labels.gyro = labels.gyro(~inds);
%     labels.baro = labels.baro(~inds);
% end


%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = F(:,4); %labels for classification
LB=(L<9);%+1;  %binary labeling
if ~class
    L=LB;
end
subjid = F(:,1:3);  %[subjid, location subjcode]
F = F(:,5:end); %the feature matrix

if split
    subj=unique(subjid(:,1));
    folds_nr = length(subj)-1; %for glmnet
else
    subj={'CF000'};
end

%LOSO CV
for indCV=1:length(subj)

    if split % LOOCV

        test_subj=subj(indCV);
        indtrain = subjid(:,1)~=test_subj;
        indtest = ~indtrain;

    else %train with all data
        indtrain = logical(ones(size(F,1),1));
        indtest=[];
    end   
    
    %
%     Ftest = F(indtest,5:end); Ftrain = F(indtrain,5:end);
    fold_id = subjid(indtrain,1); %subjtrain - to indicate the subject for each row in F
    fold_id(fold_id > test_subj) = fold_id(fold_id > test_subj)-1;
   
    %% Standardize (z-score) features
    
    fvar.std= std(F(indtrain,:));
    fvar.mu= mean(F(indtrain,:));
    fvar.eps= epsF;
    nzstd = (fvar.std./mean(fvar.std) > fvar.eps);
    fvar.nzstd= nzstd; %features w nonzero std dev

    FNZ = F(:,fvar.nzstd);
    FN = (FNZ - repmat(fvar.mu(fvar.nzstd),[size(FNZ,1),1])) ./ repmat(fvar.std(fvar.nzstd),[size(FNZ,1),1]); %features normalized
    
    if no_baro
        pre_baro=sum(nzstd(1:end-nbaro_features));
        fvar.nzstd(end-nbaro_features+1:end)=0;
        FN=FN(:,1:pre_baro);
    end
    
    FSz= size(FN,2); %feature size
    DSz= size(FN,1); %dataset size

%% Train GLMnet

    alpha=linspace(0.5, 1, nalpha); %elastic net hyperparam
    alpha=alpha(end:-1:1);
    d=cell(nalpha, 1);
    opts=glmnetSet;
    opts.nlambda= nlambda;
    min_err_iter= zeros(nalpha, 1);
    lam_min= zeros(nalpha, 1);
    lam_ind =zeros(nalpha, 1);

%loop over alpha - optimize on LOSOCV 
    for i=1: nalpha
        %rng(200);
        opts.alpha= alpha(i);

        if class
            d{i}= cvglmnet(FN(indtrain,:),L(indtrain),'multinomial',opts,'class',folds_nr,fold_id); %fall classification
        else
            d{i}= cvglmnet(FN(indtrain,:),L(indtrain),'binomial',opts,'class',folds_nr,fold_id); %fall detection model
        end;

        err_now= min(d{i}.cvm) %min cv error (cvm is the mean CV error, a vector of length(lambda))
        alp_now= alpha(i)      %current alpha value
        min_err_iter(i)= err_now;
        lam_min(i)= d{i}.lambda_min; %value of lambda that gives minimum cvm for current itearation i (alpha)
%         lam_ind(i)= find(d{i}.lambda==d{i}.lambda_min); %index of lambda that corresponding to min lambda?
        [~, lam_ind(i)]= min(d{i}.lambda);
    end

    [min_err, alp_ind]= min(min_err_iter);
    %[min_err,alp_ind]=  min(m_col);

    accuracy=(1-min_err)*100
    stdcvm= d{alp_ind}.cvsd(lam_ind(alp_ind))*100 %cvsd = estimate of standard error of cvm.
    alp_opt= alpha(alp_ind)
    lam_opt= d{alp_ind}.lambda_min
    sparsity= 100*(FSz-d{alp_ind}.nzero(lam_ind(alp_ind)))/FSz

    nz_ind= cvglmnetPredict(d{alp_ind},[],lam_opt,'nonzero');
    b= glmnetPredict(d{alp_ind}.glmnet_fit,[],lam_opt,'coefficients');
    b = b([true; nz_ind]);
   
    b_all{indCV}=b;
    fvar_all(indCV)=fvar;
    nz_ind_all{indCV}=nz_ind;

    if ~split
        save class_params_ACT_nobar fvar b nz_ind
    end

    %% Testing the model on Test Data (left out subject)
    if split
    
        FNtest = FN(indtest,nz_ind);

        if class
            conf(:,1)=glmval(b{1},FNtest,'logit');
            conf(:,2)=glmval(b{2},FNtest,'logit');
            conf(:,3)=glmval(b{3},FNtest,'logit');
            conf(:,4)=glmval(b{4},FNtest,'logit');
            conf(:,5)=glmval(b{5},FNtest,'logit');
            [~, pred]=max(conf,[],2);

            value=L(indtest);
%             value(value>5)=9; %nonfall=9
            confmat=confusionmat(pred,value);

            figure; imagesc(confmat./repmat(sum(confmat,2),[1 5])); colorbar; caxis([0 1])
        else
            conf= glmval(b, FNtest, 'logit');
            conf_all{indCV}=conf;
            pred= ceil(conf-0.5); 

            isfall = logical(L(indtest));
            isfall_all{indCV}=isfall;
            confmat(:,:)=confusionmat(isfall,pred==1);
            confmat_all(:,:,indCV)=confmat;
            
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

%     %show # features used for accelerometer and gyro
%     load features_labels_full.mat
%     featuresLabels = featuresLabels(fvar.nzstd); %remove barometer features
%     featuresLabels=featuresLabels(nz_ind); %the features selected
%     Nacc = length(find(cellfun(@(x) strcmp(x(1:4),'acce'),featuresLabels)));
%     Ngyr = length(find(cellfun(@(x) strcmp(x(1:4),'gyro'),featuresLabels)));
%     figure, histogram([Nacc Ngyr]),        
%     set(gca,'XTickLabel',{'Accelerometer','Gyroscope'},'FontSize',14)
%     saveas(gcf,'./Figs/FeaturesUsed.jpg')
%     % end
end

[X, Y, T]=perfcurve(isfall_all(2:end), conf_all(2:end), true,'XVals',[0:0.05:1]);
figure; errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));

isfall=cell2mat(isfall_all');
conf=cell2mat(conf_all');
[X, Y, T, AUC]=perfcurve(isfall, conf, true);
figure; plot(X,Y)
AUC