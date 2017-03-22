%Generate results for   the paper

% LOSO CV on Healthy - Train/Test on Indoor falls
% Input: Features set (0:reduced, 1:full)
%Number of locations (1:waist, 2:waist+pocket, 3:all)
% Output: AUC, Sens, Spec
function [Labresults,Homeresults] = GeneratePaperResults
close all

Locations={'Waist','Pocket','Hand','All'};
C=['r','y','c','m']; % colors for AUC bar plot
model=0; % 0 - RF, 1 - GLM
rng(200)

nData=0; %number of data points for training on 1 location only;
          %Set to 0 to remove constraint

features_used = ones(18,1); %full feature set
% features_used = zeros(18,1); features_used(8) = 1; %only magnitude features
% features_used([ 1 2 5 8 9]) = 1; %expanded set

featureset=getFeatureInds(features_used);

%% Store Amputees data (all 3 locations) for testing
%Feature matrix assembled as follows
% F = [subj_id location subjcode labels Features];
if ~exist('Test_Data_Amputees.mat','file')
    X_Amp = TrainingDataSetup([], [], 10, 1,'Test_Data_Amputees');
    labels_Amp = load('Test_Data_Amputees_labels');
    labels_Amp.acce = labels_Amp.labels.acce;
    labels_Amp.value = labels_Amp.labels.value;
else
    X_Amp = load('Test_Data_Amputees');
    X_Amp = X_Amp.F;
    
    labels_Amp = load('Test_Data_Amputees_labels');
    labels_Amp.acce = labels_Amp.labels.acce;
    labels_Amp.value = labels_Amp.labels.value;
end

if ~exist('HealthyData.mat','file')
    X = TrainingDataSetup([], [], 10, 0,'HealthyData');
    labels = load('HealthyData_labels');
    labels.acce = labels.labels.acce;
    labels.value = labels.labels.value;
else
    X = load('HealthyData');
    X = X.F;
    
    labels = load('HealthyData_labels');
    labels.acce = labels.labels.acce;
    labels.value = labels.labels.value;
end

if ~exist('OutdoorFalls.mat','file')
    O = TrainingDataSetup([], [], 10, 2, 'OutdoorFalls');
    labels_O = load('OutdoorFalls_labels');
    labels_O.acce = labels_O.labels.acce;
    labels_O.value = labels_O.labels.value;
else
    O = load('OutdoorFalls.mat');
    O = O.F;
    
    labels_O = load('OutdoorFalls_labels');
    labels_O.acce = labels_O.labels.acce;
    labels_O.value = labels_O.labels.value;
end

X = [X; O];
labels.acce=[labels.acce; labels_O.acce];
labels.value=[labels.value; labels_O.value];

% [AUC{2},Sens{2},Spec{2}] = LOSOCV(X,X_Amp,1,0,0); %RF
% [AUC{4},Sens{4},Spec{4}] = LOSOCV(X,X_Amp,2,0,0);
% [AUC{6},Sens{6},Spec{6}] = LOSOCV(X,X_Amp,3,0,0);
% [AUC{1},Sens{1},Spec{1}] = LOSOCV(X,X_Amp,1,1,0); %GLMnet
% [AUC{3},Sens{3},Spec{3}] = LOSOCV(X,X_Amp,2,1,0);

labROCfig=figure;
%Train and Test on all 3 locations
cvtype = [1 2 3]; %all cv
% cvtype = 2; %H-A only
[AUC,Sens,Spec,AUCErr,SpecCI,FPR,FNR,bootstat] = LOSOCV(X,X_Amp,1:3,1:3,nData,model,featureset,cvtype,0,labROCfig,labels,labels_Amp);
results.AUC = AUC;   %mean and SEM
results.AUCErr = AUCErr; %bootstrap CI
results.Sens = Sens; %mean and SEM
results.Spec = Spec;  %mean and SEM
results.SpecCI = SpecCI; %bootstrap CI @90% Sens
results.FPR = {FPR{1} FPR{2}(1:3) FPR{3}};
results.FNR = {FNR{1} FNR{2}(1:3) FNR{3}};
results.AUCboot = mean(bootstat(:,2)); %bootstrapped mean
results.Specboot = mean(bootstat(:,1));
results.mAUC = cellfun(@nanmean,AUC,'UniformOutput',false);


display('HA vs HH')
[h,p] = ttest2(results.AUC{2},results.AUC{1},'VarType','unequal')
display('HA vs AA')
[h,p] = ttest(results.AUC{2},results.AUC{3})

sprintf('\nH-H Mean AUC %.3f +- %.3f',nanmean(results.AUC{1}),1.96*nanstd(results.AUC{1})/sqrt(sum(~isnan(results.AUC{1}))))
sprintf('\nA-A Mean AUC %.3f +- %.3f',nanmean(results.AUC{3}),1.96*nanstd(results.AUC{3})/sqrt(sum(~isnan(results.AUC{3}))))
sprintf('\nH-A Mean AUC %.3f +- %.3f',nanmean(results.AUC{1}),1.96*nanstd(results.AUC{2})/sqrt(sum(~isnan(results.AUC{2}))))

sprintf('\nH-H Mean Sens %.3f +- %.3f',nanmean(results.Sens{1}),1.96*nanstd(results.Sens{1})/sqrt(sum(~isnan(results.Sens{1}))))
sprintf('\nA-A Mean Sens %.3f +- %.3f',nanmean(results.Sens{3}),1.96*nanstd(results.Sens{3})/sqrt(sum(~isnan(results.Sens{3}))))
sprintf('\nH-A Mean Sens %.3f +- %.3f',nanmean(results.Sens{1}),1.96*nanstd(results.Sens{2})/sqrt(sum(~isnan(results.Sens{2}))))

sprintf('\nH-H Mean SpeC %.3f +- %.3f',nanmean(results.Spec{1}),1.96*nanstd(results.Spec{1})/sqrt(sum(~isnan(results.Spec{1}))))
sprintf('\nA-A Mean Spec %.3f +- %.3f',nanmean(results.Spec{3}),1.96*nanstd(results.Spec{3})/sqrt(sum(~isnan(results.Spec{3}))))
sprintf('\nH-A Mean Spec %.3f +- %.3f',nanmean(results.Spec{1}),1.96*nanstd(results.Spec{2})/sqrt(sum(~isnan(results.Spec{2}))))


% %save figures
% for f = 1:4
%     figure(f)
%     ax = gca;
%     ax.FontSize = 16;
%     fig = gcf;
%     print(fig,['./Figs/Paper/Fig',num2str(f)],'-depsc','-r300')
% end

%% Train on 1 location and test on 3
cvtype = 2; %H-A only

% Test on waist
[wAUC,wSens,wSpec,AUCErr,SpecCI,FPR,FNR,bootstat] = LOSOCV(X,X_Amp,1:3,1,nData,model,featureset,cvtype,0,labROCfig,labels,labels_Amp);
results.waist.AUC = wAUC;
results.waist.AUCErr = AUCErr;
results.waist.Sens = wSens;
results.waist.Spec = wSpec;
results.waist.SpecCI = SpecCI;
results.waist.FPR = {FPR{1} FPR{2}(1:3) FPR{3}};
results.waist.FNR = {FNR{1} FNR{2}(1:3) FNR{3}};
results.waist.AUCboot = mean(bootstat(:,2)); %bootstrapped mean
results.waist.Specboot = mean(bootstat(:,1));
results.waist.mAUC = cellfun(@nanmean,wAUC,'UniformOutput',false);

% Test on pocket
[pAUC,pSens,pSpec,AUCErr,SpecCI,FPR,FNR,bootstat] = LOSOCV(X,X_Amp,1:3,2,nData,model,featureset,cvtype,0,labROCfig,labels,labels_Amp);
results.pock.AUC = pAUC;
results.pock.AUCErr = AUCErr;
results.pock.Sens = pSens;
results.pock.Spec = pSpec;
results.pock.SpecCI = SpecCI;
results.pock.FPR = {FPR{1} FPR{2}(1:3) FPR{3}};
results.pock.FNR = {FNR{1} FNR{2}(1:3) FNR{3}};
results.pock.AUCboot = mean(bootstat(:,2)); %bootstrapped mean
results.pock.Specboot = mean(bootstat(:,1));
results.pock.mAUC = cellfun(@nanmean,pAUC,'UniformOutput',false);

% Test on hand
[hAUC,hSens,hSpec,AUCErr,SpecCI,FPR,FNR,~] = LOSOCV(X,X_Amp,1:3,3,nData,model,featureset,cvtype,0,labROCfig,labels,labels_Amp);
results.hand.AUC = hAUC;
results.hand.AUCErr = AUCErr;
results.hand.Sens = hSens;
results.hand.Spec = hSpec;
results.hand.SpecCI = SpecCI;
results.hand.FPR = {FPR{1} FPR{2}(1:3) FPR{3}};
results.hand.FNR = {FNR{1} FNR{2}(1:3) FNR{3}};
results.hand.AUCboot = mean(bootstat(:,2)); %bootstrapped mean
results.hand.Specboot = mean(bootstat(:,1));
results.hand.mAUC = cellfun(@nanmean,hAUC,'UniformOutput',false);

Labresults=results;

display('Waist vs. All')
[h,p] = ttest(results.AUC{2},results.waist.AUC{2})
display('Pocket vs. All')
[h,p] = ttest(results.AUC{2},results.pock.AUC{2})
display('Hand vs. All')
[h,p] = ttest(results.AUC{2},results.hand.AUC{2})


%% Plot location results (Healthy-Amputee)
%need to add error bars
HAinds=~isnan(results.AUC{2});
HAlen=sum(HAinds);
figure, hold on
d=[results.waist.mAUC{2} results.pock.mAUC{2} results.hand.mAUC{2} results.mAUC{2}];
% for i=1:4
% bar(i,d(i),C(i))
% end
% figauc = errorbar(1:4,[results.waist.mAUC{2} results.pock.mAUC{2} results.hand.mAUC{2} results.mAUC{2}],...
%     [nanstd(results.waist.AUC{2})/sqrt(HAlen) nanstd(results.pock.AUC{2})/sqrt(HAlen) nanstd(results.hand.AUC{2})/sqrt(HAlen) nanstd(results.AUC{2})/sqrt(HAlen)],...
%     'linewidth',1.5,'linestyle','none','color','k');
% h = gca;
% h.YLim = [0.8 1];
% h.XTick = 1:4;
% h.XTickLabel = {'Waist','Pocket','Hand','All'};
% set(get(gca,'Xlabel'),'FontSize',16);
% set(gca,'FontSize',16);
% title('mean AUC')

% Add AUC values for ROC plots
labAUCstr=cell(4,1);
figure(labROCfig), hold on
for i=1:4
    subplot(2,2,i)
    labAUCstr{i}=sprintf('AUC = %0.3f',d(i));
    text(.4,.45,labAUCstr{i},'FontSize',24);
end

% [~,obj]=legend(labAUCstr);
% legend('boxoff')
% set(obj(2),'visible','off');

% plot FPR by location
FPR=[results.waist.FPR{2}; results.pock.FPR{2}; results.hand.FPR{2}; results.FPR{2}];

figure, hold on
imagesc(FPR);
M=max(max(FPR));
for i=1:3
    if results.waist.FPR{2}(i)<.2*M; Color='w'; else Color='k'; end
    text(i,1,sprintf('%0.2f',results.waist.FPR{2}(i)*100),'Color',Color)
    if results.pock.FPR{2}(i)<.2*M; Color='w'; else Color='k'; end
    text(i,2,sprintf('%0.2f',results.pock.FPR{2}(i)*100),'Color',Color)
    if results.hand.FPR{2}(i)<.2*M; Color='w'; else Color='k'; end
    text(i,3,sprintf('%0.2f',results.hand.FPR{2}(i)*100),'Color',Color)
    if results.FPR{2}(i)<.2*M; Color='w'; else Color='k'; end
    text(i,4,sprintf('%0.2f',results.FPR{2}(i)*100),'Color',Color)
end
set(gca,'YDir','reverse')
set(gca,'YTick',1:4)
set(gca,'XTick',1:4)
set(gca,'XTickLabel',{'Waist', 'Pocket', 'Hand'})
set(gca,'YTickLabel',{'Waist', 'Pocket', 'Hand', 'All'})
title('FPR')

FPRLab=table(FPR(:,1), FPR(:,2),FPR(:,3),'VariableNames',{'Waist','Pocket','Hand'},'RowNames',{'Waist','Pocket','Hand','All'});

% plot FNR by location
FNR=[results.waist.FNR{2}; results.pock.FNR{2}; results.hand.FNR{2}; results.FNR{2}];

figure, hold on
imagesc(FNR);
M=max(max(FNR));
for i=1:3
    if results.waist.FNR{2}(i)<.2*M; Color='w'; else Color='k'; end
    text(i,1,sprintf('%0.2f',results.waist.FNR{2}(i)*100),'Color',Color)
    if results.pock.FNR{2}(i)<.2*M; Color='w'; else Color='k'; end
    text(i,2,sprintf('%0.2f',results.pock.FNR{2}(i)*100),'Color',Color)
    if results.hand.FNR{2}(i)<.2*M; Color='w'; else Color='k'; end
    text(i,3,sprintf('%0.2f',results.hand.FNR{2}(i)*100),'Color',Color)
    if results.FNR{2}(i)<.2*M; Color='w'; else Color='k'; end
    text(i,4,sprintf('%0.2f',results.FNR{2}(i)*100),'Color',Color)
end
set(gca,'YDir','reverse')
set(gca,'YTick',1:4)
set(gca,'XTick',1:4)
set(gca,'XTickLabel',{'Waist', 'Pocket', 'Hand'})
set(gca,'YTickLabel',{'Waist', 'Pocket', 'Hand', 'All'})
title('FNR')

FNRLab=table(FNR(:,1),FNR(:,2),FNR(:,3),'VariableNames',{'Waist','Pocket','Hand'},'RowNames',{'Waist','Pocket','Hand','All'});

StatsLab=table({getMeanSEMStr(results.waist.AUC{2}); getMeanSEMStr(results.waist.Sens{2}); getMeanSEMStr(results.waist.Spec{2})},...
    {getMeanSEMStr(results.pock.AUC{2}); getMeanSEMStr(results.pock.Sens{2}); getMeanSEMStr(results.pock.Spec{2})},...
    {getMeanSEMStr(results.hand.AUC{2}); getMeanSEMStr(results.hand.Sens{2}); getMeanSEMStr(results.hand.Spec{2})},...
    {getMeanSEMStr(results.AUC{2}); getMeanSEMStr(results.Sens{2}); getMeanSEMStr(results.Spec{2})},...
    'VariableNames',{'Waist','Pocket','Hand','All'},'RowNames',{'AUC','Sensitivity','Specificity'});

StatsPop=table({getMeanSEMStr(results.AUC{1}); getMeanSEMStr(results.Sens{1}); getMeanSEMStr(results.Spec{1})},...
    {getMeanSEMStr(results.AUC{2}); getMeanSEMStr(results.Sens{2}); getMeanSEMStr(results.Spec{2})},...
    {getMeanSEMStr(results.AUC{3}); getMeanSEMStr(results.Sens{3}); getMeanSEMStr(results.Spec{3})},...
    'VariableNames',{'Healthy_Healthy', 'Healthy_Amputee', 'Amputee_Amputee'},...
    'RowNames',{'AUC','Sensitivity','Specificity'});

Labresults.FPRTable=FPRLab; Labresults.FNRTable=FNRLab; Labresults.StatsTable=StatsLab; Labresults.StatsPop=StatsPop;

writetable(Labresults.FPRTable,'./Figs/Paper/LabFPR.xlsx','WriteRowNames',true)
writetable(Labresults.FNRTable,'./Figs/Paper/LabFNR.xlsx','WriteRowNames',true)
writetable(Labresults.StatsTable,'./Figs/Paper/LocStats.xlsx','WriteRowNames',true)
writetable(Labresults.StatsPop,'./Figs/Paper/PopStats.xlsx','WriteRowNames',true)

%% HOME DATA ANALYSIS 
cvtype = 2;

filespath = 'Z:/Amputee Phones-R01/Home Data Collection/Amputees/';
l = load([filespath 'HomeDataAmp.mat']);
F = l.F;
sprintf('Data length = %.2f h',size(F,1)*5/60/60)

inds= any(bsxfun(@eq,X_Amp(:,1),unique(F(:,1))'),2) & X_Amp(:,4)<5; % Use falls from subjects fwith Home Data 

X_Amp = [X_Amp(inds,:);F(randperm(size(F,1),500),:)];

homeROCfig=figure;

% % Test on waist
% [wAUC,wSens,wSpec,AUCErr,SpecCI,FPR,FNR,bootstat] = LOSOCV(X,X_Amp,1:3,1,nData,model,featureset,cvtype,0,homeROCfig);
% results.waist.AUC = wAUC;
% results.waist.AUCErr = AUCErr;
% results.waist.Sens = wSens;
% results.waist.Spec = wSpec;
% results.waist.SpecCI = SpecCI;
% results.waist.FPR = {FPR{1} FPR{2}(4) FPR{3}};
% results.waist.FNR = {FNR{1} FNR{2}(4) FNR{3}};
% results.waist.AUCboot = mean(bootstat(:,2)); %bootstrapped mean
% results.waist.Specboot = mean(bootstat(:,1));
% results.waist.bootstat = bootstat;
% 
% % Test on pocket
% [pAUC,pSens,pSpec,AUCErr,SpecCI,FPR,FNR,bootstat] = LOSOCV(X,X_Amp,1:3,2,nData,model,featureset,cvtype,0,homeROCfig);
% results.pock.AUC = pAUC;
% results.pock.AUCErr = AUCErr;
% results.pock.Sens = pSens;
% results.pock.Spec = pSpec;
% results.pock.SpecCI = SpecCI;
% results.pock.FPR = {FPR{1} FPR{2}(4) FPR{3}};
% results.pock.FNR = {FNR{1} FNR{2}(4) FNR{3}};
% results.pock.AUCboot = mean(bootstat(:,2)); %bootstrapped mean
% results.pock.Specboot = mean(bootstat(:,1));
% results.pock.bootstat = bootstat;
% 
% % Test on hand
% [hAUC,hSens,hSpec,AUCErr,SpecCI,FPR,FNR,bootstat] = LOSOCV(X,X_Amp,1:3,3,nData,model,featureset,cvtype,0,homeROCfig);
% results.hand.AUC = hAUC;
% results.hand.AUCErr = AUCErr;
% results.hand.Sens = hSens;
% results.hand.Spec = hSpec;
% results.hand.SpecCI = SpecCI;
% results.hand.FPR = {FPR{1} FPR{2}(4) FPR{3}};
% results.hand.FNR = {FNR{1} FNR{2}(4) FNR{3}};
% results.hand.AUCboot = mean(bootstat(:,2)); %bootstrapped mean
% results.hand.Specboot = mean(bootstat(:,1));
% results.hand.bootstat = bootstat;

% 3 Locations
[AUC,Sens,Spec,AUCErr,SpecCI,FPR,FNR,bootstat] = LOSOCV(X,X_Amp,1:3,1:3,nData,model,featureset,cvtype,1,homeROCfig,labels,labels_Amp);
results.AUC = AUC;
results.AUCErr = AUCErr;
results.Sens = Sens;
results.Spec = Spec;
results.SpecCI = SpecCI;
results.FPR = {FPR{1} FPR{2}(4) FPR{3}};
results.FNR = {FNR{1} FNR{2}(4) FNR{3}};
results.AUCboot = mean(bootstat(:,2)); %bootstrapped mean
results.Specboot = mean(bootstat(:,1));
results.bootstat = bootstat;

Homeresults = results;
%% Plot location results (Healthy-Amputee)
%need to add error bars
figure, hold on
d=[results.waist.AUCboot results.pock.AUCboot results.hand.AUCboot results.AUCboot];
% for i=1:4
% bar(i,d(i),C(i))
% end
% figauc = errorbar(1:4,[results.waist.AUCboot results.pock.AUCboot results.hand.AUCboot results.AUCboot],...
%     [results.waist.AUCboot-results.waist.AUCErr{cvtype}(1) results.pock.AUCboot-results.pock.AUCErr{cvtype}(1) results.hand.AUCboot-results.hand.AUCErr{cvtype}(1) results.AUCboot-results.AUCErr{cvtype}(1)],...
%     [results.waist.AUCErr{cvtype}(2)-results.waist.AUCboot results.pock.AUCErr{cvtype}(2)-results.pock.AUCboot results.hand.AUCErr{cvtype}(2)-results.hand.AUCboot results.AUCErr{cvtype}(2)-results.AUCboot],...
%     'linewidth',1.5,'linestyle','none','color','k');h = gca;
% h.YLim = [0.8 1];
% h.XTick = 1:4;
% h.XTickLabel = {'Waist','Pocket','Hand','All'};
% set(get(gca,'Xlabel'),'FontSize',16);
% set(gca,'FontSize',16);
% title('mean AUC')

% Add AUC in legend for ROC plots
homeAUCstr=cell(4,1);
figure(homeROCfig), hold on
homeAUCstr{4}=sprintf('AUC = %0.3f',d(4));
text(.4,.45,homeAUCstr{4},'FontSize',24);

% % plot FPR by location
% FPR=[results.waist.FPR{2}; results.pock.FPR{2}; results.hand.FPR{2}; results.FPR{2}];
% 
% figure, hold on
% imagesc(FPR);
% M=max(max(FPR));
% i = 1;
% if results.waist.FPR{2}(i)<.2*M; Color='w'; else Color='k'; end
% text(i,1,sprintf('%0.2f',results.waist.FPR{2}(i)*100),'Color',Color)
% if results.pock.FPR{2}(i)<.2*M; Color='w'; else Color='k'; end
% text(i,2,sprintf('%0.2f',results.pock.FPR{2}(i)*100),'Color',Color)
% if results.hand.FPR{2}(i)<.2*M; Color='w'; else Color='k'; end
% text(i,3,sprintf('%0.2f',results.hand.FPR{2}(i)*100),'Color',Color)
% if results.FPR{2}(i)<.2*M; Color='w'; else Color='k'; end
% text(i,4,sprintf('%0.2f',results.FPR{2}(i)*100),'Color',Color)
% 
% set(gca,'YDir','reverse')
% set(gca,'YTick',1:4)
% set(gca,'XTick',1:4)
% set(gca,'XTickLabel',{'Waist', 'Pocket', 'Hand'})
% set(gca,'YTickLabel',{'Waist', 'Pocket', 'Hand', 'All'})
% title('FPR')
% 
% FPRHome=table(FPR(:,1),'VariableNames',{'All'},'RowNames',{'Waist','Pocket','Hand','All'});
% 
% % plot FNR by location
% FNR=[results.waist.FNR{2}; results.pock.FNR{2}; results.hand.FNR{2}; results.FNR{2}];
% 
% figure, hold on
% imagesc(FNR);
% M=max(max(FNR));
% if results.waist.FNR{2}(i)<.2*M; Color='w'; else Color='k'; end
% text(i,1,sprintf('%0.2f',results.waist.FNR{2}(i)*100),'Color',Color)
% if results.pock.FNR{2}(i)<.2*M; Color='w'; else Color='k'; end
% text(i,2,sprintf('%0.2f',results.pock.FNR{2}(i)*100),'Color',Color)
% if results.hand.FNR{2}(i)<.2*M; Color='w'; else Color='k'; end
% text(i,3,sprintf('%0.2f',results.hand.FNR{2}(i)*100),'Color',Color)
% if results.FNR{2}(i)<.2*M; Color='w'; else Color='k'; end
% text(i,4,sprintf('%0.2f',results.FNR{2}(i)*100),'Color',Color)
% set(gca,'YDir','reverse')
% set(gca,'YTick',1:4)
% set(gca,'XTick',1:4)
% set(gca,'XTickLabel',{'Waist', 'Pocket', 'Hand'})
% set(gca,'YTickLabel',{'Waist', 'Pocket', 'Hand', 'All'})
% title('FNR')
% 
% FNRHome=table(FNR(:,1),'VariableNames',{'All'},'RowNames',{'Waist','Pocket','Hand','All'});
% 
% StatsHome=table(getMeanCIStr(results.waist.bootstat,results.waist.AUCErr{2},results.waist.SpecCI{2}),...
%     getMeanCIStr(results.pock.bootstat,results.pock.AUCErr{2},results.pock.SpecCI{2}),...
%     getMeanCIStr(results.hand.bootstat,results.hand.AUCErr{2},results.hand.SpecCI{2}),...
%     getMeanCIStr(results.bootstat,results.AUCErr{2},results.SpecCI{2}),...
%     'VariableNames',{'Waist','Pocket','Hand','All'},'RowNames',{'AUC','Sensitivity','Specificity'});
% 
% Homeresults.FPRTable=FPRHome; Homeresults.FNRTable=FNRHome; Homeresults.StatsTable=StatsHome;
% writetable(Homeresults.FPRTable,'./Figs/Paper/HomeFPR.xlsx','WriteRowNames',true)
% writetable(Homeresults.FNRTable,'./Figs/Paper/HomeFNR.xlsx','WriteRowNames',true)
% writetable(Homeresults.StatsTable,'./Figs/Paper/HomeStats.xlsx','WriteRowNames',true)

end



%cvtype is a vector with the cases
%[1 2 3] = H-H, H-A, A-A
function [AUC,Sens,Spec,AUCErr,SpecCI,FPR,FNR,bootstat] = LOSOCV(X,X_test,locations_train,locations_test,nData,model,featureset,cvtype,HomeData,ROCfig,HLabels,ALabels)

rng(400)

Locations={'Waist','Pocket','Hand','All'};

% if HomeData
%     l=load('Z:/Amputee Phones-R01/Home Data Collection/Healthy/HomeDataHealthy.mat');
%     HomeData=l.F;
%     HomeData = HomeData(:,featureset);
%     HomeL = zeros(size(HomeData,1),1); %all labels are non-fall
% end

%convert labels to binary (1,4=falls, 9=activities/nonfalls)
% Assign fall categories as 1 (falls) or 0 (non-fall)
L = X(:,4); %labels for classification
L=(L<9);%+1;  %binary labeling
subjid = X(:,1:3);  %[subjid, location subjcode]
subj=unique(subjid(:,1));
F=X(:,5:end); %remove metadata
F = F(:,featureset);

L_w = L(any(bsxfun(@eq,subjid(:,2),locations_train),2));
ratio = sum(L_w)/sum(~L_w);   %falls/non_falls
sprintf('Nlocs = %d, falls/non_falls = %f',length(locations_train),ratio)

conf_all=cell(1,length(subj));
isfall_all=cell(1,length(subj));
confmat_all=zeros(2,2,length(subj));
Nruns = 1;  %average multiple runs for each subject
alpha = 0.6;
lambda = 0.015; % default .015

%% LOSO CV HEALTHY
AUC_HH=zeros(Nruns,length(subj));
Sens_HH=zeros(Nruns,length(subj));
Spec_HH=zeros(Nruns,length(subj));

if length(cvtype)==3
    popFig=figure;
end
if find(cvtype==1)
    
    for indCV=1:length(subj)
        
        
        test_subj=subj(indCV);
        indtrain = subjid(:,1)~=test_subj & any(bsxfun(@eq,subjid(:,2),locations_train),2);
        indtest = subjid(:,1) == test_subj & any(bsxfun(@eq,subjid(:,2),locations_test),2);
        
        falls_ind = find(indtrain & L);
        nfalls_ind = find(indtrain & ~L);
        
        %balance the dataset - randomly sample nData instances from each class
        %repeat it 5 times and average results
        for run = 1:Nruns
            if nData > 0
                indtrain = [falls_ind(randperm(length(falls_ind),nData)); nfalls_ind(randperm(length(nfalls_ind),nData))];
            else
                indtrain = [falls_ind; nfalls_ind];
            end
                
            if model
                % Train GLMnet                
                
                [fvar_si,b,nz_ind]=Modeltrain(F(indtrain,:),L(indtrain),alpha,lambda,0);
                % Testing the model on Test Data (left out subject)
                [pred,conf,confmat] = Modeleval(F(indtest,:),L(indtest),fvar_si,nz_ind,b,0.5,0);
                conf_all{indCV}=conf;
                confmat_all(:,:,indCV)=confmat;
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                [~, ~, ~, AUC_HH(run,indCV)]=perfcurve(isfall, conf, true);
                Sens_HH(run,indCV) = sum(pred & isfall)/sum(isfall);
                Spec_HH(run,indCV) = sum(~pred & ~isfall)/sum(~isfall);
            else
                RFModel=TreeBagger(100,F(indtrain,:),L(indtrain));
                [pred,conf]=predict(RFModel,F(indtest,:));
                conf = conf(:,2); %posterior prob of a fall
                conf_all{indCV}=conf;
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                [~, ~, ~, AUC_HH(run,indCV)]=perfcurve(isfall, conf, true);
                
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                pred=cellfun(@(x) logical(str2double(x)),pred);
                
                confmat=confusionmat(isfall,pred);
                confmat_all(:,:,indCV)=confmat;
                
                Sens_HH(run,indCV) = sum(pred & isfall)/sum(isfall);
                Spec_HH(run,indCV) = sum(~pred & ~isfall)/sum(~isfall);
                
            end
            
            Train.acce=HLabels.acce(indtrain);
            Train.value=HLabels.value(indtrain);

            Test.acce=HLabels.acce(indtest);
            Test.value=HLabels.value(indtest);

            [FPR_Thresh(indCV), FNR_Thresh(indCV)] = ThresholdDetection(Train,Test);
            
        end
    end
    
    confmat=sum(confmat_all,3);
    [X, Y, T, AUC]=perfcurve(isfall_all, conf_all, true,'TVals',[0:0.05:1]); %conf bounds with CV
        %     [X, Y, T, AUC]=perfcurve(cell2mat(isfall_all'), cell2mat(conf_all'), true,'Nboot',1000,'XVals',[0:0.05:1]); %cb with Bootstrap

    if length(cvtype) == 3
        plotConfmat(confmat,'Healthy-Healthy',1)
        figure(popFig), subplot(2,2,4), hold on
        e = errorbar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
        e.LineWidth = 2; e.Marker = 'o';
        xlabel('False positive rate')
        ylabel('True positive rate')
        
        plot(nanmean(FPR_Thresh),nanmean(1-FNR_Thresh),'x')
    else
        plotConfmat(confmat,'Healthy-Healthy');
        %plot ROC curves w confidence bounds
        figroc = figure;
        e = errorbar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
        e.LineWidth = 2; e.Marker = 'o';
        xlabel('False positive rate')
        ylabel('True positive rate')
    end

end
%% Train on Healthy, Test on Amputees
AUC_HA=NaN(Nruns,length(unique(X_test(:,1))),'double');
Sens_HA=NaN(Nruns,length(unique(X_test(:,1))),'double');
Spec_HA=NaN(Nruns,length(unique(X_test(:,1))),'double');

if any(cvtype == 2)
    
    isfall_all = {};
    conf_all = {};
    Ft = X_test(:,5:end);
    Ft = Ft(:,featureset);
    
    indtrain = any(bsxfun(@eq,subjid(:,2),locations_train),2);
    
    falls_ind = find(indtrain & L);
    nfalls_ind = find(indtrain & ~L);
    
    for run = 1:Nruns
        if nData > 0
            indtrain = [falls_ind(randperm(length(falls_ind),nData)); nfalls_ind(randperm(length(nfalls_ind),nData))];
        else
            indtrain = [falls_ind; nfalls_ind];
        end
        
        if model
            [fvar_si,b,nz_ind]=Modeltrain(F(indtrain,:),L(indtrain),alpha,lambda,0);
            
            % Testing the model on each amputee subj (external set)
            locsdata = X_test(:,2);
            S=unique(X_test(:,1));
            for i = 1:length(unique(X_test(:,1)))
                subj=S(i);
                rowid = (X_test(:,1)==subj) & any(bsxfun(@eq,X_test(:,2),locations_test),2);
                rowid_all{i}=find(rowid);
                [pred,conf,confmat] = Modeleval(Ft(rowid,:),X_test(rowid,4)<5,fvar_si,nz_ind,b,0.5,0);
                conf_all{i}=conf;
                confmat_all(:,:,subj)=confmat;
                isfall = X_test(rowid,4)<5;
                isfall_all{i}=isfall;
                if length(unique(isfall)) <2 %subject AF06 with missing activities - TO FIX
                    isfall_all{i} = [];  conf_all{i}=[]; rowid_all{i}=[];
                    locsdata(rowid) = 0;
                else
                    [~, ~, ~, AUC_HA(run,i)]=perfcurve(isfall, conf, true);
                    Sens_HA(run,i) = sum(pred & isfall)/sum(isfall);
                    Spec_HA(run,i) = sum(~pred & ~isfall)/sum(~isfall);
                end
                
                Train.acce=HLabels.acce(indtrain);
                Train.value=HLabels.value(indtrain);

                Test.acce=ALabels.acce(rowid);
                Test.value=ALabels.value(rowid);

                [FPR_Thresh(i), FNR_Thresh(i)] = ThresholdDetection(Train,Test);
                
            end
        else
            RFModel=TreeBagger(100,F(indtrain,:),L(indtrain));
            
            % Testing the model on each amputee subj (external set)
            locsdata = X_test(:,2);
            S=unique(X_test(:,1));
            
            for i = 1:length(unique(X_test(:,1)))
                subj=S(i);
                rowid = (X_test(:,1)==subj) & any(bsxfun(@eq,X_test(:,2),locations_test),2);
                rowid_all{i}=find(rowid);
                [pred,conf]=predict(RFModel,Ft(rowid,:));
                conf = conf(:,2); %posterior prob of a fall
                conf_all{i}=conf;
                isfall = X_test(rowid,4)<5;
                isfall_all{i}=isfall;
                
                pred=cellfun(@(x) logical(str2double(x)),pred);
                
                if length(unique(isfall)) <2 %subject AF06 with missing activities - TO FIX
                    isfall_all{i} = [];  conf_all{i}=[]; rowid_all{subj}=[];
                    locsdata(rowid) = 0;
                else
                    [~, ~, ~, AUC_HA(run,i)]=perfcurve(isfall, conf, true);
                    Sens_HA(run,i) = sum(pred & isfall)/sum(isfall);
                    Spec_HA(run,i) = sum(~pred & ~isfall)/sum(~isfall);
                    confmat=confusionmat(X_test(rowid,4)<9,logical(pred));
                    confmat_all(:,:,i)=confmat;
                end
                  
                Train.acce=HLabels.acce(indtrain);
                Train.value=HLabels.value(indtrain);

                Test.acce=ALabels.acce(rowid);
                Test.value=ALabels.value(rowid);

                [FPR_Thresh(i), FNR_Thresh(i)] = ThresholdDetection(Train,Test);
                
            end
            
        end
    end
    
    confmat=sum(confmat_all,3);   
    
    rowid_all = rowid_all(~cellfun(@isempty,isfall_all));
    conf_all = conf_all(~cellfun(@isempty,isfall_all));
    isfall_all = isfall_all(~cellfun(@isempty,isfall_all)); %remove empty cell from missing subj
        
    if length(unique(X_test(:,1)))>4 % Do bootstrapping for smaller number of subjects
        [X, Y, T, AUC]=perfcurve(isfall_all, conf_all, true,'TVals',[0:0.05:1]);
    else
        [X, Y, T, AUC]=perfcurve(cell2mat(isfall_all'), cell2mat(conf_all'), true,'Nboot',1000,'TVals',[0:0.05:1]); %cb with Bootstrap
    end
        
    %plot cmat and ROC
    if length(cvtype) == 3
        plotConfmat(confmat,'Healthy-Amputee',3)
        figure(popFig), subplot(2,2,4), hold on,
        e = errorbar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
        e.LineWidth = 2; e.Marker = 'o';
        xlabel('False positive rate')
        ylabel('True positive rate')
        
        plot(nanmean(FPR_Thresh),nanmean(1-FNR_Thresh),'x')
        
        figure(ROCfig), hold on
        subplot(2,2,max(locations_train)+1)
        e = errorbar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
        e.LineWidth = 2; e.Marker = 'o';
        xlabel('False positive rate')
        ylabel('True positive rate')
        set(gca,'FontSize',16);
        title(Locations{4});
        ylim([-0.1 1.1]), xlim([-0.1 1.1])
    else
        plotConfmat(confmat,'Healthy-Amputee');
        %plot ROC curves w confidence bounds    
%         if ~exist('figroc','var')
%             figroc = figure;
%         end
%         figure(figroc), hold on  
%         e = errorbar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
%         e.LineWidth = 2; e.Marker = 'o';
%         xlabel('False positive rate')
%         ylabel('True positive rate')
        
        figure(ROCfig), hold on
        if ~HomeData
            if length(locations_test)==1
                subplot(2,2,locations_test(1))
            else
                subplot(2,2,4)
            end
            e = errorbar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
            e.LineWidth = 2; e.Marker = 'o';
            xlabel('False positive rate')
            ylabel('True positive rate')
            set(gca,'FontSize',16);
            ylim([-0.1 1.1]), xlim([-0.1 1.1])
            if length(locations_test)==1
                title(Locations{locations_test(1)});
            else
                title(Locations{4});
            end
            
            plot(nanmean(FPR_Thresh),nanmean(1-FNR_Thresh),'x')
        else
            e = errorbar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
            e.LineWidth = 2; e.Marker = 'o';
            xlabel('False positive rate')
            ylabel('True positive rate')
            set(gca,'FontSize',16);
            ylim([-0.1 1.1]), xlim([-0.1 1.1])
            title('Lab+Home Results')
            
            plot(nanmean(FPR_Thresh),nanmean(1-FNR_Thresh),'x')
        end
        
    end

    sens = 0.9;
    %95% CI on specificity for fixed sensitivity
    Handle=@(x) specAUC(x,sens);
    TL = [cell2mat(conf_all') cell2mat(cellfun(@double,isfall_all','UniformOutput',false))];
    bootstat = bootstrp(1000,Handle,TL,'Options',statset('UseParallel',true));
    ci = bootci(1000,{Handle,TL},'Options',statset('UseParallel',true));
    ci_spec = ci(:,[1 3:4]); % Spec@90 Sens Spec
    ci_AUC = ci(:,2); % AUC
    
%     Spec_HA = mean(bootstat(:,1)); %the spec at 90% sensitivity
%     AUC_HA = mean(bootstat(:,2));  %mean AUC
    
    e.LineWidth = 2; e.Marker = 'o';
    AUCerr_HA=ci_AUC;
    SpecCI_HA = ci_spec;

    %FPR and FNR by location
    locsdata=locsdata(cell2mat(rowid_all'));
    indFN = (cell2mat(conf_all')<.5 & cell2mat(isfall_all')==1);
    indFP = (cell2mat(conf_all')>=.5 & cell2mat(isfall_all')==0);
    locsFP = locsdata(indFP);
    locsFN = locsdata(indFN);
    hFP=sum(locsFP==3); wFP=sum(locsFP==1); pFP = sum(locsFP==2);
    hFN=sum(locsFN==3); wFN=sum(locsFN==1); pFN = sum(locsFN==2);
%     figure, bar([wFP/sum(cell2mat(isfall_all') & locsdata==1) wFN/sum(~cell2mat(isfall_all') & locsdata==1); ...
%         pFP/sum(cell2mat(isfall_all') & locsdata==2) pFN/sum(~cell2mat(isfall_all') & locsdata==2); ... 
%         hFP/sum(cell2mat(isfall_all') & locsdata==3) hFN/sum(~cell2mat(isfall_all') & locsdata==3)]), legend('FP','FN')
    
    FPR_HA=[wFP/sum(~cell2mat(isfall_all') & locsdata==1) pFP/sum(~cell2mat(isfall_all') & locsdata==2) hFP/sum(~cell2mat(isfall_all') & locsdata==3) length(locsFP)/sum(~cell2mat(isfall_all'))];
    FNR_HA=[wFN/sum(cell2mat(isfall_all') & locsdata==1) pFN/sum(cell2mat(isfall_all') & locsdata==2) hFN/sum(cell2mat(isfall_all') & locsdata==3) length(locsFN)/sum(cell2mat(isfall_all'))];
end

%% LOSO on Amputees

AUC_AA=NaN(Nruns,length(subj),'double');
Sens_AA=NaN(Nruns,length(subj),'double');
Spec_AA=NaN(Nruns,length(subj),'double');
if any(cvtype == 3)
    
    
    isfall_all = {};
    conf_all = {};
    
    F = X_test(:,5:end);
    F = F(:,featureset);
    L = X_test(:,4); %labels for classification
    L=(L<9);%+1;  %binary labeling
    subjid = X_test(:,1:3);  %[subjid, location subjcode]
    subj=unique(subjid(:,1));
    
    for indCV=1:length(subj)
        
        test_subj=subj(indCV);
        indtrain = subjid(:,1)~=test_subj & any(bsxfun(@eq,subjid(:,2),locations_train),2);
        indtest = subjid(:,1) == test_subj & any(bsxfun(@eq,subjid(:,2),locations_test),2);
        
        falls_ind = find(indtrain & L);
        nfalls_ind = find(indtrain & ~L);
        
        for run = 1:Nruns
            %balance the dataset
            if nData > 0
                indtrain = [falls_ind(randperm(length(falls_ind),nData)); nfalls_ind(randperm(length(nfalls_ind),nData))];
            else
                indtrain = [falls_ind; nfalls_ind];
            end
            
            if model
                % Train GLMnet
                
                [fvar_si,b,nz_ind]=Modeltrain(F(indtrain,:),L(indtrain),alpha,lambda,0);
                % Testing the model on Test Data (left out subject)
                [pred,conf,confmat] = Modeleval(F(indtest,:),L(indtest),fvar_si,nz_ind,b,0.5,0);
                conf_all{indCV}=conf;
                confmat_all(:,:,indCV)=confmat;
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                if length(unique(isfall)) <2 %subject AF06 with missing activities - TO FIX
                    isfall_all{indCV} = [];  conf_all{indCV}=[];
                else
                    [~, ~, ~, AUC_AA(run,indCV)]=perfcurve(isfall, conf, true);
                    Sens_AA(run,indCV) = sum(pred & isfall)/sum(isfall);
                    Spec_AA(run,indCV) = sum(~pred & ~isfall)/sum(~isfall);
                end
            else
                RFModel=TreeBagger(100,F(indtrain,:),L(indtrain));
                [pred,conf]=predict(RFModel,F(indtest,:));
                conf = conf(:,2); %posterior prob of a fall
                conf_all{indCV}=conf;
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                
                
                isfall = logical(L(indtest));
                isfall_all{indCV}=isfall;
                pred=cellfun(@(x) logical(str2double(x)),pred);
                
                confmat=confusionmat(isfall,pred);
                confmat_all(:,:,indCV)=confmat;
                
                if length(unique(isfall)) <2 %subject AF06 with missing activities - TO FIX
                    isfall_all{indCV} = [];  conf_all{indCV}=[];
                else
                    [~, ~, ~, AUC_AA(run,indCV)]=perfcurve(isfall, conf, true);
                    Sens_AA(run,indCV) = sum(pred & isfall)/sum(isfall);
                    Spec_AA(run,indCV) = sum(~pred & ~isfall)/sum(~isfall);
                end
            end
            
            Train.acce=ALabels.acce(indtrain);
            Train.value=ALabels.value(indtrain);

            Test.acce=ALabels.acce(indtest);
            Test.value=ALabels.value(indtest);

            [FPR_Thresh(indCV), FNR_Thresh(indCV)] = ThresholdDetection(Train,Test);
            
        end
    end
    
    confmat=sum(confmat_all,3);
    [X, Y, T, AUC]=perfcurve(isfall_all(~cellfun(@isempty,isfall_all)), conf_all(~cellfun(@isempty,isfall_all)), true,'TVals',[0:0.05:1]);
%     [X, Y, T, AUC]=perfcurve(cell2mat(isfall_all'),cell2mat(conf_all'),true,'Nboot',1000,'XVals',[0:0.05:1]);

    %plot cmat and ROC
    if length(cvtype) == 3
        figure(popFig)
        plotConfmat(confmat,'Amputee-Amputee',2)
        figure(popFig), subplot(2,2,4),  hold on
        e = errorbar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
        e.LineWidth = 2; e.Marker = 'o';
        xlabel('False positive rate')
        ylabel('True positive rate')
        legend('Healthy-Healthy','Healthy-Amputee','Amputee-Amputee')
        xlim([-0.1 0.6]), ylim([0 1.05]), axis square
        
        plot(nanmean(FPR_Thresh),nanmean(1-FNR_Thresh),'x')
    else
        plotConfmat(confmat,'Amputee-Amputee');
        %plot ROC curves w confidence bounds    
        if ~exist('figroc','var')
            figroc = figure;
        end
        figure(figroc), hold on  
        e = errorbar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1));
        e.LineWidth = 2; e.Marker = 'o';
        xlabel('False positive rate')
        ylabel('True positive rate')
        legend('Healthy-Healthy','Healthy-Amputee','Amputee-Amputee')
    end

end


AUC = {AUC_HH;AUC_HA;AUC_AA};
Sens = {Sens_HH;Sens_HA;Sens_AA};
Spec = {Spec_HH;Spec_HA;Spec_AA};
AUCErr = {[];AUCerr_HA;[]};
SpecCI = {[];SpecCI_HA;[]};
FPR = {[];FPR_HA;[]};
FNR = {[];FNR_HA;[]};


end

function plotConfmat(confmat,figtitle,varargin)
if length(varargin) < 1
    figure;
else
    figure(gcf), subplot(2,2,varargin{1})
end

activities = {'Non-Fall','Fall'};
imagesc(confmat./repmat(sum(confmat,2),[1 2]));
confmat = confmat./repmat(sum(confmat,2),[1 2]);
% [cmin,cmax] = caxis;
caxis([0,1]) %set colormap to 0 1
ax = gca;
ax.XTick = 1:size(activities,2);
ax.YTick = 1:size(activities,2);
set(gca,'XTickLabel',activities,'FontSize',14)
set(gca,'YTickLabel',activities,'FontSize',14)
ax.XTickLabelRotation = 45;
axis square
if isempty(figtitle)
    title('Fall Detection')
else
    title(figtitle);
end
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
end

function string=getMeanSEMStr(X) % returns string with Mean +/- SEM for a vector
    string=sprintf('%0.3f +/- %0.3f',nanmean(X),1.96*nanstd(X)/sqrt(sum(~isnan(X))));
end

function cellstring=getMeanCIStr(Stat,AUC_CI,S_CI)
    cellstring=cell(3,1);
    cellstring{1}=sprintf('%0.3f [%0.3f - %0.3f]',nanmean(Stat(:,2)),AUC_CI(1),AUC_CI(2));
    cellstring{2}=sprintf('%0.3f [%0.3f - %0.3f]',nanmean(Stat(:,3)),S_CI(1,2),S_CI(2,2));
    cellstring{3}=sprintf('%0.3f [%0.3f - %0.3f]',nanmean(Stat(:,4)),S_CI(1,3),S_CI(2,3));
end
