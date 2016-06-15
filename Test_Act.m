%Test trained GLMnet model with activities data
function Test_Act
close all
%% Load the classifier parameters  
% load class_params_dec.mat
%  List of variables in class_params:
%   b:  coefficients of logistic regression classifier
%   fvar: a structure containing the necessary vaiables for feature
%   nz_ind: indices of selected features extraction

%the list of jittered clips
files = dir('labels_plus_data_ACT.mat');

clear id conf
filename = files.name;
load(filename), disp(filename)

%sensors data, each k corresponds to 1 data point
skips=[];
for k = 1:length(labels.acce)
    accel_data = labels.acce{k};
    gyro_data = labels.gyro{k};
    baro_data = labels.baro{k};
    if size(gyro_data,1)>=100 && size( accel_data ,1)>=100 && size( baro_data,1 )>=10
        [id(k), conf(k)]= fall_model_eval(accel_data, gyro_data,  baro_data);
    else
        skips=[skips k];
        continue
    end      
end

id(skips)=[];
conf(skips)=[];
labels.value(skips)=[];

isfall = labels.value < 5;
fall_err = sum(~id(isfall))/sum(isfall);
nfall_err = sum(id(~isfall))/sum(~isfall);
err = (sum(id ~= isfall))/length(id); %error rate
confmat(:,:)=confusionmat(id==1,isfall);
    
confmat
% err = err([1 2 4 3 6 5 7]); %last corresponds to uniform distribution
% stdwin = [.01 .1 1 1.5  2 2.5];
% 
% %plot acc with unif jittered clips
% 
% figure, plot(stdwin,err(1:6),'o-')
% 
% figure, plot([stdwin(1) 5], err([1 7]),'o-')
% 
% figure, imagesc(confmat(:,:,1)./repmat(sum(confmat(:,:,1),2),[1 2])), colorbar
% 
% caxis([0 1])
% set(gca,'XTick', [1 2]), set(gca, 'XTickLabels', {'No Fall', 'Fall'})
% set(gca,'YTick', [1 2]), set(gca, 'YTickLabels', {'No Fall', 'Fall'})
% 
% figure, imagesc(confmat(:,:,7)./repmat(sum(confmat(:,:,7),2),[1 2])), colorbar
% 
% caxis([0 1])
% set(gca,'XTick', [1 2]), set(gca, 'XTickLabels', {'No Fall', 'Fall'})
% set(gca,'YTick', [1 2]), set(gca, 'YTickLabels', {'No Fall', 'Fall'})


