%% Threshold Method for classification

% Uses Thresholds on the acceleration magnitude to classify falls
% UFT - Upper Fall Threshold 
% LFT - Lower Fall Threshold
% data - cell array of acceleration clips

function [FPR,FNR]=ThresholdDetection(TrainData,TestData)

    UFT = min(cellfun(@(x) sqrt(max(sum(x(:,2:end).^2,2))),TrainData.acce(TrainData.value<9)));
    LFT = max(cellfun(@(x) sqrt(min(sum(x(:,2:end).^2,2))),TrainData.acce(TrainData.value<9)));

    MaxAcc = cellfun(@(x) sqrt(max(sum(x(:,2:end).^2,2))),TestData.acce);
    MinAcc = cellfun(@(x) sqrt(min(sum(x(:,2:end).^2,2))),TestData.acce);

    isfall=MaxAcc>UFT & MinAcc<LFT;
    
    FNR = sum(~isfall(TestData.value<9))/sum(TestData.value<9);
    FPR = sum(isfall(TestData.value==9))/sum(TestData.value==9);
    
end
    