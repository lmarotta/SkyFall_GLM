%% Threshold Method for classification

% Uses Thresholds on the acceleration magnitude to classify falls
% UFT - Upper Fall Threshold 
% LFT - Lower Fall Threshold
% data - cell array of acceleration clips

function isfall=ThresholdDetection(UFT,LFT,data)

    MaxAcc = cellfun(@(x) sqrt(max(sum(x(:,2:end).^2,2))),data);
    MinAcc = cellfun(@(x) sqrt(min(sum(x(:,2:end).^2,2))),data);

    isfall=MaxAcc>UFT & MinAcc<LFT;
end
    