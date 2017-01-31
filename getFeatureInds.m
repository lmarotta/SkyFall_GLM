
% INPUT: groups
%   18x1 logical vector (1 - use feature group, 0 - don't use)
%   feature groups:
%       1 - Raw Signal Statistics
%       2 - Raw Signal Correlation Coefficients
%       3 - Raw Signal 5s FFT bins
%       4 - Raw Signal 1s FFT bins
%		5 - Derivative Statistics
%		6 - Derivative 5s FFT bins
%		7 - Derivative 1s FFT bins
%		8 - Resultant Vector and Magnitude 
%		9 - Angle Statistics (ArcTan)
%		10 - Entropies
%		11 - Raw Signal Cross Products
%		12 - Derivative Cross Products
%		13 - Raw Signal Statistics on 1s FFT bins
%		14 - Raw Signal Entropies on 1s FFT bins
%		15 - Raw Signal Statistics on 1s binned signal energy
%		16 - Derivative Statistics on 1s FFT bins
%		17 - Derivative Entropies on 1s FFT bins
%		18 - Barometer

function featureInds=getFeatureInds(groups)

    grouplengths=[15 6 75 75 15 75 75 20 21 12 24 24 60 15 15 60 15 7];

    groups = logical(groups);
    featureInds=[];
    
    % Gyro
    for i=1:length(groups)-1
        featureInds=[featureInds; repmat(groups(i),grouplengths(i),1)];
    end
    % Acce
    for i=1:length(groups)-1
        featureInds=[featureInds; repmat(groups(i),grouplengths(i),1)];
    end
    % Baro
    featureInds=[featureInds; repmat(groups(end),grouplengths(end),1)];
    
end