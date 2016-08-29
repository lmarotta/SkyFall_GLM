%% RotateData.m
% Rotates the accelerometer and gyro data randomly along the z-axis

load labels_plus_data_ACT

X=rand(size(labels.acce,2))*180-90; % random angles between +/-90

for i=1:size(labels.acce,2)
    if ~isempty(labels.acce{i}) && ~isempty(labels.gyro{i})
        rot=[cosd(X(i)) sind(X(i)) 0; -sind(X(i)) cosd(X(i)) 0; 0 0 1];
        labels.acce{i}=[labels.acce{i}(:,1) labels.acce{i}(:,2:end)*rot];
        labels.gyro{i}=[labels.gyro{i}(:,1) labels.gyro{i}(:,2:end)*rot];
    end
end

save labels_plus_data_ACT_rot labels