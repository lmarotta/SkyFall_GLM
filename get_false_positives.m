load FallProbe_TestData;
load class_params_ACT_nobar;

plot = 1;

features = labels.features;
features = features(:,2:end); %delete first column with all zeros
% standardizing phone features
% if baro features are used by model, 
% comment the following line out
nz_ind_plus_baro = [nz_ind; zeros(78,1)];

mu = fvar.mu(logical(nz_ind_plus_baro));
std = fvar.std(logical(nz_ind_plus_baro));
features = (features - repmat(mu,size(features,1),1))./repmat(std,size(features,1),1);

% predicted values
conf= glmval(b, features, 'logit');
pred= ceil(conf-0.5);

false_positives_ind = logical(pred);
data.acce = labels.acce(false_positives_ind);
data.gyro = labels.gyro(false_positives_ind);
data.baro = labels.baro(false_positives_ind);
data.acce_interp = labels.acce_inerp(false_positives_ind);
data.gyro_interp = labels.gyro_inerp(false_positives_ind);
data.features = labels.features(false_positives_ind, true(1,size(labels.features,2)));
data.timestampSTART_END = labels.timestampSTART_END(false_positives_ind,logical([1 1]));

save FalsePositives data

if plot
    % plot data with predictions
    for i=1:length(labels.acce)
       labels.acce{i}(:,1) = labels.acce{i}(:,1) + labels.timestampSTART_END(i,1)/1000;
    end
    acce_data = vertcat(labels.acce{:});
    t = acce_data(:,1) - labels.timestampSTART_END(1,1)/1000;
    figure, plot(t,acce_data(:,2), t, acce_data(:,3), t,acce_data(:,4))
    hold on
    t1 = linspace(0,t(end),length(labels.acce));
    plot(t1,pred,'ro');

    % plot false positives
    for i=1:length(data.acce)
        clip = data.acce{i};
        t = clip(:,1);
       figure, plot(t,clip(:,2), t,clip(:,3), t,clip(:,4)); 
    end
end


