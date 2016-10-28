load PhoneProbe_data

t=0:.02:4.98;

% clip with large amount of gyroscope samples
gyro_large_clip = labels.gyro{20};
% remove duplicates
gyro_large_clip=sortrows(gyro_large_clip);
gyro_large_clip(diff(gyro_large_clip(:,1))==0,:)=[];
gyro_origin_interp=[t' spline(gyro_large_clip(:,1)'-gyro_large_clip(1,1),gyro_large_clip(:,2:end)',t')'];
if any(diff(gyro_large_clip(:,1))>0.25)
    disp('gap larger than 0.25 sec detected in the original clip')
end
% t1 = 1:size(gyro_large_clip,1);
%t1 = gyro_large_clip(:,1);

%figure, plot(t1,gyro_large_clip(:,2), t1,gyro_large_clip(:,3), t1,gyro_large_clip(:,4)), legend('X','Y','Z')
%xlabel('Time'), ylabel('Sensor Value')
%title('Gyro Clip with 489 Samples')

gyro_large_clip_resampled = datasample(gyro_large_clip, 200, 'Replace',false);
gyro_large_clip_resampled_sorted = sortrows(gyro_large_clip_resampled,1);
gyro_resampl_interp=[t' spline(gyro_large_clip_resampled_sorted(:,1)'-gyro_large_clip_resampled_sorted(1,1),gyro_large_clip_resampled_sorted(:,2:end)',t')'];

if any(diff(gyro_large_clip_resampled_sorted(:,1))>0.25)
    disp('gap larger than 0.25 sec detected in resampled clip');
end 

%t2 = gyro_large_clip_resampled_sorted(:,1);

%figure, plot(t2,gyro_large_clip_resampled_sorted(:,2), t2,gyro_large_clip_resampled_sorted(:,3), t2,gyro_large_clip_resampled_sorted(:,4)), legend('X','Y','Z')
%xlabel('Time'), ylabel('Sensor Value')
%title('Gyro Clip 489 Resampled to 200')

figure
subplot(3,1,1), plot(t,gyro_origin_interp(:,2), t,gyro_resampl_interp(:,2)), legend('Original', 'Resampled')
xlabel('Time'), ylabel('Sensor Value'), title('X')

subplot(3,1,2), plot(t,gyro_origin_interp(:,3), t,gyro_resampl_interp(:,3)), legend('Original', 'Resampled')
xlabel('Time'), ylabel('Sensor Value'), title('Y')

subplot(3,1,3), plot(t,gyro_origin_interp(:,4), t,gyro_resampl_interp(:,4)), legend('Original', 'Resampled')
xlabel('Time'), ylabel('Sensor Value'), title('Z')

% features test
accel = labels.acce{20};
accel=sortrows(accel);
accel(diff(accel(:,1))==0,:)=[];
accel_interp=[t' spline(accel(:,1)'-accel(1,1),accel(:,2:end)',t')'];

gyro_origin = gyro_origin_interp;
gyro_resampl = gyro_resampl_interp;

baro = labels.baro{20};
baro=sortrows(baro);
baro(diff(baro(:,1))==0,:)=[];
baro_interp=[t' spline(baro(:,1)'-baro(1,1),baro(:,2:end)',t')'];

load class_params_ACT_nobar

FN_origin = extract_feature_test_phone(accel_interp, gyro_origin, baro_interp, fvar, 1);
FN_resampl = extract_feature_test_phone(accel_interp, gyro_resampl, baro_interp, fvar, 1);
FN_origin_norm = (FN_origin-fvar.mu)./fvar.std;
FN_resampl_norm = (FN_resampl-fvar.mu)./fvar.std;

fn_diff = FN_origin_norm - FN_resampl_norm;

ratios = FN_origin.*(FN_resampl.^-1);
figure, bar(ratios)%, ylim([-1 3])
title('Ratios Features from Original and Resampled Gyro Data Clip')
xlabel('Features')
ylabel('Feature Value')

errors = (FN_origin-FN_resampl)./FN_origin;
figure, bar(errors)%, ylim([-1 3])
title('Errors Features from Original and Resampled Gyro Data Clip')
xlabel('Features')
ylabel('Error Value')

ratios_norm = FN_origin_norm.*(FN_resampl_norm.^-1);
figure, bar(ratios_norm)%, ylim([-1 3])
title('Ratios Norm Features from Original and Resampled Gyro Data Clip')
xlabel('Features')
ylabel('Feature Value')

errors_norm = (FN_origin_norm-FN_resampl_norm)./FN_origin_norm;
figure, bar(errors_norm)%, ylim([-1 3])
title('Errors Norm Features from Original and Resampled Gyro Data Clip')
xlabel('Features')
ylabel('Error Value')

load features_labels_full

fn_labels = featuresLabels(fvar.nzstd);
bad_fn = ratios >5 | ratios <-5;
bad_labels = fn_labels(bad_fn)';

bad_norm_fn = fn_diff > 1 | fn_diff < -1;
bad_norm_labels = fn_labels(bad_norm_fn)';

% resample 100 times

gyro_resampl_set = {};
FN_resampl_set = [];

for i=1:100 
   gyro_large_clip_resampled = datasample(gyro_large_clip, 200, 'Replace',false);
    gyro_large_clip_resampled_sorted = sortrows(gyro_large_clip_resampled,1);
    gyro_resampl_interp=[t' spline(gyro_large_clip_resampled_sorted(:,1)'-gyro_large_clip_resampled_sorted(1,1),gyro_large_clip_resampled_sorted(:,2:end)',t')'];

    if any(diff(gyro_large_clip_resampled_sorted(:,1))>0.25)
        disp('gap larger than 0.25 sec detected in resampled clip');
    end 
    gyro_resampl_set{i} = gyro_resampl_interp;
    FN_resampl = extract_feature_test_phone(accel_interp, gyro_resampl_interp, baro_interp, fvar, 1);
    FN_resampl_norm = (FN_resampl-fvar.mu)./fvar.std;
    FN_resampl_set = [FN_resampl_set; FN_resampl_norm];
end

fn_diff_100 = FN_origin_norm - mean(FN_resampl_set);
figure, plot(fn_diff_100);

bad_aver_resampl_fn = fn_diff_100 > 1 | fn_diff_100 < -1;
bad_aver_resampl_labels = fn_labels(bad_aver_resampl_fn)';




