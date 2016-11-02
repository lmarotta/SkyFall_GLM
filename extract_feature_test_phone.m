function FN= extract_feature_test_phone(accel_data, gyro_data, baro_data, fvar, new_FFT)
%----------------------------------------------------------------------------------------------------
% extracts features from  5 seconds interval of  sensor recordings (accelerometer, gyro and barometer)
%----------------------------------------------------------------------------------------------------
% INPUTS
%         gyro_data : Tg * 4 matrix, Tg number of time stamps, 1st column: time stamps, 2nd, 3rd and 4th are the the sensor readings (normalized time stamps)
%         accel_data: Ta * 4 matrix, Ta number of time stamps, 1st column: time stamps, 2nd, 3rd and 4th are the the sensor readings (normalized time stamps)
%         baro_data : Tb * 3 matrix, Tb number of time stamps, 1st column: time stamps, 2nd and 3rd are  the sensor readings (normalized time stamps)
%
%         fvar: contains the required vaiables for buliding  features  with
%         the following fields
%                fvar.eps: the resolution threshold
%                fvar.mean: mean of the training set data
%                fvar.std:  standard deviation of training set data
%                fvar.nzstd: columns with non-zero standard deviation in training set
%
% OUTPUT
%         FN: 1*p vector containing  Features, where p is the size of feature set is returned empty


%%  returns empty FN if there is not enough (<2) data points in  the data structure.
if size(gyro_data,1)>=100 && size( accel_data ,1)>=100 && size( baro_data,1 )>=10
    
    data{1}= gyro_data(:, 2:end);
    data{2}= accel_data(:, 2:end);
    data{3}= baro_data(:, 2:end);
    %%
    stamp{1}= gyro_data(:, 1);
    stamp{2}= accel_data(:, 1);
    stamp{3}= baro_data(:, 1);
    %%
    FC=CalcFeatures(data,stamp,new_FFT, fvar.eps);
    
    F= FC(fvar.nzstd);
    %% standardizing the data by subtracting the mean and dividing it by variance
    %% ( mean and variance are computed on the training set)
    FN= (F- fvar.mu)./ fvar.std;
else
    FN=[];
    
end;

return
