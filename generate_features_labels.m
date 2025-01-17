featuresLabels = cell(1, 1781);

sensors = {'gyro_', 'acce_','baro_'};

%features for gyro and accelerometer
fullFeturesSetForOneSensor = {'mean', 'median', 'real_val_40comp_FFT', 'imag_val_40comp_FFT', 'abs_val_40comp_FFT', 'fit_func_121', 'fit_func_131', ...
    'fit_func_141', 'fit_func_122', 'fit_func_132', 'fit_func_142', 'fit_func_123', 'fit_func_133', 'fit_func_143', 'fit_func_124', ...
    'fit_func_134', 'fit_func_144', 'fit_func_125', 'fit_func_135', 'fit_func_145', 'corr_coef', 'std', 'skewness', 'kurtosis', ...
    'mean_of_deriv', 'med_of_deriv', 'corr_coef_of_deriv', 'std_of_deriv', 'skewness_of_deriv', 'kurtosis_of_deriv', ...
    'real_val_40comp_FFT_deriv', 'imag_val_40comp_FFT_deriv', 'abs_val_40comp_FFT_deriv', 'fit_func_121_deriv', 'fit_func_131_deriv', 'fit_func_141_deriv', ...
    'fit_func_122_deriv', 'fit_func_132_deriv', 'fit_func_142_deriv', 'fit_func_123_deriv', 'fit_func_133_deriv', 'fit_func_143_deriv'};

numOfRepFeatures = [4 3 120 120 120 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 6 4 4 4 3 3 3 3 3 3 120 120 120 2 2 2 3 3 3 4 4 4];
% 1 - to add timestamp or axis to a label; 2- to add axis to a label; 0 - no info about axis
startAxisIndex =   [1 2 2   2   2   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 2 2 2 2 2 2 2   2   2   0 0 0 0 0 0 0 0 0]; 
listOfAxis = {'_timestamp', '_x', '_y', '_z'};

%features for barometer
featuresSetForBaro = {'mean', 'median', 'corr_coef', 'fit_func_121', 'fit_func_131', 'fit_func_122', 'fit_func_132', 'fit_func_123', 'fit_func_133', ...
    'fit_func_124', 'fit_func_134', 'std', 'skewness', 'kurtosis', 'mean_of_deriv', 'med_of_deriv', 'std_of_deriv', 'real_val_5comp_FFT', 'imag_val_5comp_FFT', ...
    'abs_val_5comp_FFT'};

numOfRepFeaturesBaro = [3 2 3 2 2 3 3 4 4 5 5 3 3 3 2 2 2 10 10 10];
% 1 - to add timestamp or axis to a label; 2- to add axis to a label; 0 - no info about axis
startAxisIndexBaro   = [1 2 1 0 0 0 0 0 0 0 0 1 1 1 2 2 2 2  2  2];
listOfAxisBaro = {'_timestamp', '_pressure', '_altitude'};

%generates vector with features labels
currIndex = 1;
for j=1:length(sensors)
    %gyro and accelerometer
    if j<=2
        featuresSet = fullFeturesSetForOneSensor;
        numOfRep = numOfRepFeatures;
        axis = listOfAxis;
        startIndecies = startAxisIndex;
    %barometer    
    else
        featuresSet = featuresSetForBaro;
        numOfRep = numOfRepFeaturesBaro;
        axis = listOfAxisBaro;
        startIndecies = startAxisIndexBaro;
    end
    for i=1:1:length(numOfRep)
        if startIndecies(1,i)==0
            for k=1:numOfRep(1,i)
                featuresLabels(currIndex) = strcat(sensors(1,j),featuresSet(1,i));
                currIndex=currIndex+1;
            end
        else
            axisStart = startIndecies(1,i);
            rep = numOfRep(1,i);
            repPerAxis = rep/(length(axis)-(axisStart-1));
            currAxis = axisStart;
            while currAxis<=length(axis)
                for l=1:repPerAxis
                    featuresLabels(currIndex) = strcat(sensors(1,j),featuresSet(1,i),axis(1,currAxis));
                    currIndex=currIndex+1;
                end
                currAxis=currAxis+1;
            end
        end
    end
end

%load class_params_ACT_nobar

%featuresLabels = featuresLabels(:,fvar.nzstd);

save features_labels_full featuresLabels