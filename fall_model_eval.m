function [id, conf]= fall_model_eval(accel_data, gyro_data,  baro_data)

%------------------------------------------------------------------------------------
% Fall prediction model: predicts fall using a generalized linear model
% trained by elastic net logisitic regression
%-----------------------------------------------------------------------------------
%
% INPUT ARGUMENTS:
%  gyro_data : Tg * 4 matrix, Tg number of time stamps, 1st column: time stamps, 2nd, 3rd and 4th are the the sensor readings (normalized time stamps)  
%  accel_data: Ta * 4 matrix, Ta number of time stamps, 1st column: time stamps, 2nd, 3rd and 4th are the the sensor readings (normalized time stamps)
%  baro_data : Tb * 3 matrix, Tb number of time stamps, 1st column: time stamps, 2nd and 3rd are  the sensor readings (normalized time stamps)
%
% OUTPUT ARGUMENTS:                 
%  id: 0 for no fall,  1 for fall
%  conf: 0 no fall,   0.5  toss up,   1 fall  



%% Load the classifier parameters  
%  List of variables in class_params:
%   b:  coefficients of logistic regression classifier
%   fvar: a structure containing the necessary vaiables for feature
%   nz_ind: indices of selected features  
% extraction

load class_params
%load class_params_dec  %original model parameters 

%% Extract features 
% FN: Sparse features extracted for the data

% FN= extract_feature_test(accel_data, gyro_data, baro_data, fvar);
FN= extract_feature_test_dec(accel_data, gyro_data, baro_data, fvar);


% remove the redundant features
FN_nz= FN(:, nz_ind);


%% Run the model
%-----------------------------------------------------------------------------------------------------------------------------------
% runs a standard matlab classifier model( multinomial  logistic regression model)  with the set of parameters w obtained  by glmnet  
%----------------------------------------------------------------------------------------------------------------------------------

conf= glmval(b, FN_nz, 'logit');
id= ceil(conf-0.5);

return;