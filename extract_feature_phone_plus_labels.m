function [FN, fl, L, fold, mu, stdF, eps, nzstd]= extract_feature_phone_plus_labels(labels, fold_nr, indtrain)

eps=1e-6;
labelsNum= numel(labels.value);

F=[];
L=[];
sub_ind_vect=[];
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



% for i=1:labelsNum
% 
% 
% 
% %%  returns empty FN if there is not enough (<2) data points in  the data structure.
%  if size(labels.gyro{i},1)>=100 && size( labels.acce{i} ,1)>=100 && size( labels.baro{i},1 )>=10 %&& cis
%         
%         data{1}= labels.gyro{i}(:, 2:end);
%         data{2}= labels.acce{i}(:, 2:end);
%         data{3}= labels.baro{i}(:, 2:end);
%         stamp{1}= labels.gyro{i}(:, 1);
%         stamp{2}= labels.acce{i}(:, 1);
%         stamp{3}= labels.baro{i}(:, 1);
    %%
%     FC=[];
%     FC_label={};
%     sensors = {'gyro_', 'acce_','baro_'};
    %%
%     for j=1:3
%         if j<=2
%             %---------------------------------------------------------------------------------------------------------------------------------
%             %Extracting Features for accelerometer and Gyro
%             %---------------------------------------------------------------------------------------------------------------------------------
%             %% Sorting data in terms of the time stamps
%             [stamp{j}, ind]= sort(stamp{j}, 'ascend');
%             data{j}= data{j}(ind, :);
%             %% First column: Elapsed time from the first recording
%             DC=[stamp{j}-stamp{j}(1),data{j}];
%             %% mean
%             DCM= mean(DC,1);
%             DCM_label = repmat({strcat(sensors(1,j),'mean')},1,length(DCM));
%             %% median
%             DCMed= median(data{j},1);
%             DCMed_label = repmat({strcat(sensors(1,j),'median')},1,length(DCMed));
%             
%             %%  the real, imaginary and absolute value of the first 40 components   of fast fourier transform
%             DCFFT= fft(data{j},40);
%             DCFFT_re= real(DCFFT(:)');
%             DCFFT_im= imag(DCFFT(:)');
%             DCFFT_abs= abs(DCFFT(:)');
%             
%             DCFFT_re_label = repmat({strcat(sensors(1,j),'real_val_40comp_FFT')},1,length(DCFFT_re));
%             DCFFT_im_label = repmat({strcat(sensors(1,j),'imag_val_40comp_FFT')},1,length(DCFFT_im));
%             DCFFT_abs_label = repmat({strcat(sensors(1,j),'abs_val_40comp_FFT')},1,length(DCFFT_abs));
%             
%             %% Fitting the recodings as the function of time
% %             bnum=100;
% %             alp=1e-4;
% %             mu= linspace(0,5,bnum);
% %             var=(20/bnum)^2;
% %             dsz=size(DC,1);
% %             Tmat= repmat( DC(:,1), 1,bnum);
% %             Mmat= repmat( mu, dsz, 1);
% %             Phi= [exp(-(Tmat-Mmat).^2./var),ones(dsz,1)];
% %             Wmat= pinv(Phi'*Phi+alp*eye(bnum+1))*Phi'*DC(:,2:4);
% %             Wvect=Wmat(:)';            
%             
%             DCfit121=  polyfit(DC(:,1),DC(:,2), 1);
%             DCfit131=  polyfit(DC(:,1),DC(:,3), 1);
%             DCfit141=  polyfit(DC(:,1),DC(:,4), 1);            
%             DCfit122= polyfit(DC(:,1),DC(:,2), 2);
%             DCfit132= polyfit(DC(:,1),DC(:,3), 2);
%             DCfit142= polyfit(DC(:,1),DC(:,4), 2);
%             DCfit123= polyfit(DC(:,1),DC(:,2), 3);
%             DCfit133= polyfit(DC(:,1),DC(:,3), 3);
%             DCfit143= polyfit(DC(:,1),DC(:,4), 3);
%             DCfit124= polyfit(DC(:,1),DC(:,2), 4);
%             DCfit134= polyfit(DC(:,1),DC(:,3), 4);
%             DCfit144= polyfit(DC(:,1),DC(:,4), 4);
%             DCfit125= polyfit(DC(:,1),DC(:,2), 5);
%             DCfit135= polyfit(DC(:,1),DC(:,3), 5);
%             DCfit145= polyfit(DC(:,1),DC(:,4), 5);
%             
%             DCfit=[DCfit121,  DCfit131, DCfit141, DCfit122,  DCfit132, DCfit142, DCfit123,...
%                 DCfit133, DCfit143, DCfit124,  DCfit134, DCfit144, DCfit125,  DCfit135, DCfit145];
%             
%             DCfit121_label = repmat({strcat(sensors(1,j),'fit_func_121')},1,length(DCfit121));
%             DCfit131_label = repmat({strcat(sensors(1,j),'fit_func_131')},1,length(DCfit131));
%             DCfit141_label = repmat({strcat(sensors(1,j),'fit_func_141')},1,length(DCfit141));
%             DCfit122_label = repmat({strcat(sensors(1,j),'fit_func_122')},1,length(DCfit122));
%             DCfit132_label = repmat({strcat(sensors(1,j),'fit_func_132')},1,length(DCfit132));
%             DCfit142_label = repmat({strcat(sensors(1,j),'fit_func_142')},1,length(DCfit142));
%             DCfit123_label = repmat({strcat(sensors(1,j),'fit_func_123')},1,length(DCfit123));
%             DCfit133_label = repmat({strcat(sensors(1,j),'fit_func_133')},1,length(DCfit133));
%             DCfit143_label = repmat({strcat(sensors(1,j),'fit_func_143')},1,length(DCfit143));
%             DCfit124_label = repmat({strcat(sensors(1,j),'fit_func_124')},1,length(DCfit124));
%             DCfit134_label = repmat({strcat(sensors(1,j),'fit_func_134')},1,length(DCfit134));
%             DCfit144_label = repmat({strcat(sensors(1,j),'fit_func_144')},1,length(DCfit144));
%             DCfit125_label = repmat({strcat(sensors(1,j),'fit_func_125')},1,length(DCfit125));
%             DCfit135_label = repmat({strcat(sensors(1,j),'fit_func_135')},1,length(DCfit135));
%             DCfit145_label = repmat({strcat(sensors(1,j),'fit_func_145')},1,length(DCfit145));
%             
%             DCfit_label = {DCfit121_label{:}, DCfit131_label{:}, DCfit141_label{:}, DCfit122_label{:}, DCfit132_label{:}, DCfit142_label{:}, ...
%                 DCfit123_label{:}, DCfit133_label{:}, DCfit143_label{:}, DCfit124_label{:}, DCfit134_label{:}, DCfit144_label{:}, ...
%                 DCfit125_label{:}, DCfit135_label{:}, DCfit145_label{:}};
%             
%             %% computing the correlation coefficients between different channels of the sensors
%             DCCorr= corrcoef(DC);
%             DCCorr(isnan( DCCorr))=0;
%             DCCorr=DCCorr+eps;
%             DCCorrup= triu(DCCorr,1);
%             DCCorrInd= find(abs(DCCorrup)>0);
%             DCCorrV= DCCorrup(DCCorrInd)';
%             DCCorrV_label = repmat({strcat(sensors(1,j),'corr_coef')},1,length(DCCorrV));
%             
%             %% standard deviation of channels
%             DCstd= std(DC,0,1);
%             DCstd_label = repmat({strcat(sensors(1,j),'std')},1,length(DCstd));
%             %% skewness od data (3rd order moment)
%             DCS= skewness(DC);
%             DCS_label = repmat({strcat(sensors(1,j),'skewness')},1,length(DCS));
%             %% Kurtosis of the dataset
%             DCK= kurtosis(DC);
%             DCK_label = repmat({strcat(sensors(1,j),'kurtosis')},1,length(DCK));
%             %%
%             %--------------------------------------------------------------------------------------------------------------------------------------
%             %Computing the statistics of the derivative of data
%             %--------------------------------------------------------------------------------------------------------------------------------------
%             %%computing the  derivates (DSLN)
%             DDC= diff(data{j},1);
%             stamp_diff= diff(stamp{j});
%             IDSL= (stamp_diff>0);
%             DSLN= DDC(IDSL,:)./repmat(stamp_diff(IDSL),1,3);
%             %% Computing the mean of derivatives
%             DDCM= mean(DSLN,1);
%             DDCM_label = repmat({strcat(sensors(1,j),'mean_of_deriv')},1,length(DDCM));
%             %% Computing the median of derivatives
%             DDCMed= median(DSLN,1);
%             DDCMed_label = repmat({strcat(sensors(1,j),'med_of_deriv')},1,length(DDCMed));
%             %%  The real, imaginary and absolute value of the first 40 components   of fast fourier transform  of derivatives
%             DDCFFT= fft(DSLN,40);
%             DDCFFT_re= real(DDCFFT(:)');
%             DDCFFT_im= imag(DDCFFT(:)');
%             DDCFFT_abs= abs(DDCFFT(:)');
%             DDCFFT_re_label = repmat({strcat(sensors(1,j),'real_val_40comp_FFT_deriv')},1,length(DDCFFT_re));
%             DDCFFT_im_label = repmat({strcat(sensors(1,j),'imag_val_40comp_FFT_deriv')},1,length(DDCFFT_im));
%             DDCFFT_abs_label = repmat({strcat(sensors(1,j),'abs_val_40comp_FFT_deriv')},1,length(DDCFFT_abs));
%             
%             %% Fitting the derivatives of recording signals as a function of time
%             TD= stamp{j}(IDSL)-stamp{j}(ind(1));
% %             bnum=100;
% %             alp=1e-5;
% %             muD= linspace(0,5,bnum);
% %             varD=(20/bnum)^2;
% %             dDsz=size(DSLN,1);
% %             TDmat= repmat( TD, 1,bnum);
% %             MDmat= repmat( muD, dDsz, 1);
% %             Phi= [exp(-(TDmat-MDmat).^2./varD),ones(dDsz,1)];
% %             WDmat= pinv(Phi'*Phi+alp*eye(bnum+1))*Phi'*DSLN;
% %             WDvect=WDmat(:)';
%             
%             DDCfit121= polyfit(TD,DSLN(:,1), 1);
%             DDCfit131= polyfit(TD,DSLN(:,2), 1);
%             DDCfit141= polyfit(TD,DSLN(:,3), 1);
%             DDCfit122= polyfit(TD,DSLN(:,1), 2);
%             DDCfit132= polyfit(TD,DSLN(:,2), 2);
%             DDCfit142= polyfit(TD,DSLN(:,3), 2);
%             DDCfit123= polyfit(TD,DSLN(:,1), 3);
%             DDCfit133= polyfit(TD,DSLN(:,2), 3);
%             DDCfit143= polyfit(TD,DSLN(:,3), 3);
%             DDCfit= [DDCfit121, DDCfit131,  DDCfit141, DDCfit122, DDCfit132,  DDCfit142,...
%                 DDCfit123, DDCfit133,  DDCfit143];
%             
%             DDCfit121_label = repmat({strcat(sensors(1,j),'fit_func_121_deriv')},1,length(DDCfit121));
%             DDCfit131_label = repmat({strcat(sensors(1,j),'fit_func_131_deriv')},1,length(DDCfit131));
%             DDCfit141_label = repmat({strcat(sensors(1,j),'fit_func_141_deriv')},1,length(DDCfit141));
%             DDCfit122_label = repmat({strcat(sensors(1,j),'fit_func_122_deriv')},1,length(DDCfit122));
%             DDCfit132_label = repmat({strcat(sensors(1,j),'fit_func_132_deriv')},1,length(DDCfit132));
%             DDCfit142_label = repmat({strcat(sensors(1,j),'fit_func_142_deriv')},1,length(DDCfit142));
%             DDCfit123_label = repmat({strcat(sensors(1,j),'fit_func_123_deriv')},1,length(DDCfit123));
%             DDCfit133_label = repmat({strcat(sensors(1,j),'fit_func_133_deriv')},1,length(DDCfit133));
%             DDCfit143_label = repmat({strcat(sensors(1,j),'fit_func_143_deriv')},1,length(DDCfit143));
%             DDCfit_label = {DDCfit121_label{:}, DDCfit131_label{:}, DDCfit141_label{:}, DDCfit122_label{:}, DDCfit132_label{:}, DDCfit142_label{:}, ...
%                 DDCfit123_label{:}, DDCfit133_label{:}, DDCfit143_label{:}};
%             
%             %% correlation coefficients between different channels of the sensors
%             DDCCorr= corrcoef(DSLN)+eps;
%             DDCCorr(isnan( DDCCorr))=0;
%             DDCCorr= DDCCorr+eps;
%             DDCCorrup= triu(DDCCorr,1);
%             DDCCorrInd= find(abs(DDCCorrup)>0);
%             DDCCorrV= DDCCorrup(DDCCorrInd)';
%             DDCCorrV_label = repmat({strcat(sensors(1,j),'corr_coef_deriv')},1,length(DDCCorrV));
%             
%             %%  standard deviations of derivatives
%             DDCstd= std(DSLN,0,1);
%             DDCstd_label = repmat({strcat(sensors(1,j),'std_deriv')},1,length(DDCstd));
%             %% skewness of derivatives
%             DDCS= skewness(DSLN);
%             DDCS_label = repmat({strcat(sensors(1,j),'skewness_deriv')},1,length(DDCS));
%             %% kurtosis of derivatives
%             DDCK= kurtosis(DSLN);
%             DDCK_label = repmat({strcat(sensors(1,j),'kurtosis_deriv')},1,length(DDCK));
%             
%             %% stacking the features of the  sensors (Gyro or accelerometer) together
%             FC= [FC, DCM,  DCMed, DCFFT_re, DCFFT_im,  DCFFT_abs, DCfit,  DCCorrV, DCstd DCS, DCK, DDCM, DDCMed, DDCCorrV, DDCstd, DDCS, DDCK, DDCFFT_re, DDCFFT_im,  DDCFFT_abs, DDCfit];
%             FC_label = {FC_label{:}, DCM_label{:},  DCMed_label{:}, DCFFT_re_label{:}, DCFFT_im_label{:},  DCFFT_abs_label{:}, DCfit_label{:},  DCCorrV_label{:}, ...
%                 DCstd_label{:}, DCS_label{:}, DCK_label{:}, DDCM_label{:}, DDCMed_label{:}, DDCCorrV_label{:}, DDCstd_label{:}, DDCS_label{:}, DDCK_label{:}, DDCFFT_re_label{:},...
%                 DDCFFT_im_label{:},  DDCFFT_abs_label{:}, DDCfit_label{:}};
%             %%
%         else
%             %--------------------------------------------------------------------------------------------------------------------------------------
%             %  Extracting  features for barometer
%             %--------------------------------------------------------------------------------------------------------------------------------------
%             
% %             continue % skip barometer features
%             %% sorting the  data in terms of their time stamps
%             [stamp{j}, ind]= sort(stamp{j},'ascend');
%             data{j}= data{j}(ind, :);
%             %% First column: Elapsed time from the first recording
%             DC=[stamp{j}-stamp{j}(1),data{j}];
%             %%  Mean
%             DCM= mean(DC,1);
%             DCM_label = repmat({strcat(sensors(1,j),'mean')},1,length(DCM));
%             %% Median
%             DCMed=median(data{j},1);
%             DCMed_label = repmat({strcat(sensors(1,j),'median')},1,length(DCMed));
%             %% Correlation coefficients between the sensors and time
%             DCCorr= corrcoef(DC);
%             DCCorr(isnan( DCCorr))=0;
%             DCCorr= DCCorr+eps;
%             DCCorrup= triu(DCCorr,1);
%             DCCorrInd= find(abs(DCCorrup)>0);
%             DCCorrV= DCCorrup(DCCorrInd)';
%             DCCorrV_label = repmat({strcat(sensors(1,j),'corr_coef')},1,length(DCCorrV));
%             %% The real, imaginary and absolute value of the first 5 components   of fast fourier transform
%             DCFFT= fft(data{j},5);
%             DCFFT_re= real(DCFFT(:)');
%             DCFFT_im= imag(DCFFT(:)');
%             DCFFT_abs= abs(DCFFT(:)');
%             
%             DCFFT_re_label = repmat({strcat(sensors(1,j),'real_val_5comp_FFT')},1,length(DCFFT_re));
%             DCFFT_im_label = repmat({strcat(sensors(1,j),'imag_val_5comp_FFT')},1,length(DCFFT_im));
%             DCFFT_abs_label = repmat({strcat(sensors(1,j),'abs_val_5comp_FFT')},1,length(DCFFT_abs));
%             
%             %% Fitting the  recording as the function of time
%             DCfit121= polyfit(DC(:,1),DC(:,2), 1);
%             DCfit131= polyfit(DC(:,1),DC(:,3), 1);
%             DCfit122= polyfit(DC(:,1),DC(:,2), 2);
%             DCfit132= polyfit(DC(:,1),DC(:,3), 2);
%             DCfit123= polyfit(DC(:,1),DC(:,2), 3);
%             DCfit133= polyfit(DC(:,1),DC(:,3), 3);
%             DCfit124= polyfit(DC(:,1),DC(:,2), 4);
%             DCfit134= polyfit(DC(:,1),DC(:,3), 4);
%                    
%             DCfit=[DCfit121,  DCfit131, DCfit122,  DCfit132,  DCfit123,...
%                 DCfit133,  DCfit124,  DCfit134];
%             
%             DCfit121_label = repmat({strcat(sensors(1,j),'fit_func_121')},1,length(DCfit121));
%             DCfit131_label = repmat({strcat(sensors(1,j),'fit_func_131')},1,length(DCfit131));
%             DCfit122_label = repmat({strcat(sensors(1,j),'fit_func_122')},1,length(DCfit122));
%             DCfit132_label = repmat({strcat(sensors(1,j),'fit_func_132')},1,length(DCfit132));
%             DCfit123_label = repmat({strcat(sensors(1,j),'fit_func_123')},1,length(DCfit123));
%             DCfit133_label = repmat({strcat(sensors(1,j),'fit_func_133')},1,length(DCfit133));
%             DCfit124_label = repmat({strcat(sensors(1,j),'fit_func_124')},1,length(DCfit124));
%             DCfit134_label = repmat({strcat(sensors(1,j),'fit_func_134')},1,length(DCfit134));
%             
%             DCfit_label = {DCfit121_label{:}, DCfit131_label{:}, DCfit122_label{:}, DCfit132_label{:}, ...
%                 DCfit123_label{:}, DCfit133_label{:}, DCfit124_label{:}, DCfit134_label{:}};
%             
%             %% standard deviation
%             DCstd = std(DC,0,1);
%             DCstd_label = repmat({strcat(sensors(1,j),'std')},1,length(DCstd));
%             %% 
%             DCskew= skewness(DC);
%             DCskew_label = repmat({strcat(sensors(1,j),'skewness')},1,length(DCskew));
%             %% Kurtosis of data
%             DCkurt= kurtosis(DC);
%             DCkurt_label = repmat({strcat(sensors(1,j),'kurtosis')},1,length(DCkurt));
%             
%             %% derivative of data
%             DDC= diff(data{j},1);
%             
%             stamp_diff= diff(stamp{j});
%             
%             
%             IDSL= (stamp_diff>0);
%             %% checks whether there is enough data to compute the derivative
%             if any(IDSL>0)
%                 DSLN= DDC(IDSL,:)./ repmat(stamp_diff(IDSL),1,2);
%                 %% mean
%                 DDCM= mean(DSLN,1);                
%                 %% median
%                 DDCMed= median(DSLN,1);                
%                 %% standard deviation
%                 DDCstd= std(DSLN,0,1);                
%             else
%                 %% mean
%                 DDCM= zeros(1,2);
%                 %% median
%                 DDCstd= zeros(1,2);
%                 %% standard deviation
%                 DDCMed=zeros(1,2);
%             end
%             DDCM_label = repmat({strcat(sensors(1,j),'mean_of_deriv')},1,length(DDCM));
%             DDCMed_label = repmat({strcat(sensors(1,j),'med_of_deriv')},1,length(DDCMed));
%             DDCstd_label = repmat({strcat(sensors(1,j),'std_deriv')},1,length(DDCstd));
%             
%             %% Stacking  all the features of barometer together
%             FC= [FC, DCM, DCMed DCCorrV, DCfit, DCstd, DCskew, DCkurt, DDCM, DDCMed, DDCstd, DCFFT_re, DCFFT_im, DCFFT_abs];
%             FC_label= {FC_label{:}, DCM_label{:}, DCMed_label{:}, DCCorrV_label{:}, DCfit_label{:}, DCstd_label{:}, DCskew_label{:}, DCkurt_label{:},...
%                 DDCM_label{:}, DDCMed_label{:}, DDCstd_label{:}, DCFFT_re_label{:}, DCFFT_im_label{:}, DCFFT_abs_label{:}};
%         end;
%         
%     end;

F=csvread('full_feature_set.csv');
F=F(indtrain,1:1781);

L=labels.value;

% find zeroed rows in F
z_inds=max(abs(F).')==0;

F(z_inds,:)=[];
L(z_inds)=[];
subjects=labels.subject;
subjects(z_inds)=[];

for i=1:length(subjects)

    Amp_flag =regexp(subjects{i},'AF','ignorecase');

    if isempty(Amp_flag)
        sub_ind=str2num(subjects{i}(3:4))+7;
    else
        sub_ind=str2num(subjects{i}(3:4));
    end;

    sub_ind_vect = [sub_ind_vect; sub_ind];

end




% F=F(1:end-50,:);
% sub_ind_vect=sub_ind_vect(1:end-50);
% L=L(1:end-50);
dsz= size(F,1);
v1= ones(dsz,1);
mu= mean(F);
stdF= std(F);


%%%%Standardize the result

nzstd=stdF>eps*mean(stdF);
% remove barometer features
nzstd(1701:1781)=zeros(1,81);
mu=mu(nzstd);
F=F(:,nzstd);
fl = [];
stdF=stdF(nzstd);

FN=(F-v1*mu)./(v1*stdF);

fold=zeros(dsz,1);


[~,id]=sort(sub_ind_vect,'ascend');
FN= FN(id, :);
L= L(id);

%fold_ind= linspace(1,dsz+1,folds_nr);
%fold_sz= floor(dsz/fold_nr);

fold_grid=round(linspace(1,dsz,fold_nr+1));



for i=1:fold_nr
 fold(fold_grid(i):fold_grid(i+1))=i;
end;


return;