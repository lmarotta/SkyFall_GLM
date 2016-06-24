function FN= extract_feature_test_dec(accel_data, gyro_data, baro_data, fvar)
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
    FC=[];
    %%
    for j=1:3
        if j<=2
            %---------------------------------------------------------------------------------------------------------------------------------
            %Extracting Features for accelerometer and Gyro
            %---------------------------------------------------------------------------------------------------------------------------------
            %% Sorting data in terms of the time stamps
            [stamp{j}, ind]= sort(stamp{j}, 'ascend');
            data{j}= data{j}(ind, :);
            %% First column: Elapsed time from the first recording
            DC=[stamp{j}-stamp{j}(1),data{j}];
            %% mean
            DCM= mean(DC,1);
            %% median
            DCMed= median(data{j},1);
            %%  the real, imaginary and absolute value of the first 40 components   of fast fourier transform
            DCFFT= fft(data{j},40);
            DCFFT_re= real(DCFFT(:)');
            DCFFT_im= imag(DCFFT(:)');
            DCFFT_abs= abs(DCFFT(:)');
            
            %% Fitting the recodings as the function of time
%             bnum=100;
%             alp=1e-4;
%             mu= linspace(0,5,bnum);
%             var=(20/bnum)^2;
%             dsz=size(DC,1);
%             Tmat= repmat( DC(:,1), 1,bnum);
%             Mmat= repmat( mu, dsz, 1);
%             Phi= [exp(-(Tmat-Mmat).^2./var),ones(dsz,1)];
%             Wmat= pinv(Phi'*Phi+alp*eye(bnum+1))*Phi'*DC(:,2:4);
%             Wvect=Wmat(:)';
            
            
            
            
            DCfit121=  polyfit(DC(:,1),DC(:,2), 1);
            DCfit131=  polyfit(DC(:,1),DC(:,3), 1);
            DCfit141=  polyfit(DC(:,1),DC(:,4), 1);
            
            DCfit122= polyfit(DC(:,1),DC(:,2), 2);
            DCfit132= polyfit(DC(:,1),DC(:,3), 2);
            DCfit142= polyfit(DC(:,1),DC(:,4), 2);
            DCfit123= polyfit(DC(:,1),DC(:,2), 3);
            DCfit133= polyfit(DC(:,1),DC(:,3), 3);
            DCfit143= polyfit(DC(:,1),DC(:,4), 3);
            DCfit124= polyfit(DC(:,1),DC(:,2), 4);
            DCfit134= polyfit(DC(:,1),DC(:,3), 4);
            DCfit144= polyfit(DC(:,1),DC(:,4), 4);
            DCfit125= polyfit(DC(:,1),DC(:,2), 5);
            DCfit135= polyfit(DC(:,1),DC(:,3), 5);
            DCfit145= polyfit(DC(:,1),DC(:,4), 5);
            
            DCfit=[DCfit121,  DCfit131, DCfit141, DCfit122,  DCfit132, DCfit142, DCfit123,...
                DCfit133, DCfit143, DCfit124,  DCfit134, DCfit144, DCfit125,  DCfit135, DCfit145];
            
            %% computing the correlation coefficients between different channels of the sensors
            DCCorr= corrcoef(DC);
            DCCorr(isnan( DCCorr))=0;
            DCCorr=DCCorr+fvar.eps;
            DCCorrup= triu(DCCorr,1);
            DCCorrInd= find(abs(DCCorrup)>0);
            DCCorrV= DCCorrup(DCCorrInd)';
            
            %% standard deviation of channels
            DCstd= std(DC,0,1);
            %% skewness od data (3rd order moment)
            DCS= skewness(DC);
            %% Kurtosis of the dataset
            DCK= kurtosis(DC);
            %%
            %--------------------------------------------------------------------------------------------------------------------------------------
            %Computing the statistics of the derivative of data
            %--------------------------------------------------------------------------------------------------------------------------------------
            %%computing the  derivates (DSLN)
            DDC= diff(data{j},1);
            stamp_diff= diff(stamp{j});
            IDSL= (stamp_diff>0);
            DSLN= DDC(IDSL,:)./repmat(stamp_diff(IDSL),1,3);
            %% Computing the mean of derivatives
            DDCM= mean(DSLN,1);
            %% Computing the median of derivatives
            DDCMed= median(DSLN,1);
            %%  The real, imaginary and absolute value of the first 40 components   of fast fourier transform  of derivatives
            DDCFFT= fft(DSLN,40);
            DDCFFT_re= real(DDCFFT(:)');
            DDCFFT_im= imag(DDCFFT(:)');
            DDCFFT_abs= abs(DDCFFT(:)');
            %% Fitting the derivatives of recording signals as a function of time
            TD= stamp{j}(IDSL)-stamp{j}(ind(1));
%             bnum=100;
%             alp=1e-5;
%             muD= linspace(0,5,bnum);
%             varD=(20/bnum)^2;
%             dDsz=size(DSLN,1);
%             TDmat= repmat( TD, 1,bnum);
%             MDmat= repmat( muD, dDsz, 1);
%             Phi= [exp(-(TDmat-MDmat).^2./varD),ones(dDsz,1)];
%             WDmat= pinv(Phi'*Phi+alp*eye(bnum+1))*Phi'*DSLN;
%             WDvect=WDmat(:)';
            
            DDCfit121= polyfit(TD,DSLN(:,1), 1);
            DDCfit131= polyfit(TD,DSLN(:,2), 1);
            DDCfit141= polyfit(TD,DSLN(:,3), 1);
            DDCfit122= polyfit(TD,DSLN(:,1), 2);
            DDCfit132= polyfit(TD,DSLN(:,2), 2);
            DDCfit142= polyfit(TD,DSLN(:,3), 2);
            DDCfit123= polyfit(TD,DSLN(:,1), 3);
            DDCfit133= polyfit(TD,DSLN(:,2), 3);
            DDCfit143= polyfit(TD,DSLN(:,3), 3);
            DDCfit= [DDCfit121, DDCfit131,  DDCfit141, DDCfit122, DDCfit132,  DDCfit142,...
                DDCfit123, DDCfit133,  DDCfit143];
            
            %% correlation coefficients between different channels of the sensors
            DDCCorr= corrcoef(DSLN)+fvar.eps;
            DDCCorr(isnan( DDCCorr))=0;
            DDCCorr= DDCCorr+fvar.eps;
            DDCCorrup= triu(DDCCorr,1);
            DDCCorrInd= find(abs(DDCCorrup)>0);
            DDCCorrV= DDCCorrup(DDCCorrInd)';
            %%  standard deviations of derivatives
            DDCstd= std(DSLN,0,1);
            %% skewness of derivatives
            DDCS= skewness(DSLN);
            %% kurtosis of derivatives
            DDCK= kurtosis(DSLN);
            
            %% stacking the features of the  sensors (Gyro or accelerometer) together
            FC= [FC, DCM,  DCMed, DCFFT_re, DCFFT_im,  DCFFT_abs, DCfit,  DCCorrV, DCstd DCS, DCK, DDCM, DDCMed, DDCCorrV, DDCstd, DDCS, DDCK, DDCFFT_re, DDCFFT_im,  DDCFFT_abs, DDCfit];
            
            %%
        else
            %--------------------------------------------------------------------------------------------------------------------------------------
            %  Extracting  features for barometer
            %--------------------------------------------------------------------------------------------------------------------------------------
            continue % skip barometer features
            %% sorting the  data in terms of their time stamps
            [stamp{j}, ind]= sort(stamp{j},'ascend');
            data{j}= data{j}(ind, :);
            %% First column: Elapsed time from the first recording
            DC=[stamp{j}-stamp{j}(1),data{j}];
            %%  Mean
            DCM= mean(DC,1);
            %% Median
            DCMed=median(data{j},1);
            %% Correlation coefficients between the sensors and time
            DCCorr= corrcoef(DC);
            DCCorr(isnan( DCCorr))=0;
            DCCorr= DCCorr+fvar.eps;
            DCCorrup= triu(DCCorr,1);
            DCCorrInd= find(abs(DCCorrup)>0);
            DCCorrV= DCCorrup(DCCorrInd)';
            %% The real, imaginary and absolute value of the first 5 components   of fast fourier transform
            DCFFT= fft(data{j},5);
            DCFFT_re= real(DCFFT(:)');
            DCFFT_im= imag(DCFFT(:)');
            DCFFT_abs= abs(DCFFT(:)');
            
            %% Fitting the  recording as the function of time
            DCfit121= polyfit(DC(:,1),DC(:,2), 1);
            DCfit131= polyfit(DC(:,1),DC(:,3), 1);
            DCfit122= polyfit(DC(:,1),DC(:,2), 2);
            DCfit132= polyfit(DC(:,1),DC(:,3), 2);
            DCfit123= polyfit(DC(:,1),DC(:,2), 3);
            DCfit133= polyfit(DC(:,1),DC(:,3), 3);
            DCfit124= polyfit(DC(:,1),DC(:,2), 4);
            DCfit134= polyfit(DC(:,1),DC(:,3), 4);
    
                   
            DCfit=[DCfit121,  DCfit131, DCfit122,  DCfit132,  DCfit123,...
                DCfit133,  DCfit124,  DCfit134];
            %% standard deviation
            DCstd= std(DC,0,1);
            %% 
            DCskew= skewness(DC);
            %% Kurtosis of data
            DCkurt= kurtosis(DC);
            %% derivative of data
            DDC= diff(data{j},1);
            
            stamp_diff= diff(stamp{j});
            
            
            IDSL= (stamp_diff>0);
            %% checks whether there is enough data to compute the derivative
            if any(IDSL>0)
                DSLN= DDC(IDSL,:)./ repmat(stamp_diff(IDSL),1,2);
                %% mean
                DDCM= mean(DSLN,1);
                %% median
                DDCMed= median(DSLN,1);
                %% standard deviation
                DDCstd= std(DSLN,0,1);
            else
                %% mean
                DDCM= zeros(1,2);
                %% median
                DDCstd= zeros(1,2);
                %% standard deviation
                DDCMed=zeros(1,2);
            end;
            %% Stacking  all the features of barometer together
            FC= [FC, DCM, DCMed DCCorrV, DCfit, DCstd, DCskew, DCkurt, DDCM, DDCMed, DDCstd, DCFFT_re, DCFFT_im, DCFFT_abs]; 
            
        end;
        
    end;
    
    F= FC(fvar.nzstd);
    %% standardizing the data by subtracting the mean and dividing it by variance
    %% ( mean and variance are computed on the training set)
    FN= (F- fvar.mu)./ fvar.std;
else
    FN=[];
    
end;

return
