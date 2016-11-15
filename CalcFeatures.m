%% CalcFeatures

function FC=CalcFeatures(data,stamp)
fvar.eps=1e-6; %threshold to prevent correlation coefficients from being 0
nbins=8; %for FFT
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
        DCM= mean(data{j},1);
        %% median
        DCMed= median(data{j},1);
        %%  the real, imaginary and absolute value of the first 40 components   of fast fourier transform
        
        DCFFT= fft(data{j});
        DCFFT_re=[];
        DCFFT_im=[];
        DCFFT_abs=[];
        for ii=1:40
            for jj=1:3
                DCFFT_re= [DCFFT_re trapz(real(DCFFT(floor(length(DCFFT)/40*(ii-1)+1):floor(length(DCFFT)/40*ii),jj)'))];
                DCFFT_im= [DCFFT_im trapz(imag(DCFFT(floor(length(DCFFT)/40*(ii-1)+1):floor(length(DCFFT)/40*ii),jj)'))];
                DCFFT_abs= [DCFFT_abs trapz(abs(DCFFT(floor(length(DCFFT)/40*(ii-1)+1):floor(length(DCFFT)/40*ii),jj)'))];
            end
        end
        for ii=1:5
            d=data{j}(floor(length(data{j})/5*(ii-1)+1):floor(length(data{j})/5*ii),:);
            %                     m=mean(d);
            m=zeros(1, size(d,2));
            DCFFT=fft(d-repmat(m,[size(d,1) 1]));
            for jj=1:nbins
                for kk=1:3 %axis of sensor
                    DCFFT_re=[DCFFT_re trapz(real(DCFFT(floor(length(DCFFT)/nbins*(jj-1)+1):floor(length(DCFFT)/nbins*jj),kk)'))];
                    DCFFT_im=[DCFFT_im trapz(imag(DCFFT(floor(length(DCFFT)/nbins*(jj-1)+1):floor(length(DCFFT)/nbins*jj),kk)'))];
                    DCFFT_abs=[DCFFT_abs trapz(abs(DCFFT(floor(length(DCFFT)/nbins*(jj-1)+1):floor(length(DCFFT)/nbins*jj),kk)'))];
                end
            end
        end
        
        
        
        % resultant and max features
        res=sum(data{j}.^2,2);
        DCRM=mean(res);
        DCRMed=median(res);
        DCRSD=[std(res) skewness(res) kurtosis(res)];
        DCMAX=[max(abs(data{j})) max(res) max(diff(res)) min(diff(res))];
        DCIQR=[iqr(data{j}) iqr(res) iqr(diff(res))];
        
        DCRANGE=range(data{j});
        
        Angle1=atan2(data{j}(:,1),data{j}(:,2));
        Angle2=atan2(data{j}(:,1),data{j}(:,3));
        Angle3=atan2(data{j}(:,2),data{j}(:,3));
        
        % low-pass filter on angles
        
        %                [B,A]=butter(1,5/125);
        %
        %                Angle1=filtfilt(B,A,Angle1);
        %                Angle2=filtfilt(B,A,Angle2);
        %                Angle3=filtfilt(B,A,Angle3);
        
        DCAMEAN=[mean(Angle1) mean(Angle2) mean(Angle3)];
        DCARANGE=[range(Angle1) range(Angle2) range(Angle3)];
        DCAIQR=[iqr(Angle1) iqr(Angle2) iqr(Angle3)];
        DCAMAX=[max(Angle1) max(Angle2) max(Angle3)];
        DCAMIN=[min(Angle1) min(Angle2) min(Angle3)];
        DCASTD=[std(Angle1) std(Angle2) std(Angle3) ...
            skewness(Angle1) skewness(Angle2) skewness(Angle3) ...
            kurtosis(Angle1) kurtosis(Angle2) kurtosis(Angle3)];
        
        %% Change in orientation after impulse (0 if no impulse)
        [M,I]=max(res.^.5);
        if M>2*9.8 && j==2 % threshold in g's
            indstart=max(I-50,1);
            indend=min(I+50,length(Angle1));
            
            d1=Angle1(indend)-Angle1(indstart);
            d2=Angle2(indend)-Angle2(indstart);
            d3=Angle3(indend)-Angle3(indstart);
            dall=dot(data{j}(indend),data{j}(indstart));
            
            DCADELTA=[abs(d1) abs(d2) abs(d3) dall dall/norm(dall)];
        else
            DCADELTA=[0 0 0 0 0];
        end
        
        %Autocorr Features
        
        X=xcorr(data{j});
        XM=mean(X,1);
        XSD=[std(X) skewness(X) kurtosis(X)];
        XMed=[median(X) iqr(X) range(X) max(X) min(X)];
        
        X=diff(X);
        XDM=mean(X,1);
        XDSD=[std(X) skewness(X) kurtosis(X)];
        XDMed=[median(X) iqr(X) range(X) max(X) min(X)];
        
        NEWFEAT=[DCRM DCRMed DCRSD DCMAX DCRANGE DCIQR DCAMEAN DCARANGE DCAIQR DCAMAX DCAMIN DCASTD DCADELTA XM XSD XMed XDM XDSD XDMed];
        
        
        
        
        
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
        %last line could be updated with DCCorrV= DCCorrup
        
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
        
        DDCFFT= fft(DSLN);
        DDCFFT_re=[];
        DDCFFT_im=[];
        DDCFFT_abs=[];
        for ii=1:40
            for jj=1:3
                DDCFFT_re= [DDCFFT_re trapz(real(DDCFFT(floor(length(DCFFT)/40*(ii-1)+1):floor(length(DDCFFT)/40*ii),jj)'))];
                DDCFFT_im= [DDCFFT_im trapz(imag(DDCFFT(floor(length(DCFFT)/40*(ii-1)+1):floor(length(DDCFFT)/40*ii),jj)'))];
                DDCFFT_abs= [DDCFFT_abs trapz(abs(DDCFFT(floor(length(DCFFT)/40*(ii-1)+1):floor(length(DDCFFT)/40*ii),jj)'))];
            end
        end
        for ii=1:5
            d=DSLN(floor(length(DSLN)/5*(ii-1)+1):floor(length(DSLN)/5*ii),:);
            %                     m=mean(d);
            m=zeros(1, size(d,2));
            DDCFFT=fft(d-repmat(m,[size(d,1) 1]));
            for jj=1:nbins
                for kk=1:3
                    DDCFFT_re=[DDCFFT_re trapz(real(DDCFFT(floor(length(DDCFFT)/nbins*(jj-1)+1):floor(length(DDCFFT)/nbins*jj),kk)'))];
                    DDCFFT_im=[DDCFFT_im trapz(imag(DDCFFT(floor(length(DDCFFT)/nbins*(jj-1)+1):floor(length(DDCFFT)/nbins*jj),kk)'))];
                    DDCFFT_abs=[DDCFFT_abs trapz(abs(DDCFFT(floor(length(DDCFFT)/nbins*(jj-1)+1):floor(length(DDCFFT)/nbins*jj),kk)'))];
                end
            end
        end
        
        
        %% Fitting the derivatives of recording signals as a function of time
        TD= stamp{j}(IDSL)-stamp{j}(ind(1));
        
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
        FC= [FC, DCM,  DCMed, DCFFT_re, DCFFT_im,  DCFFT_abs, DCfit,  DCCorrV, DCstd DCS, DCK, DDCM, DDCMed, DDCCorrV, DDCstd, DDCS, DDCK, DDCFFT_re, DDCFFT_im,  DDCFFT_abs, DDCfit, NEWFEAT];
        
        %%
    else
        %--------------------------------------------------------------------------------------------------------------------------------------
        %  Extracting  features for barometer
        %--------------------------------------------------------------------------------------------------------------------------------------
        %             continue % skip barometer features
        %% sorting the  data in terms of their time stamps
        [stamp{j}, ind]= sort(stamp{j},'ascend');
        data{j}= data{j}(ind, :);
        %% First column: Elapsed time from the first recording
        DC=[stamp{j}-stamp{j}(1),data{j}];
        %             DC=[stamp{j}-stamp{j}(1),data{j}-mean(data{j})];
        %%  Mean
        DCM= mean(data{j},1);
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
        
        DCFFT= fft(data{j});
        DCFFT_re=[];
        DCFFT_im=[];
        DCFFT_abs=[];
        %FFT over 5 sec clips
        for ii=1:5
            for jj=1:2
                DCFFT_re= [DCFFT_re trapz(real(DCFFT(floor(length(DCFFT)/40*(ii-1)+1):floor(length(DCFFT)/40*ii),jj)'))];
                DCFFT_im= [DCFFT_im trapz(imag(DCFFT(floor(length(DCFFT)/40*(ii-1)+1):floor(length(DCFFT)/40*ii),jj)'))];
                DCFFT_abs= [DCFFT_abs trapz(abs(DCFFT(floor(length(DCFFT)/40*(ii-1)+1):floor(length(DCFFT)/40*ii),jj)'))];
            end
        end
        %FFT over 1 sec windows for the 5 sec clips
        for ii=1:5
            d=data{j}(floor(length(data{j})/5*(ii-1)+1):floor(length(data{j})/5*ii),:);
            %                     m=mean(d);
            m=zeros(1, size(d,2));
            DCFFT=fft(d-repmat(m,[size(d,1) 1]));
            for kk=1:2
                DCFFT_re=[DCFFT_re trapz(real(DCFFT(:,kk)'))];
                DCFFT_im=[DCFFT_im trapz(imag(DCFFT(:,kk)'))];
                DCFFT_abs=[DCFFT_abs trapz(abs(DCFFT(:,kk)'))];
            end
        end
        
        
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