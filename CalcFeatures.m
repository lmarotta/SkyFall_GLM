%% CalcFeatures

function FC=CalcFeatures(data,stamp)
fvar.eps=1e-6; %threshold to prevent correlation coefficients from being 0
nbins_large=50;
nbins=10; %for FFT
nwin=10;
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
        %% standard deviation of channels
        DCstd= std(data{j},0,1);
        %% skewness od data (3rd order moment)
        DCS= skewness(data{j});
        %% Kurtosis of the dataset
        DCK= kurtosis(data{j});
        %% Entropy
        DCENT=Entropy(floor(zscore(data{j})/.2));
        
        %%  the real, imaginary and absolute value of the first 40 components   of fast fourier transform
        
        DCFFT= fft(data{j});
        DCFFT_re=[];
        DCFFT_im=[];
        DCFFT_abs=[];
        for ii=1:nbins_large/2
            for jj=1:3
                DCFFT_re= [DCFFT_re trapz(real(DCFFT(floor(length(DCFFT)/nbins_large*(ii-1)+1):floor(length(DCFFT)/nbins_large*ii),jj)'))];
                DCFFT_im= [DCFFT_im trapz(imag(DCFFT(floor(length(DCFFT)/nbins_large*(ii-1)+1):floor(length(DCFFT)/nbins_large*ii),jj)'))];
                DCFFT_abs= [DCFFT_abs trapz(abs(DCFFT(floor(length(DCFFT)/nbins_large*(ii-1)+1):floor(length(DCFFT)/nbins_large*ii),jj)'))];
            end
        end
        
        FFTENT=Entropy(floor(zscore(abs(DCFFT))/.2));
        
        DC_E=[];
        E=[];
        
        for ii=1:nwin-(nwin/10-1)
            winsize=floor(length(data{j})/10);
            d=data{j}(floor(length(data{j})/nwin*(ii-1)+1):floor(length(data{j})/nwin*(ii-1)+winsize),:);
            % m=mean(d);
            m=zeros(1, size(d,2));
            DCFFT=fft(d-repmat(m,[size(d,1) 1]));
            E_temp=[];
            for jj=1:5
                for kk=1:3 %axis of sensor
                    temp(1,1,kk)=trapz(abs(DCFFT(floor(length(DCFFT)/nbins*(jj-1)+1):floor(length(DCFFT)/nbins*jj),kk)'));
                    
                    if mod(ii-1,nwin/10)==0
                        DCFFT_abs=[DCFFT_abs temp(kk)];
                    end  
                end
                E_temp=[E_temp temp];
            end
            DC_E=[DC_E; E_temp];
            E=[E; sum(d.^2)];
        end
        
        Ent_temp=[];
        
        for ii=1:5
            for jj=1:3
                Ent_temp(:,ii,jj)=Entropy(floor(zscore(DC_E(:,ii,jj))/.2));
            end
        end
        
        DCENERGY=[mean(DC_E) std(DC_E) skewness(DC_E) kurtosis(DC_E) Ent_temp];
        DCENERGY=DCENERGY(:)';
%         DCENERGY=[DCENERGY(:)' mean(E) std(E) skewness(E) kurtosis(E) Entropy(floor(zscore(E)/.2))];
        
        %vector magnitude (norm) and max features
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
        
        DCAMEAN=[mean(Angle1) mean(Angle2) mean(Angle3)];
        DCARANGE=[range(Angle1) range(Angle2) range(Angle3)];
        DCAIQR=[iqr(Angle1) iqr(Angle2) iqr(Angle3)];
        DCAMAX=[max(Angle1) max(Angle2) max(Angle3)];
        DCAMIN=[min(Angle1) min(Angle2) min(Angle3)];
        DCASTD=[std(Angle1) std(Angle2) std(Angle3) ...
            skewness(Angle1) skewness(Angle2) skewness(Angle3) ...
            kurtosis(Angle1) kurtosis(Angle2) kurtosis(Angle3)];
              
        % SampEntropy
        
        DCENT=[DCENT Entropy(floor(zscore(data{j})/.2))];
        
        % X Products
        XY=data{j}(:,1).*data{j}(:,2);
        XZ=data{j}(:,1).*data{j}(:,3);
        YZ=data{j}(:,2).*data{j}(:,3);
        
        XProd=[XY XZ YZ];
        
        XP=[abs(mean(XProd)) std(XProd) skewness(XProd) kurtosis(XProd) median(XProd) iqr(XProd) range(XProd) max(XProd) min(XProd)];
        
        %Autocorr Features
               
        X=[xcorr(data{j}(:,1)) xcorr(data{j}(:,2)) xcorr(data{j}(:,3))];
        XM=mean(X,1);
        XSD=[std(X) skewness(X) kurtosis(X)];
        XMed=[median(X) iqr(X) range(X) max(X) min(X)];
        
        NEWFEAT=[DCRM DCRMed DCRSD DCMAX DCRANGE DCIQR DCAMEAN DCARANGE DCAIQR DCAMAX DCAMIN DCASTD XM XSD XMed];    
        
        %% computing the correlation coefficients between different channels of the sensors
        DCCorr= corrcoef(DC);
        DCCorr(isnan( DCCorr))=0;
        DCCorr=DCCorr+fvar.eps;
        DCCorrup= triu(DCCorr,1);
        DCCorrInd= find(abs(DCCorrup)>0);
        DCCorrV= DCCorrup(DCCorrInd)';
        %last line could be updated with DCCorrV= DCCorrup
        
        
        %%
        %--------------------------------------------------------------------------------------------------------------------------------------
        %Computing the statistics of the derivative of data
        %--------------------------------------------------------------------------------------------------------------------------------------
        %%computing the  derivates (DSLN)
        DDC= diff(data{j},1);
        stamp_diff= diff(stamp{j});
        IDSL= (stamp_diff>0);
        DSLN= DDC(IDSL,:)./repmat(stamp_diff(IDSL),1,3);
        stamp_diff=stamp_diff(IDSL);
        %% Computing the mean of derivatives
        DDCM= mean(DSLN,1);
        %% Computing the median of derivatives
        DDCMed= median(DSLN,1);
        %%  The real, imaginary and absolute value of the first 40 components   of fast fourier transform  of derivatives
        
        DDCFFT= fft(DSLN);
        DDCFFT_re=[];
        DDCFFT_im=[];
        DDCFFT_abs=[];
        for ii=1:nbins_large/2
            for jj=1:3
                DDCFFT_abs= [DDCFFT_abs trapz(abs(DDCFFT(floor(length(DCFFT)/nbins_large*(ii-1)+1):floor(length(DDCFFT)/nbins_large*ii),jj)'))];
            end
        end
        
        DFFTENT=Entropy(floor(zscore(abs(DDCFFT))/.2));
        
        FFT_E=[];
        E=[];
        
        for ii=1:nwin-(nwin/10-1)
            winsize=floor(length(data{j})/10);
            d=data{j}(floor(length(data{j})/nwin*(ii-1)+1):floor(length(data{j})/nwin*(ii-1)+winsize),:);
            %                     m=mean(d);
            m=zeros(1, size(d,2));
            DDCFFT=fft(d-repmat(m,[size(d,1) 1]));
           
            E_temp=[];
            for jj=1:5
                for kk=1:3
                    temp(1,1,kk)=trapz(abs(DDCFFT(floor(length(DDCFFT)/nbins*(jj-1)+1):floor(length(DDCFFT)/nbins*jj),kk)'));
                    
                    if mod(ii-1,nwin/10)==0
                        DCFFT_abs=[DCFFT_abs temp(kk)];
                    end
                end
                E_temp=[E_temp temp];
            end
            FFT_E=[FFT_E; E_temp];
            E=[E; sum(d.^2)];
        end
        
        Ent_temp=[];
        
        for ii=1:5
            for jj=1:3
                Ent_temp(:,ii,jj)=Entropy(floor(zscore(FFT_E(:,ii,jj))/.2));
            end
        end
        
        FFTENERGY=[mean(FFT_E) std(FFT_E) skewness(FFT_E) kurtosis(FFT_E) Ent_temp];
        FFTENERGY=[FFTENERGY(:)' mean(E) std(E) skewness(E) kurtosis(E) Entropy(floor(zscore(E)/.2))];
        
        % SampEntropy
        
        DCENT=[DCENT Entropy(floor(zscore(DDC)/.2))];
        
        % X Products
        XY=DSLN(:,1).*DSLN(:,2);
        XZ=DSLN(:,1).*DSLN(:,3);
        YZ=DSLN(:,2).*DSLN(:,3);
        
        XProd=[XY XZ YZ];
        
        XP=[XP abs(mean(XProd)) std(XProd) skewness(XProd) kurtosis(XProd) median(XProd) iqr(XProd) range(XProd) max(XProd) min(XProd)];

        %%  standard deviations of derivatives
        DDCstd= std(DSLN,0,1);
        %% skewness of derivatives
        DDCS= skewness(DSLN);
        %% kurtosis of derivatives
        DDCK= kurtosis(DSLN);
        
        %% stacking the features of the  sensors (Gyro or accelerometer) together
        FC= [FC, DCM, DCMed, DCstd, DCS, DCK, DCCorrV, DCFFT_abs, DDCM, DDCMed, DDCstd, DDCS, DDCK, DDCFFT_abs, NEWFEAT, FFTENT, DFFTENT, DCENT, XP, DCENERGY, FFTENERGY];
        
    else
        %--------------------------------------------------------------------------------------------------------------------------------------
        %  Extracting  features for barometer
        %--------------------------------------------------------------------------------------------------------------------------------------
        %             continue % skip barometer features
        %% sorting the  data in terms of their time stamps
        [stamp{j}, ind]= sort(stamp{j},'ascend');
        data{j}= data{j}(ind, 1);
        data{j}=data{j}-ones(length(data{j}),1)*mean(data{j});
        %% First column: Elapsed time from the first recording
        DC=[stamp{j}-stamp{j}(1),data{j}]; % use only pressure (no altitude)

        %% Correlation coefficients between the sensors and time
        DCCorr= corrcoef(DC);
        DCCorr(isnan( DCCorr))=0;
        DCCorr= DCCorr+fvar.eps;
        DCCorrup= triu(DCCorr,1);
        DCCorrInd= find(abs(DCCorrup)>0);
        DCCorrV= DCCorrup(DCCorrInd)';
        %% The real, imaginary and absolute value of the first 5 components   of fast fourier transform
        
        DCFFT= fft(data{j});
        DCFFT_abs=[];

        for ii=1:5
            d=data{j}(floor(length(data{j})/5*(ii-1)+1):floor(length(data{j})/5*ii),:);
            %                     m=mean(d);
            m=zeros(1, size(d,2));
            DCFFT=fft(d-repmat(m,[size(d,1) 1]));
            for kk=1:1
                DCFFT_abs=[DCFFT_abs trapz(abs(DCFFT(:,kk)'))];
            end
        end
        
        %% standard deviation
        DCstd= std(data{j},0,1);

        %% Stacking  all the features of barometer together
        FC= [FC, DCCorrV, DCstd, DCFFT_abs];
        
    end;
    
end;