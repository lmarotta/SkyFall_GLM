clear
close all

load('RealWorldFallsJuly.mat');

Ts=0.02;
fc1=0.5; %to remove low frequency components (falls are charactherized in higher frequency domain)
nfc1=2*Ts*fc1; % normalized freq: nfc=fc/fnyq
[B,A]=butter(4,nfc1,'high'); % second order

% extractfeatures
for i=1:length(mylabels.acce)
    
    %acc_norm and gyro
    z0=cell2mat(mylabels.acce(i,1));
    z1=filtfilt(B,A,z0);
    z2=sqrt((z0(:,2).^2)+(z0(:,3).^2)+(z0(:,4).^2));
    z=filtfilt(B,A,z2);
    
    
    p1=cell2mat(mylabels.gyro(i,1));
    p=sqrt((p1(:,2).^2)+(p1(:,3).^2)+(p1(:,4).^2));
    
    %x,y,z
    b1=z1(:,2);
    b2=z1(:,3);
    b3=z1(:,4);
    x_val=mean(b1);
    y_val=mean(b2);
    z_val=mean(b3);
    
    
    d1=p1(:,2);
    d2=p1(:,3);
    d3=p1(:,4);
    
    %Variance and Angle
    
    myVariance(i,:)=var(z);

    %angle bins
    [my_max,my_ind]=max(abs(z));
    new_b1=b1(my_ind:end);
    new_b2=b2(my_ind:end);
    new_b3=b3(my_ind:end);
    
    first_b1=b1(1:my_ind);
    first_b2=b2(1:my_ind);
    first_b3=b3(1:my_ind);
    
    x_val1=mean(new_b1);
    y_val1=mean(new_b2);
    z_val1=mean(new_b3);
    
    x_val2=mean(first_b1);
    y_val2=mean(first_b2);
    z_val2=mean(first_b3);
    
    vectcross=[x_val1; y_val1; z_val1];
    vect=[x_val2; y_val2; z_val2];
    
    angolino=dot(vectcross,vect)./(norm(vectcross).*norm(vect));
    myAngle(i,:)=acosd(angolino)';
    %stat_featurea
    
    k_a(i,:)=kurtosis(z);
    s_a(i,:)=skewness(z);
    m1_a(i,:)=mean(z2);
    m2_a(i,:)=median(z2);
    i_r_a(i,:)=iqr(z2);
    s_d_a(i,:)=std(z);
    
    k_g(i,:)=kurtosis(p);
    s_g(i,:)=skewness(p);
    m1_g(i,:)=mean(p);
    m2_g(i,:)=median(p);
    i_r_g(i,:)=iqr(p);
    s_d_g(i,:)=std(p);
    
    kx_8(i,:)=kurtosis(b1);
    sx_8(i,:)=skewness(b1);
    i_rx_8(i,:)=iqr(b1);    
    m1x_8(i,:)=mean(b1);
    m2x_8(i,:)=median(b1);
    s_dx_8(i,:)=std(b1);
    
    ky_8(i,:)=kurtosis(b2);
    sy_8(i,:)=skewness(b2);
    i_ry_8(i,:)=iqr(b2);
    m1y_8(i,:)=mean(b2);
    m2y_8(i,:)=median(b2);
    s_dy_8(i,:)=std(b2);
    
    
    kz_8(i,:)=kurtosis(b3);
    sz_8(i,:)=skewness(b3);
    i_rz_8(i,:)=iqr(b3);
    m1z_8(i,:)=mean(b3);
    m2z_8(i,:)=median(b3);
    s_dz_8(i,:)=std(b3);
    
    kx_0(i,:)=kurtosis(d1);
    sx_0(i,:)=skewness(d1);
    i_rx_0(i,:)=iqr(d1);
    m1x_0(i,:)=mean(d1);
    m2x_0(i,:)=median(d1);
    s_dx_0(i,:)=std(d1);
   
    ky_0(i,:)=kurtosis(d2);
    sy_0(i,:)=skewness(d2);
    i_ry_0(i,:)=iqr(d2);
    m1y_0(i,:)=mean(d2);
    m2y_0(i,:)=median(d2);
    s_dy_0(i,:)=std(d2);
    
    kz_0(i,:)=kurtosis(d3);
    sz_0(i,:)=skewness(d3);
    i_rz_0(i,:)=iqr(d3);
    m1z_0(i,:)=mean(d3);
    m2z_0(i,:)=median(d3);
    s_dz_0(i,:)=std(d3);
    
    %energy and rms
    rootms(i,:)=rms(z);
    energiax(i,:)=sum(b1.^2);
    energiay(i,:)=sum(b2.^2);
    energiaz(i,:)=sum(b3.^2);
    energiag(i,:)=sum(p.^2);
    
    %derivative_features
    der=diff(z);
    steepa0(i,:)= min(der);
    
    der1=diff(b1);
    vertx8(i,:)= min(der1);
    
    der2=diff(b2);
    verty8(i,:)= min(der2);
     
    der3=diff(b3);
    vertz8(i,:)= min(der3);
    
    der4=diff(p);
    steepg0(i,:)= min(der4);
    
    %frequency domain features
    
    [Sw,fw]=pwelch(z,[],[],[],50);
    [m(i,:),ind]=max(Sw);
    sk(i,:)=skewness(Sw);
    kr(i,:)=kurtosis(Sw);
    s_e(i,:)=sampleEntropy(z,1,0.2*std(z),1);
    
    %peak frequency
    max_f(i,:)=fw(ind);
    
    %orientation
    accx = atand(b2./b3);
    accy = atand(b1./b3);
    accz = atand(b1./b2);
    
    m_ox(i,:) = max(abs(accx));
    v_ox(i,:) = var(accx);
    
    m_oy(i,:) = max(abs(accy));
    v_oy(i,:) = var(accy);
    
    m_oz(i,:) = max(abs(accz));
    v_oz(i,:) = var(accz);
    
end




% put features in table
d=mylabels.phone;
RealFalls=table(d,myVariance,myAngle,k_a,s_a,s_d_a,m1_a,m2_a,i_r_a,k_g,s_g,s_d_g,m1_g,m2_g,i_r_g,rootms,energiax,energiay,energiaz,energiag,steepa0,steepg0,vertx8,verty8,vertz8,kx_8,sx_8,i_rx_8,ky_8,sy_8,i_ry_8,kz_8,sz_8,i_rz_8,kx_0,sx_0,i_rx_0,ky_0,sy_0,i_ry_0,kz_0,sz_0,i_rz_0,max_f,m,sk,kr,s_e,m1x_8,m2x_8,s_dx_8,m1y_8,m2y_8,s_dy_8,m1z_8,m2z_8,s_dz_8,m1x_0,m2x_0,s_dx_0,m1y_0,m2y_0,s_dy_0,m1z_0,m2z_0,s_dz_0,m_ox,v_ox,m_oy,v_oy,m_oz,v_oz);
a1={'Subject','Variance','Angle','Kurtosis_a','Skewness_a','SD_a','Mean_a','Median_a','IQR_a','Kurtosis_g','Skewness_g','SD_g','Mean_g','Median_g','IQR_g','RMS','EnergyX','EnergyY','EnergyZ','EnergyG','Acc_steepness_afterpeak','Gyro_steepness_afterpeak','Acc_steepness_afterpeak_X','Acc_steepness_afterpeak_Y','Acc_steepness_afterpeak_Z','Kurtosis_x_a','Skewness_x_a','IQR_x_a','Kurtosis_y_a','Skewness_y_a','IQR_y_a','Kurtosis_z_a','Skewness_z_a','IQR_z_a','Kurtosis_x_g','Skewness_x_g','IQR_x_g','Kurtosis_y_g','Skewness_y_g','IQR_y_g','Kurtosis_z_g','Skewness_z_g','IQR_z_g','Max_f','Periodogram_maxf','Skewness_fft','Kurtosis_fft','S_Entropy','NF1','NF2','NF3','NF4','NF5','NF6','NF7','NF8','NF9','NF10','NF11','NF12','NF13','NF14','NF15','NF16','NF17','NF18','maxOrienx','varOrienx','maxOrieny','varOrieny','maxOrienz','varOrienz'};
RealFalls.Properties.VariableNames = a1; 