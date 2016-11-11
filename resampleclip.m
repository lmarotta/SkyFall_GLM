function xres = resampleclip(x,sensor)

    dt = diff(x(:,1));     %remove duplicated timestamps and associated data
    indrep = dt == 0;
    x = x(~indrep,:); 
    if strcmp(sensor,'a') %accelerometer and gyro
        sf = 50;    %resampling freq [Hz]
        xr = resample(x(:,2:end),x(:,1),sf); %resample at 50 Hz and filter signal
        t = linspace(0,length(xr)/sf,length(xr));
    else
        sf = 6;
        t = linspace(0,x(end,1),x(end,1)*sf);
        xr = interp1(x(:,1),x(:,2:end),t,'nearest');
    end
    xres = [t' xr];
    

 