%Build new training dataset by loading all jittered clips, removing clips
%with less than 100 samples and resamples to a constant sampling frequency
clear all
files = dir('labels_plus_data_jittered_*');
data = struct();
for f = 1:length(files)
    disp(files(f).name)
    labels = load(files(f).name);
    labels = labels.labels;
    
    %remove empty clips
    ind1 = cellfun(@isempty, labels.acce);
    ind2 = cellfun(@isempty, labels.gyro);
    ind3 = cellfun(@isempty, labels.baro);
    inds = ind1+ind2+ind3;
    inds = logical(inds);
    labels.timestamp(inds) =[];
    labels.value(inds) = [];
    labels.subject(inds) = [];
    labels.acce(inds) = [];
    labels.baro(inds) = [];
    labels.gyro(inds) = [];
    labels = rmfield(labels,'text');
    acce_clipsize = cellfun(@length,labels.acce);
    gyro_clipsize = cellfun(@length,labels.gyro);
    baro_clipsize = cellfun(@length,labels.baro);
    ind1 = acce_clipsize < 100; %| acce_clipsize > 1000;
    ind2 = gyro_clipsize < 100; %| gyro_clipsize > 1000;
    ind3 = baro_clipsize < 10;
    inds = logical(ind1+ind2+ind3);
    
    %sort timestamps in each clip
    % cellfun(@(x) sortrows(x,1), labels.acce,'UniformOutput',false);
    labels.acce = cellfun(@sortrows, labels.acce,'UniformOutput',false);
    labels.gyro = cellfun(@sortrows, labels.gyro,'UniformOutput',false);
    labels.baro = cellfun(@sortrows, labels.baro,'UniformOutput',false);
    %% offset time to 0
    for i =1:length(labels.acce)
        labels.acce{i}(:,1) = labels.acce{i}(:,1)-labels.acce{i}(1,1);
        labels.gyro{i}(:,1) = labels.gyro{i}(:,1)-labels.gyro{i}(1,1);
        labels.baro{i}(:,1) = labels.baro{i}(:,1)-labels.baro{i}(1,1);
    end
    
    
    %     %% %distribution of clip sizes
    %     acce_clipsize = cellfun(@length,labels.acce);
    %     gyro_clipsize = cellfun(@length,labels.gyro);
    %     baro_clipsize = cellfun(@length,labels.baro);
    %     figure,
    %     subplot(311), histogram(acce_clipsize)
    %     subplot(312), histogram(gyro_clipsize)
    %     subplot(313), histogram(baro_clipsize)
    %
    
    %% Resample all clips to 50Hz and 6Hz
    
    labels.acce = cellfun(@(x) resampleclip(x,'a'),labels.acce,'UniformOutput',false);
    labels.gyro = cellfun(@(x) resampleclip(x,'a'),labels.gyro,'UniformOutput',false);
    labels.baro = cellfun(@(x) resampleclip(x,'b'),labels.baro,'UniformOutput',false);
    
    %% Aggregate data
    if f == 1
        data = labels;
    else
        data.timestamp = [data.timestamp labels.timestamp];
        data.value = [data.value labels.value];
        data.subject = [data.subject labels.subject];
        data.acce = [data.acce labels.acce];
        data.gyro = [data.gyro labels.gyro];
        data.baro = [data.baro labels.baro];
    end
    
    
end
data = rmfield(data,'jitterstd');

%% most subjects contain only a subsets of falls or fall like 
% Remove subjects with less than 4 types of events (fall and fall like)
labels = data;
% show distribution of data across subjects
us = unique(labels.subject);
labels.subjcode = zeros(length(labels.subject),1);
for s = 1:length(us)
    us(s);
    ind = find(strcmp(labels.subject,us(s)));
    labels.subjcode(ind) = s;
end
figure, histogram(labels.subjcode)

%show distribution of fall types across subjects
nact = [];
figure, hold on
for i = 1:max(labels.subjcode)
    ind = labels.subjcode == i;
    nact = [nact length(unique(labels.value(ind)))];
    nu = unique(labels.value(ind));
    plot(i*ones(length(nu),1),nu,'o')
end
figure, bar(nact)

%remove subjects with less than 4 types of events (fall and fall like)
subjremove = find(nact < 4);
for i = 1:length(subjremove)
    ind = labels.subjcode == subjremove(i);
    labels.timestamp(ind) = [];
    labels.value(ind) = [];
    labels.subject(ind) = [];
    labels.acce(ind) = [];
    labels.baro(ind) = [];
    labels.gyro(ind) = [];
    labels.subjcode(ind) = [];
end

%reindex subjects
us = unique(labels.subject);
labels.subjcode = zeros(length(labels.subject),1);
for s = 1:length(us)
    us(s);
    ind = find(strcmp(labels.subject,us(s)));
    labels.subjcode(ind) = s;
end
figure, histogram(labels.subjcode)
%show distribution of fall types across subjects
nact = [];
figure, hold on
for i = 1:max(labels.subjcode)
    ind = labels.subjcode == i;
    nact = [nact length(unique(labels.value(ind)))];
    nu = unique(labels.value(ind));
    plot(i*ones(length(nu),1),nu,'o')
end
figure, bar(nact)

data = labels;
save('DataFallsResampled.mat','data')

%% distribution of clip sizes
acce_clipsize = cellfun(@length,data.acce);
gyro_clipsize = cellfun(@length,data.gyro);
baro_clipsize = cellfun(@length,data.baro);
figure,
subplot(311), histogram(acce_clipsize)
subplot(312), histogram(gyro_clipsize)
subplot(313), histogram(baro_clipsize)


%% Visualization
% %% plot 5 random clips
% close all
% inds = randperm(length(labels.acce));
% inds = inds(1:5);
% x = labels.acce;
% % x = labels.gyro;
% for c = 1:length(inds)
%     figure, plot(x{inds(c)}(:,2:end))
% end
%
% %% plot random clips barom
% close all
% inds = randperm(length(labels.acce));
% inds = inds(1:5);
% x = labels.baro;
% for c = 1:length(inds)
%     figure, plot(x{inds(c)}(:,end))
% end
%
%
% %% test resample on accelerometer clips
% ind = acce_clipsize > 500; %clips with repeating data (logical)
% N = length(find(ind))
% % inds = find(ind); %clip indices
% %--plot clips
% acc = labels.acce(ind);
% clip = 3;
% figure, plot(acc{clip}(:,2:end))
% %plot all clips
% % for clip = 1:length(find(ind))
% %     figure, plot(acc{clip}(:,2:end))
% % end
%
% %Test Interpolate data to 250 samples
% %remove duplicated timestamps and associated data
% dt = diff(acc{clip}(:,1));
% indrep = dt == 0;
% acc{clip} = acc{clip}(~indrep,:);
% x = acc{clip};
%
% xres = resample(x(:,2:end),x(:,1),50);
% xres = interp1(x(:,1),x(:,2:end),linspace(0,5,250),'spline');
% figure, plot(xres)
% %%
%
% ind = baro_clipsize > 40; %clips with repeating data (logical)
% N = length(find(ind))
% % inds = find(ind); %clip indices
% %--plot clips
% b = labels.baro(ind);
% x = b{1};
%
% sf = 6;
% dt = diff(x(:,1));
% indrep = dt == 0;
% x = x(~indrep,:);
% % xr = resample(x(:,2:end),x(:,1),sf);
% t = linspace(0,x(end,1),x(end,1)*sf);
% xr = interp1(x(:,1),x(:,2:end),t,'nearest')
% xres = [t' xr];
%
% % xres = resampleclip(x,6)
% figure, plot(x(:,1),x(:,end),'o-',t,xres(:,end),'o-')
