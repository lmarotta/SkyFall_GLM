%%% create jittered version of the original (10 sec long) fall clips
function JitterClips(allclips)

%% Create clips normally jittered w different std dev
stdwin = [.1 .5 1 1.5 2 2.5];
strwin = {'01','05','1','1_5','2','2_5'};

for l = 1:length(stdwin)
    clearvars -except stdwin strwin l allclips
    load labels_plus_dataLONG.mat
    
    if ~allclips
        load indsjit.mat %use test clips (25% of data)
        
        labels.timestamp = labels.timestamp(indsjit);
        labels.value = labels.value(indsjit);
        labels.subject = labels.subject(indsjit);
        labels.text = labels.text(indsjit);
        labels.acce = labels.acce(indsjit);
        labels.gyro = labels.gyro(indsjit);
        labels.baro = labels.baro(indsjit);      
    end

    labels.jitterstd = stdwin(l);   

    %remove clips with too few datapoints
    sizeacc = cellfun(@(x) size(x,1)>100,labels.acce);
    sizegyro = cellfun(@(x) size(x,1)>100,labels.gyro);
    sizebaro = cellfun(@(x) size(x,1)>10,labels.baro);
    keep = sizeacc & sizegyro & sizebaro;
    labels.timestamp = labels.timestamp(keep);
    labels.value = labels.value(keep);
    labels.subject = labels.subject(keep);
    labels.text = labels.text(keep);
    labels.acce = labels.acce(keep);
    labels.gyro = labels.gyro(keep);
    labels.baro = labels.baro(keep);
    
    %measure duration of each clip
    for k = 1:length(labels.acce)
        if ~isempty(labels.acce{k})
            dur(k) = max(labels.acce{k}(:,1))-min(labels.acce{k}(:,1));
        end
    end
    
    %delete clips < 7.5 secs
    id = dur>7.5;
    labels.timestamp = labels.timestamp(id);
    labels.value = labels.value(id);
    labels.subject = labels.subject(id);
    labels.text = labels.text(id);
    labels.acce = labels.acce(id);
    labels.gyro = labels.gyro(id);
    labels.baro = labels.baro(id);
    dur = dur(id)';
    
    %jitter
    N = length(labels.timestamp);   %# of clips
    starttimes = 2.5 + max(min(stdwin(l)*randn(N,1), 2.5),-2.5);  %normal with std stdwin sec
    endtimes = min(starttimes+5,dur);
    
    %find clips less than 5 secs long
    durchop = endtimes-starttimes;
    idshort = durchop < 5;
    while any(idshort)
        n = sum(idshort);
        starttimes(idshort) = 2.5 + max(min(stdwin(l)*randn(n,1),2.5),-2.5);  %normal with std 1 sec
        endtimes(idshort) =  min(starttimes(idshort)+5,dur(idshort));
        durchop = endtimes-starttimes;
        idshort = durchop < 5;
    end
    
    %cutting the clips
    starttimestamp = zeros(length(labels.acce),1);
    for c = 1:length(labels.acce)
        labels.acce{c} = sortrows(labels.acce{c});
        starttimestamp(c) = labels.acce{c}(1,1);
        labels.acce{c}(:,1) = labels.acce{c}(:,1)-starttimestamp(c);
        timediff = abs(labels.acce{c}(:,1)-starttimes(c));
        [~,indmin] = min(timediff);
        timediff2 = abs(labels.acce{c}(:,1)-endtimes(c));
        [~,indmin2] = min(timediff2);
        
        labels.acce{c} = labels.acce{c}(indmin:indmin2,:);
    end
    for c = 1:length(labels.gyro)
        labels.gyro{c} = sortrows(labels.gyro{c});
        labels.gyro{c}(:,1) = labels.gyro{c}(:,1)-starttimestamp(c);
        timediff = abs(labels.gyro{c}(:,1)-starttimes(c));
        [~,indmin] = min(timediff);
        timediff2 = abs(labels.gyro{c}(:,1)-endtimes(c));
        [~,indmin2] = min(timediff2);
        
        labels.gyro{c} = labels.gyro{c}(indmin:indmin2,:);
    end
    for c = 1:length(labels.baro)
        labels.baro{c} = sortrows(labels.baro{c});
        labels.baro{c}(:,1) = labels.baro{c}(:,1)-starttimestamp(c);
        timediff = abs(labels.baro{c}(:,1)-starttimes(c));
        [~,indmin] = min(timediff);
        timediff2 = abs(labels.baro{c}(:,1)-endtimes(c));
        [~,indmin2] = min(timediff2);
    
        labels.baro{c} = labels.baro{c}(indmin:indmin2,:);
    end
    
%     figure, histogram(starttimes)
    
    %remove clips with too few datapoints
    sizeacc = cellfun(@(x) size(x,1)>100,labels.acce);
    sizegyro = cellfun(@(x) size(x,1)>100,labels.gyro);
    sizebaro = cellfun(@(x) size(x,1)>10,labels.baro);
    keep = sizeacc & sizegyro & sizebaro;
    labels.timestamp = labels.timestamp(keep);
    labels.value = labels.value(keep);
    labels.subject = labels.subject(keep);
    labels.text = labels.text(keep);
    labels.acce = labels.acce(keep);
    labels.gyro = labels.gyro(keep);
    labels.baro = labels.baro(keep);
    
    filename = ['labels_plus_data_jittered_', strwin{l}, 's.mat'];
    disp(['saving ', filename])
    save(filename, 'labels')
    
end

%% Create Uniform jitter

    clearvars -except stdwin strwin l allclips
    load labels_plus_dataLONG.mat
    
    if ~allclips
        load indsjit.mat %use only 25% of original datasets (for testing)
        
        labels.timestamp = labels.timestamp(indsjit);
        labels.value = labels.value(indsjit);
        labels.subject = labels.subject(indsjit);
        labels.text = labels.text(indsjit);
        labels.acce = labels.acce(indsjit);
        labels.gyro = labels.gyro(indsjit);
        labels.baro = labels.baro(indsjit);      
    end   

    %remove clips with too few datapoints
    sizeacc = cellfun(@(x) size(x,1)>100,labels.acce);
    sizegyro = cellfun(@(x) size(x,1)>100,labels.gyro);
    sizebaro = cellfun(@(x) size(x,1)>10,labels.baro);
    keep = sizeacc & sizegyro & sizebaro;
    labels.timestamp = labels.timestamp(keep);
    labels.value = labels.value(keep);
    labels.subject = labels.subject(keep);
    labels.text = labels.text(keep);
    labels.acce = labels.acce(keep);
    labels.gyro = labels.gyro(keep);
    labels.baro = labels.baro(keep);
    
    %measure duration of each clip
    for k = 1:length(labels.acce)
        if ~isempty(labels.acce{k})
            dur(k) = max(labels.acce{k}(:,1))-min(labels.acce{k}(:,1));
        end
    end
    
    %delete clips < 7.5 secs
    id = dur>7.5;
    labels.timestamp = labels.timestamp(id);
    labels.value = labels.value(id);
    labels.subject = labels.subject(id);
    labels.text = labels.text(id);
    labels.acce = labels.acce(id);
    labels.gyro = labels.gyro(id);
    labels.baro = labels.baro(id);
    dur = dur(id)';
    
    %jitter
    N = length(labels.timestamp);   %# of clips
    starttimes = 5*rand(N,1);  %uniform
    endtimes = min(starttimes+5,dur);
    
    %find clips less than 5 secs long
    durchop = endtimes-starttimes;
    idshort = durchop < 5;
    while any(idshort)
        n = sum(idshort);
        starttimes(idshort) = 5*rand(n,1);  %normal with std 1 sec
        endtimes(idshort) =  min(starttimes(idshort)+5,dur(idshort));
        durchop = endtimes-starttimes;
        idshort = durchop < 5;
    end
    
    %cutting the clips
    starttimestamp = zeros(length(labels.acce),1);
    for c = 1:length(labels.acce)
        labels.acce{c} = sortrows(labels.acce{c});
        starttimestamp(c) = labels.acce{c}(1,1);
        labels.acce{c}(:,1) = labels.acce{c}(:,1)-starttimestamp(c);
        timediff = abs(labels.acce{c}(:,1)-starttimes(c));
        [~,indmin] = min(timediff);
        timediff2 = abs(labels.acce{c}(:,1)-endtimes(c));
        [~,indmin2] = min(timediff2);
        
        labels.acce{c} = labels.acce{c}(indmin:indmin2,:);
    end
    for c = 1:length(labels.gyro)
        labels.gyro{c} = sortrows(labels.gyro{c});
        labels.gyro{c}(:,1) = labels.gyro{c}(:,1)-starttimestamp(c);
        timediff = abs(labels.gyro{c}(:,1)-starttimes(c));
        [~,indmin] = min(timediff);
        timediff2 = abs(labels.gyro{c}(:,1)-endtimes(c));
        [~,indmin2] = min(timediff2);
        
        labels.gyro{c} = labels.gyro{c}(indmin:indmin2,:);
    end
    for c = 1:length(labels.baro)
        labels.baro{c} = sortrows(labels.baro{c});
        labels.baro{c}(:,1) = labels.baro{c}(:,1)-starttimestamp(c);
        timediff = abs(labels.baro{c}(:,1)-starttimes(c));
        [~,indmin] = min(timediff);
        timediff2 = abs(labels.baro{c}(:,1)-endtimes(c));
        [~,indmin2] = min(timediff2);
    
        labels.baro{c} = labels.baro{c}(indmin:indmin2,:);
    end
    
    figure, histogram(starttimes)
    
    %remove clips with too few datapoints
    sizeacc = cellfun(@(x) size(x,1)>100,labels.acce);
    sizegyro = cellfun(@(x) size(x,1)>100,labels.gyro);
    sizebaro = cellfun(@(x) size(x,1)>10,labels.baro);
    keep = sizeacc & sizegyro & sizebaro;
    labels.timestamp = labels.timestamp(keep);
    labels.value = labels.value(keep);
    labels.subject = labels.subject(keep);
    labels.text = labels.text(keep);
    labels.acce = labels.acce(keep);
    labels.gyro = labels.gyro(keep);
    labels.baro = labels.baro(keep);
    
    filename = 'labels_plus_data_jittered_Unif.mat';
    disp(['saving ', filename])
    save(filename, 'labels')
    
    
    
