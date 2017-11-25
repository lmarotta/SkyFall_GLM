
close all
t = [];
for i = 1:length(activity_detection)
t = [t;str2num(activity_detection(i).Timestamp)];
end
figure(1);
plot(datetime(1970,1,1,0,0,t)-hours(6),'ob');
figure(2);
plot(datetime(1970,1,1,0,0,labels.timestampSTART_END(:,1)/1000)-hours(6),'or');