clear 
close all
ParseProbeOutput();
load('FallProbe_TestData.mat');

%load classification tree to evaulate data
load('tree.mat'); %original tree 5 sec windows
load('treeplusdata.mat'); % tree 5 sec windows with data from Luca walking/lying
load('knn_tree'); % knn 5 sec windows 
load('tree_new.mat'); % tree 10 sec windows

% t2 = prune(tree,'level',1); %to reduce overfitting

%obtain longer segments of data
n_s=2; % n_s*5 seconds
l=1;
in=0;
index=0;
t1=0;
x1=0;
y1=0;
z1=0;
s1=datetime(1970,1,1,0,0,0);

time=labels.timestampSTART_END(:,1);
time=datetime(1970,1,1,0,0,[time/1000]);
for i = 1:floor(size(labels.acce,1)/n_s)*n_s
    acce=labels.acce(i,1);
    acce=acce{1,1};
    t{i,1}=acce(:,1);
    x{i,1}=acce(:,2);
    y{i,1}=acce(:,3);
    z{i,1}=acce(:,4);
    
    
    start{i,1}=time(i);
end
for j=1:length(t)/n_s
    i_t=0;
    for h=1:n_s
    t1=[t1;t{h+in,1}+5*i_t];
    x1=[x1;x{h+in,1}];
    y1=[y1;y{h+in,1}];
    z1=[z1;z{h+in,1}];
    %new_start{j,1}=start{1+in,1};
    s1=[s1;start{h+in,1}];
    i_t=h;
    end
    t1(1)=[];
    x1(1)=[];
    y1(1)=[];
    z1(1)=[];
    s1(1)=[];
    new_t{j,1}=t1;
    new_x{j,1}=x1;
    new_y{j,1}=y1;
    new_z{j,1}=z1;
    start_time{j,1}=s1;
    clearvars t1 x1 z1 y1 s1
    t1=0;
    x1=0;
    y1=0;
    z1=0;
    s1=datetime(1970,1,1,0,0,0);
    in=n_s*j;
end
for i=1:length(new_t)
new_acce{i,1}= [new_t{i,1} new_x{i,1} new_y{i,1} new_z{i,1}];
end

%filter to treat data to extract 3-D vectors
fs=50;
fc1=0.7;
nfc1=2*fc1/fs;
[B,A]=butter(2,nfc1,'low');

%initialize structure to save results
activities_prefall=[];
date_fall=[];
lying_time=[];
mov_a_l=[];
lyin_a_m=[];
FallAnalysis = struct('ActBeforeFall',activities_prefall,'TimeImpact',date_fall,'LyingTime',lying_time,'MoveAfterLying',mov_a_l,'LyingAfterMoving',lyin_a_m);

move_in=0; %will be useful to plot data

%to eliminate activity detection data that are not in chronological order
%(happened only in experiments, should not happen with deployment phones)
counter=[];
for i=2:size(activity_detection,2)
if  datetime(1970,1,1,0,0,[str2num(activity_detection(i).Timestamp)]) < datetime(1970,1,1,0,0,[str2num(activity_detection(i-1).Timestamp)])
    counter = i;
end
end

if isempty(counter) == 0
activity_detection(counter:end)=[];
end

%extract probability and time of the impact
P=struct2cell(activity_detection);
t_fall_indices = find(labels.is_fall);
t_fall_ind_1=min(t_fall_indices);
fprintf('\nThe number of recognized fall events (p >= 0.90 and a > 1.6g) is %d', length(t_fall_indices));
[highest_prob_t,t_highest_prob_ind] = max(labels.fall_probab);

if isempty(t_fall_ind_1) == 1
  t_fall_indices(1) = t_highest_prob_ind;
end

%results are obtained for all events with p> 0.90 (if no event reaches the
%threshold, the one with highest p is evaluated)
for p = 1:length(t_fall_indices)

%FIND ACTIVITIES BEFORE FALL
time_impact=labels.timestampSTART_END(t_fall_indices(p),1);  
tupper=datetime(1970,1,1,0,0,[time_impact/1000]); 
tlower=datetime(1970,1,1,0,0,0);
fprintf('\nThe impact happened on %s (Chicago time) \n', tupper - hours(6)); %SELECT hours(5) if light saving is off (before Nov 5), otherwise hours(6)
wholeTime=[];
wholeType=[];
wholeConf=[];
X=[]; %initialize the activities vector

%exctract data from cell
for i=1:size(P,3)
    wholeTime=[wholeTime;str2num(activity_detection(i).Timestamp)];
    wholeType(i).Type=(activity_detection(i).Type);
    wholeConf(i).conf=(activity_detection(i).Conf);
end

%select segment of data before the fall
for j= 1:size(P,3)
    tf(j,:)= isbetween(datetime(1970,1,1,0,0,[str2num(activity_detection(j).Timestamp)]),tlower,tupper);
    tf=logical(tf);
    beforeSegment=wholeTime(tf,1);
end

for k= 1:length(beforeSegment)
    beforeType(k).Type=wholeType(k).Type;
    beforeConf(k).conf=wholeConf(k).conf;
    
end

if isempty(beforeSegment) == 1
    fprintf('\nThere are no activity recognition data available before the impact\n');
else
    
%select activity data segments before the impact
n_segments=10; %20=10 min, 10=5min etc..   n is set to 10, if there are less segments the script will use less segments and report it
if length(beforeSegment)<=10
    n_segments=length(beforeSegment);
end
for h=1:n_segments
    mySegment(h,:)=beforeSegment(end+1-h);
    myActivity(h).Type=beforeType(end+1-h).Type;
    myConf(h).conf=beforeConf(end+1-h).conf;
end

%check for time differences between portions of activity data and get rid
%of those data that happen after a hole greater than one minute
for h=1:n_segments-1
    timeDiff(h,:)= datetime(1970,1,1,0,0,[mySegment(h)]) - datetime(1970,1,1,0,0,[mySegment(h+1)]); 
end    

index_segment=[];
for h=1:n_segments-1
    if  timeDiff(h) > minutes(1.5)
        index_segment(h) = h;
end
end 

if any(index_segment) == 1
time_holes=find(index_segment);
n_segments=min(time_holes);
clear mySegment
end

for h=1:n_segments
    mySegment(h,:)=beforeSegment(end+1-h);
    myActivity(h).Type=beforeType(end+1-h).Type;
    myConf(h).conf=beforeConf(end+1-h).conf;
end
if size(mySegment) == [1,1]
   predact{1,1}=['STILL'];
   fprintf('The probability of the fall was %f \n', highest_prob_t);
   fprintf('There are no available activity recognition data before the fall: we assume that the subject was STILL\n');
else
mySegment = sort(mySegment,'descend');  
m = tupper - datetime(1970,1,1,0,0,[mySegment(end)]); % m = tupper - datetime(1970,1,1,0,0,[beforeSegment(end-n_segments)]);
m1= tupper - datetime(1970,1,1,0,0,[mySegment(1)]);   % m1= tupper - datetime(1970,1,1,0,0,[beforeSegment(end-1)]);     
                                                      
%sum all the confidences of each activity type and select the three highest
%ones with percentages
activities = ['IN_VEHICLE';'ON_BICYCLE';'ON_FOOT   ';'RUNNING   ';'STILL     ';'TILTING   ';'UNKNOWN   ';'WALKING   '];
actname = cellstr(activities);
sum_walking=zeros(1,1);
sum_still=zeros(1,1);
sum_running=zeros(1,1);
sum_onfoot=zeros(1,1);
sum_invehicle=zeros(1,1);
sum_unknown=zeros(1,1);
sum_tilting=zeros(1,1);
sum_onbicycle=zeros(1,1);
for f=1:n_segments
    a=myActivity(f);
    b=myConf(f);
    a1=find(strcmp(a.Type,'WALKING'));
    a2=find(strcmp(a.Type,'STILL'));
    a3=find(strcmp(a.Type,'RUNNING'));
    a4=find(strcmp(a.Type,'ON_FOOT'));
    a5=find(strcmp(a.Type,'IN_VEHICLE'));
    a6=find(strcmp(a.Type,'UNKNOWN'));
    a7=find(strcmp(a.Type,'TILTING'));
    a8=find(strcmp(a.Type,'ON_BICYCLE'));
        if isempty(a1)
        sum_walking=sum_walking;
        else
        sum_walking=sum_walking+b.conf(a1);
        end
        if isempty(a2)
        sum_still=sum_still;
        else
        sum_still=sum_still+b.conf(a2);
        end
        if isempty(a3)
        sum_running=sum_running;
        else  
        sum_running=sum_running+b.conf(a3);
        end  
        if isempty(a4)
        sum_onfoot=sum_onfoot;
        else
        sum_onfoot=sum_onfoot+b.conf(a4);
        end
        if isempty(a5)
        sum_invehicle=sum_invehicle;
        else
        sum_invehicle=sum_invehicle+b.conf(a5);
        end
        if isempty(a6)
        sum_unknown=sum_unknown;
        else
        sum_unknown=sum_unknown+b.conf(a6);
        end
        if isempty(a7)
        sum_tilting=sum_tilting;
        else
        sum_tilting=sum_tilting+b.conf(a7);
        end
        if isempty(a8)
        sum_onbicycle=sum_onbicycle;
        else
        sum_onbicycle=sum_onbicycle+b.conf(a8);
        end
end
X=[sum_invehicle sum_onbicycle sum_onfoot sum_running  sum_still sum_tilting sum_unknown sum_walking ];
max_value=max(X);
ind=find(X == max(X(:)));
predact= cellstr(actname(ind,1));
sum_activities=sum(X);
X1=unique(X);
fprintf('The probability of the fall was %f \n', highest_prob_t);
if length(X1) >=3
second_value=X1([end-1]);
third_value=X1([end-2]);
ind1=find(X == second_value);
ind2=find(X == third_value);
predact_2= cellstr(actname(ind1,1));
predact_3= cellstr(actname(ind2,1));
Normalized_predact=[max_value/sum_activities;second_value/sum_activities;third_value/sum_activities];
if size(beforeType(end).Type,2) == 1 && size(beforeType(end-1).Type,2) == 1
   fprintf('\nThe predominant activity in the first available %.0f minutes segment before the impact was:\n %s ', n_segments*0.5, predact{1,1});
   fprintf('\nThe first segment of available activity data starts %.2f minutes and ends %.2f minutes before the fall event\n',minutes(m),minutes(m1));  
elseif size(predact,1) == 2
fprintf('\nThe predominant activities in the first available %.0f minutes segment before the impact were:\n 1.%s and %s %.2f\n 2.%s %.2f\n 3.%s %.2f\n \n', n_segments*0.5, predact{1,1},predact{2,1}, Normalized_predact(1), predact_2{1,1},Normalized_predact(2), predact_3{1,1},Normalized_predact(3));
fprintf('The first segment of available activity data starts %.2f minutes and ends %.2f minutes  before the fall event\n', minutes(m),minutes(m1));   
else
fprintf('\nThe predominant activities in the first available %.0f minutes segment before the impact were:\n 1.%s %.2f\n 2.%s %.2f\n 3.%s %.2f\n \n', n_segments*0.5, predact{1,1},Normalized_predact(1), predact_2{1,1},Normalized_predact(2), predact_3{1,1},Normalized_predact(3));
fprintf('The first segment of available activity data starts %.2f minutes and ends %.2f minutes  before the fall event\n', minutes(m),minutes(m1));       
end
elseif length(X1) == 2
fprintf('\nThe predominant activities in the first available %.0f minutes segment before the impact was:\n 1.%s %.2f\n', n_segments*0.5, predact{1,1},max_value/sum_activities);
fprintf('The first segment of available activity data starts %.2f minutes and ends %.2f minutes  before the fall event\n', minutes(m),minutes(m1));     
end
end
end

if isempty(X)==0;
figure(1);
bar(X/sum_activities);
ylabel('Probability of each activity');
xlabel('Activities before falling');
set(gca,'XTick',1:16)
names = activities;
names(1,:) = 'IN VEHICLE';
names(2,:) = 'ON BICYCLE';
names(3,:) = 'ON FOOT   ';
set(gca,'XTickLabel',names)
end


%FIND LYING TIME
%select segment of accelerometer data starting one window after the fall
tlower1=datetime(1970,1,1,0,0,[time_impact/1000]); 
startvect=start_time;  
index = find([startvect{:}] == tlower1);
index_new=ceil(index/n_s); %round up to go to the exact cell
a=cell(size(new_acce,1)-index_new,1);
 for i=1:size(new_acce,1)-index_new
     stored_a(i)={new_acce(index_new-1+i,1)}; 
     a(i,1)=stored_a{1,i};
 end
 
%calculate variance of x,y,z and angle between consecutive window vectors
for i=1:size(a,1)
    z=cell2mat(a(i,1));
    rowx=filtfilt(B,A,z(:,2));
    rowy=filtfilt(B,A,z(:,3));
    rowz=filtfilt(B,A,z(:,4));
    x_val=mean(rowx);
    y_val=mean(rowy);
    z_val=mean(rowz);
    vect_cross(:,i)=[x_val; y_val; z_val];
    norma=sqrt((z(:,2).^2)+(z(:,3).^2)+(z(:,4).^2));
    myVariance(i,:)=var(norma);
    myVariancex(i,:)=var(rowx);
    myVariancey(i,:)=var(rowy);
    myVariancez(i,:)=var(rowz);
end
vect_cross(:,1)=[];
Vectcross=vect_cross;
clear vect_cross;
myVariance(1)=[];
myVariance(end)=[];
myVariancex(1)=[];
myVariancex(end)=[];
myVariancey(1)=[];
myVariancey(end)=[];
myVariancez(1)=[];
myVariancez(end)=[];
myVariance1=myVariance;
clear myVariance;
myVariancex1=myVariancex;
clear myVariancex;
myVariancey1=myVariancey;
clear myVariancey;
myVariancez1=myVariancez;
clear myVariancez;
for u = 1:size(Vectcross,2)-1
    angolino(u)=dot(Vectcross(:,u),Vectcross(:,u+1))./(norm(Vectcross(:,u)).*norm(Vectcross(:,u+1)));
    clearvars  norm y z in i  x_val y_val z_val  m  bell a stored_a
end
angolino1=angolino;
clear angolino;
myAngle=acosd(angolino1)';

Var2=myAngle;
Var3=myVariance1;

%use classification tree (trained on experiments data) to classify if a
%segment means lying (0) or movement (1)
myFeatures=table(Var2,Var3);
[outcome,conf_inte]=predict(tree_new,myFeatures); %or t2 using prune, or treeplusdata which is richer, or knn
c_i_lying=conf_inte(:,1);
c_i_moving=conf_inte(:,2);
varxy=find(myVariancex1<0.02 & myVariancey1<0.02); %to avoid false positives we set a threshold for each variance axis
varyz=find(myVariancez1<0.02 & myVariancey1<0.02);
varxz=find(myVariancex1<0.02 & myVariancez1<0.02);
outcome(varxy)=0;
outcome(varxy)=0;
outcome(varxy)=0;

move_index=find(outcome==1, 1,'first');
if move_index ~= 1
tot_ci=(sum(c_i_lying(1:move_index-1))+c_i_moving(move_index))/move_index;
elseif isempty(move_index) == 1
move_index = length(outcome);
tot_ci=sum(c_i_lying(1:move_index))/move_index;
fprintf('\nWarning! No movement was detected for the whole available data set after the fall event, subject could lie for a time equal or greater than the one shown below');
else
tot_ci=c_i_moving(move_index);
end
lying_t=(move_index)*(n_s*5);

%evaluate the time of movement after lying
outcomeaft_l=outcome;
outcomeaft_l(1:move_index)=[];
mov_after_lying=find(outcomeaft_l==0, 1,'first');
% if isempty(mov_after_lying) == 1
%    mov_after_lying = 1;
% end
%if the movement lasts less than 10 sec it could be an outlier
move_index_2=0; %in case mov_after_lying >=2.
if mov_after_lying > 1
outcomeaft_l(1:mov_after_lying-1)=[];    
move_index_2=find(outcomeaft_l==1, 1,'first')-1;
elseif mov_after_lying <= 1
outcomeaft_l(1)=[];
move_index_2=find(outcomeaft_l==1, 1,'first');
end
if isempty(move_index_2) == 1
move_index_2=length(outcomeaft_l);
fprintf('\nWarning! No movement was detected for the whole available data set after the fall event, subject could lie for a time equal or greater than the one shown below');

end
% if mov_after_lying <= 1
% movement=mov_after_lying*(5*n_s);
% else
movement=mov_after_lying*(5*n_s);
lying_2=move_index_2*(5*n_s);

fprintf('\nThe subject did not move for %d or less seconds (CI %.4f) after the fall before getting back up.\nThe movement lasted %d seconds and after the movement the subject stopped moving again for %d seconds \n', lying_t, tot_ci, movement, lying_2 );

% end
%saving data in a structure
FallAnalysis(p).ActBeforeFall=predact{1,1};
FallAnalysis(p).TimeImpact=tupper; 
FallAnalysis(p).LyingTime=lying_t;
FallAnalysis(p).MoveAfterLying=movement;
FallAnalysis(p).LyingAfterMoving=lying_2;
%plotting accelerometer data from four windows before the fall till the
%window after restoring movement after the fall
%Fall is in window 5
if isempty(mov_after_lying) == 1
    mov_after_lying = 0;
end
mov_tot=move_index+mov_after_lying+move_index_2;
if mov_tot >= 50
   mov_tot= 50;  % matlab has trouble plotting more than 50 windows
end
if index_new <= 5
    index_new = 5;
end
for i=1:6+mov_tot
x=cell2mat(new_acce(i+index_new-5,1));
figure(1+i+move_in); %the 1+ is due to the first figure being the histogram
plot(x(:,1),x(:,2));
hold on
plot(x(:,1),x(:,3));
hold on
plot(x(:,1),x(:,4));
title('Accelerations along three axes');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
legend('x-axis','y-axis','z-axis');
end
move_in=i+move_in;
end

save FallAnalysis FallAnalysis
% clear
