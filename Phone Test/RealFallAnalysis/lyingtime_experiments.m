clear all
close all
load('FallProbe_TestData170710.mat');
load('fall_170710.mat');
fs=50;
fc1=0.7;
nfc1=2*fc1/fs;
[B,A]=butter(2,nfc1,'low');

n_s=2; % n_s*5 seconds
l=1;
in=0;
index=0;
t1=0;
x1=0;
y1=0;
z1=0;
s1=datetime(1970,1,1,0,0,0);
time=labels.timestampSTART_END(:,2);
% alloallo=datetime(1970,1,0,0,[time/1000]);
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

activities_prefall=[];
date_fall=[];
lying_time=[];

startvect=start_time;  
cut_five= find(lietime ~= 5 & lietime ~= 15 &  lietime ~= 25 & lietime ~= 35);  
Start=Start(cut_five);
lietime=lietime(cut_five);
mytotalarrow=[];
for p = 1:length(Start)
 
time_impact=Start(p);  

%find lying time

tlower1=datetime(1970,1,1,0,0,[time_impact/1000]); 
tupper1=datetime(1970,1,1,0,0,[1.999999900000000e+12/1000]);

date_diff = time - tlower1;
date_diff(date_diff < 0) = Inf;
[~,indexfall] = min(date_diff);
closest_time = time(indexfall,:);
index = find([startvect{:}] == closest_time);
index_new=ceil(index/n_s);
arrow=zeros((lietime(p)/10),1);
myarrow=[arrow;1];
mytotalarrow=[mytotalarrow;myarrow];
sub_fall=find(mytotalarrow==1);
sub_fall=sub_fall+1;
mytotalarrow(sub_fall)=2;
mytotalarrow(1)=2;
mytotalarrow(end)=[];
 
 for i=1:lietime(p)/10+1
     stored_a(i)={new_acce(index_new+i,1)}; 
     a(i,1)=stored_a{1,i};

 end
b={a};
bell=b{1,1};
for h=1:size(bell,1)
    seg=bell{h,1};
    for u=1:3  
    segm=seg(:,u+1);
    segm=segm(1:248);
    half= floor(length(segm)/2);
    first_half(:,u)=segm(1:half);
    second_half(:,u)=segm(half+1:end);
   
    end
    first_half=num2cell(first_half);
    half_1{h}=first_half;
    second_half=num2cell(second_half);
    half_2{h}=second_half;  
   clearvars first_half second_half
end
 %calculate variance of x,y,z and norm
for in=1:size(bell,1)
    z=bell{in,1};
    h1=cell2mat(half_1{1,in});
    h2=cell2mat(half_2{1,in});
    rowx1=filtfilt(B,A,h1(:,1));
    rowy1=filtfilt(B,A,h1(:,2));
    rowz1=filtfilt(B,A,h1(:,3));
    rowx2=filtfilt(B,A,h2(:,1));
    rowy2=filtfilt(B,A,h2(:,2));
    rowz2=filtfilt(B,A,h2(:,3));
    x_val1=mean(rowx1);
    y_val1=mean(rowy1);
    z_val1=mean(rowz1);
    x_val2=mean(rowx2);
    y_val2=mean(rowy2);
    z_val2=mean(rowz2);
    vect_cross1=[x_val1 y_val1 z_val1];
    vect_cross2=[x_val2 y_val2 z_val2];
    angolino(in,:)=dot(vect_cross1,vect_cross2)./(norm( vect_cross1).*norm(vect_cross2));
    angolo=acosd(angolino);
    norma=sqrt((z(:,2).^2)+(z(:,3).^2)+(z(:,4).^2));
    myVariance(in,:)=var(norma);
end
myVariance=num2cell(myVariance);
myVariance1{p}=myVariance;
angolo=num2cell(angolo);
angolone{p}=angolo;
clearvars vect_cross1 vect_cross2 myVariance norm y z in i  a stored_a x_val1 y_val1 z_val1 x_val2 y_val2 z_val2 testing  angolino m a stored_a bell h1 h2
end


myangle=[];
myvarianza=[];
for i = 1: size(angolone,2)
    myangle=[myangle;angolone{1,i}];
    myvarianza=[myvarianza;myVariance1{1,i}]; 
end
myvarianza=cell2mat(myvarianza);  
myangle=cell2mat(myangle);
table0710=[mytotalarrow,myangle,myvarianza];        
save('table0710.mat','table0710')


