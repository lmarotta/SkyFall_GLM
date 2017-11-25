clear
load('FallProbe_TestData170802.mat');
n_s=2;
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