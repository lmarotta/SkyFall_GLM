%% Jitter Clips
% Create n copies of fall data withj uniformly distributed start time
% from 0 to 5

n=10; % number of resamplings

load falls_data_10sec
m=length(data.acce);

X=5*rand(n, 1);

acce=cell(m*10,1);
gyro=cell(m*10,1);
baro=cell(m*10,1);

for i=1:n
    acce((i-1)*m+1:i*m)=cellfun(@(x) x(X(i)<(x(:,1)/1000-x(1,1)/1000) & (x(:,1)/1000-x(1,1)/1000)<X(i)+5,:),...
        data.acce,'UniformOutput',false);
    gyro((i-1)*m+1:i*m)=cellfun(@(x) x(X(i)<(x(:,1)/1000-x(1,1)/1000) & (x(:,1)/1000-x(1,1)/1000)<X(i)+5,:),...
        data.gyro,'UniformOutput',false);
    baro((i-1)*m+1:i*m)=cellfun(@(x) x(X(i)<(x(:,1)/1000-x(1,1)/1000) & (x(:,1)/1000-x(1,1)/1000)<X(i)+5,:),...
        data.baro,'UniformOutput',false);
end

data.acce=acce;
data.gyro=gyro;
data.baro=baro;

data.type_str=repmat(data.type_str,[n 1]);
data.subject=repmat(data.subject,[n 1]);
data.location=repmat(data.location,[n 1]);

save falls_jittered data
