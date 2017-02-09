%% Jitter Clips
% Create n copies of fall data withj uniformly distributed start time
% from 0 to 5

function data=JitterClips(data,n)

m=length(data.acce);
rng(200)
X=3*rand(n, 1)+0;

acce=cell(m*n,1);
gyro=cell(m*n,1);
baro=cell(m*n,1);

for i=1:n
    acce((i-1)*m+1:i*m)=cellfun(@(x) x(X(i)<(x(:,1)/1000-x(1,1)/1000) & (x(:,1)/1000-x(1,1)/1000)<X(i)+5,:),...
        data.acce,'UniformOutput',false);
    gyro((i-1)*m+1:i*m)=cellfun(@(x) x(X(i)<(x(:,1)/1000-x(1,1)/1000) & (x(:,1)/1000-x(1,1)/1000)<X(i)+5,:),...
        data.gyro,'UniformOutput',false);
    baro((i-1)*m+1:i*m)=cellfun(@(x) x(X(i)<(x(:,1)/1000-x(1,1)/1000) & (x(:,1)/1000-x(1,1)/1000)<X(i)+5,:),...
        data.baro,'UniformOutput',false);
end

acce=cellfun(@(x) [(x(:,1)-min(x(:,1)))/1000 x(:,2:end)],acce,'UniformOutput',false);
gyro=cellfun(@(x) [(x(:,1)-min(x(:,1)))/1000 x(:,2:end)],gyro,'UniformOutput',false);
baro=cellfun(@(x) [(x(:,1)-min(x(:,1)))/1000 x(:,2:end)],baro,'UniformOutput',false);

data.acce=acce;
data.gyro=gyro;
data.baro=baro;

data.value=repmat(data.value,[n 1]);
data.type_str=repmat(data.type_str,[n 1]);
data.subject=repmat(data.subject,[n 1]);
data.location=repmat(data.location,[n 1]);