%% Jitter Clips
% Create n copies of fall data withj uniformly distributed start time
% from 0 to 5

function data=JitterClips(data,n)

jittersize=3; % range to jitter around (max 5s)
eps=10^-3;
m=length(data.acce);
rng(200)
if n>1
    X=jittersize*rand(n, 1)-.5*jittersize;
else
    X=0;
end
    
acce=cell(m*n,1);
gyro=cell(m*n,1);
baro=cell(m*n,1);

LongInd=cellfun(@(x) length(x)>400,data.acce);
AccMaxT=cell(length(data.acce),1);
temp=cellfun(@(x) x(find((sum(x(:,2:end).^2,2)<max(sum(x(.1*range(x(:,1))<x(:,1)-x(1,1) & x(:,1)-x(1,1)<.7*range(x(:,1)),2:end).^2,2))+eps) ...
    & (sum(x(:,2:end).^2,2)>max(sum(x(.1*range(x(:,1))<x(:,1)-x(1,1) & x(:,1)-x(1,1)<.7*range(x(:,1)),2:end).^2,2))-eps),1),1)-x(1,1),data.acce(LongInd),'UniformOutput',false);
AccMaxT(LongInd)=temp(:);
temp=cellfun(@(x) x(find((sum(x(:,2:end).^2,2)<max(sum(x(:,2:end).^2,2))+eps) ...
    & (sum(x(:,2:end).^2,2)>max(sum(x(:,2:end).^2,2))-eps),1),1)-x(1,1),data.acce(~LongInd),'UniformOutput',false);
AccMaxT(~LongInd)=temp(:);

for i=1:n
    acce((i-1)*m+1:i*m)=cellfun(@(x,T) x(T/1000-2.5-X(i)<(x(:,1)/1000-x(1,1)/1000) ...
        & (x(:,1)/1000-x(1,1)/1000)<T/1000-2.5-X(i)+5,:),...
        data.acce,AccMaxT,'UniformOutput',false);
    gyro((i-1)*m+1:i*m)=cellfun(@(x,T) x(T/1000-2.5-X(i)<(x(:,1)/1000-x(1,1)/1000) ...
        & (x(:,1)/1000-x(1,1)/1000)<T/1000-2.5-X(i)+5,:),...
        data.gyro,AccMaxT,'UniformOutput',false);
    baro((i-1)*m+1:i*m)=cellfun(@(x,T) x(T/1000-2.5-X(i)<(x(:,1)/1000-x(1,1)/1000) ...
        & (x(:,1)/1000-x(1,1)/1000)<T/1000-2.5-X(i)+5,:),...
        data.baro,AccMaxT,'UniformOutput',false);
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