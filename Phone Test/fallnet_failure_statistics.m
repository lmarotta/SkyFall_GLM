load FallProbe_TestData

fallnet_success = labels.evalstart(:,1);
fallnet_success(:,2) = 1;
fprintf('fallnet evaluated %i clips\n',  size(fallnet_success,1));

% remove duplicates
fallnet_success=sortrows(fallnet_success);
fallnet_success(diff(fallnet_success(:,1))==0,:)=[];
fprintf('fallnet evaluated %i distinct clips\n',  size(fallnet_success,1));

fallnet_failure = labels.failure.evalstart(:,1);
fallnet_failure(:,2) = 0;
fprintf('fallnet failed to evaluate %i clips\n',  size(fallnet_failure,1));

% remove duplicates
fallnet_failure=sortrows(fallnet_failure);
fallnet_failure(diff(fallnet_failure(:,1))==0,:)=[];
fprintf('fallnet failed to evaluate %i distinct clips\n',  size(fallnet_failure,1));

fallnet = [fallnet_success; fallnet_failure];

% remove duplicates
fallnet=sortrows(fallnet);
fallnet(diff(fallnet(:,1))==0,:)=[];

duration = fallnet(end,1) - fallnet(1,1);
fprintf('the data was collected during %i minutes\n\n', floor(duration/1000/60));

t = fallnet(:,1)-repmat(fallnet(1,1),size(fallnet(:,1)));
figure, plot(t/1000,fallnet(:,2),'o-'); ylim([-0.1 1.1])

%%  statistics

%histogram of duration of clips
td = (labels.timestampSTART_END(:,2)-labels.timestampSTART_END(:,1));
figure, histogram(td), xlabel('Clip Duration [s]'), ylabel('Frequency of clips')
figure, hold on, subplot(311), histogram(labels.sensor_counts(:,1)), title('acc'), xlabel('# of samples in clip')
subplot(312), histogram(labels.sensor_counts(:,2)), title('gyr')
subplot(313), histogram(labels.sensor_counts(:,3)), title('bar')


% histogram of time for successful prediction
figure, hold on, subplot(411), histogram(labels.duration(:,1),10), title('Success Preparation'), xlabel('duration')
subplot(412), histogram(labels.duration(:,2),10), title('Success Verification')
subplot(413), histogram(labels.duration(:,3),10), title('Success Prediction')
subplot(414), histogram(labels.failure.duration(:,1),10), title('Failures Preparation')

%histogram of fallnet failures
reasons = labels.failure.reason;
unique_reasons =unique(reasons,'stable');
%reasons_count=cellfun(@(x) sum(ismember(reasons,x)),unique_reasons,'un',0);

B = categorical(reasons,unique_reasons);
figure, histogram(B,'BarWidth',0.5)


%% plot start and end timestamps for failure and success 

figure
plot(labels.timestampSTART_END(:,1)-labels.timestampSTART_END(1,1),zeros(length(labels.timestampSTART_END(:,1)),1),'ob'), hold on
plot(labels.timestampSTART_END(:,2)-labels.timestampSTART_END(1,1),zeros(length(labels.timestampSTART_END(:,2)),1),'xb')

plot(labels.failure.timestampSTART_END(labels.failure.timestampSTART_END(:,1)<10^100,1)-labels.failure.timestampSTART_END(1,1),ones(sum(labels.failure.timestampSTART_END(:,1)<10^100),1),'or'), hold on
plot(labels.failure.timestampSTART_END(labels.failure.timestampSTART_END(:,1)<10^100,2)-labels.failure.timestampSTART_END(1,1),ones(sum(labels.failure.timestampSTART_END(:,1)<10^100),1),'xr')
ylim([-0.5 1.5])

figure, subplot(2,1,1), hold on, title('success')
histogram(labels.timestampSTART_END(:,2)-labels.timestampSTART_END(:,1))
subplot(212), hold on, title('Failed')
histogram(labels.failure.timestampSTART_END(:,2)-labels.failure.timestampSTART_END(:,1))

%% plots fallnet failures values
%{
for i=1:size(unique_reasons)
    reason = unique_reasons{i,1};
    failures_ind = strcmp(reasons,reason);
    failures_count = labels.failure.value(failures_ind);
    figure,  histogram(failures_count), title(reason)
end
%}
