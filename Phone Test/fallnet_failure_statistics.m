load FallProbe_TestData_cut

fallnet_success = labels.timestamp(:,2);
fallnet_success(:,2) = 1;
fprintf('fallnet evaluated %i clips\n',  size(fallnet_success,1));

% remove duplicates
fallnet_success=sortrows(fallnet_success);
fallnet_success(diff(fallnet_success(:,1))==0,:)=[];
fprintf('fallnet evaluated %i distinct clips\n',  size(fallnet_success,1));

fallnet_failure = labels.failure.ts(:,2);
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
