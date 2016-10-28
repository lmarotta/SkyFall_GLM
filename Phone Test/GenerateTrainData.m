%% Add daily non-falls data to the falls data file for training model
%load unlabeled home train data
load TrainData_Home_unlab.mat
labels_home = labels;

%include labels if not there
N = length(labels_home.winsize); %number of clips
labels_home.value = 9*ones(N,1);   %only non-falls (labels >=9)


%fall labeld lab data
load labels_full.mat
inds = cellfun(@isempty,labels.subject); %empty clips
%remove empty clips
labels.timestamp(inds) =[];
labels.value(inds) = [];
labels.subject(inds) = [];
labels.acce(inds) = [];
labels.baro(inds) = [];
labels.gyro(inds) = [];
labels = rmfield(labels,'text');


%concatenate home data
labels.timestamp = [labels.timestamp labels_home.timestampSTART_END(:,1)'];
labels.value = [labels.value labels_home.value'];
labels.acce = [labels.acce labels_home.acce'];
labels.baro = [labels.baro labels_home.baro'];
labels.gyro = [labels.gyro labels_home.gyro'];

%assign subjects to home data
Subj = unique(labels.subject);
Nsubj = length(unique(labels.subject)); %5 amputees, 9 controls
Nhome = length(labels_home.winsize);
newsubj = repmat(Subj,1,round(Nhome/Nsubj));
newsubj = newsubj(1:Nhome);
labels.subject = [labels.subject newsubj];

save labels_full_whome.mat labels