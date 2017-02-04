
F = [];
lf = [];
filespath = 'Z:/Amputee Phones-R01/Home Data Collection/Amputees/';f = dir('./TrainingData/*.mat');
S=[];
for i=1:length(f)
    tempData=load(['./TrainingData/' f(i).name]);
    S=[S; unique(tempData.data.subject)];
end
S=unique(S);
f = dir([filespath,'/Raw/*.mat']);
for i = 1:length(f)
    sprintf('Filename %s',f(i).name)
    l = load([filespath '/Raw/' f(i).name]); labels = l.labels;
    subj=f(i).name(1:5);
    subjid=find(strcmp(subj,S));
    %extract features
    Fnew = HomeDataSetup(labels,1.5); %1.5g threshold for acceleration clips
    Fnew = [ones(size(Fnew,1),1)*[subjid 1 1 9] Fnew]; 
    F = [F;Fnew];
    lf = [lf;size(Fnew,1)*5/60/60]; %store duration of each day
    sprintf('Data length = %.2f h',size(F,1)*5/60/60) 
    disp(lf)
end