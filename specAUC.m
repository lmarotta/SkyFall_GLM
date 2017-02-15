%function that computes the specificity at a given level of sensitivity
%(sens)
%this fcn will be used as an input to bootci

function S_AUC = specAUC(TL,sens)

%first find the threshold at which 90% sensitivity is achieved
TL = sortrows(TL,-1);
% t = 0.01:0.01:1-0.01;
% for thr = 1:length(t)
%     Sens_t(thr) = sum(TL(TL(:,1)>=t(thr),2))/( sum(TL(TL(:,1)>=t(thr),2)) + sum(TL(TL(:,1)<t(thr),2)) );
% end
% ind90 = find(Sens_t >= sens); ind90 = ind90(end); T90 = t(ind90); %the thrshold which yields 90% Sens

inds=find(TL(:,2),ceil(sum(TL(:,2))*sens));
I=inds(end);
S=1-sum(~TL(1:I,2))/sum(~TL(:,2));

Sens=sum(TL(TL(:,1)>.5,2))/sum(TL(:,2));
Spec=sum(~TL(TL(:,1)<.5,2))/sum(~TL(:,2));

%now computes the specificity 
% S = sum(~TL(TL(:,1)<T90,2))/( sum(~TL(TL(:,1)<T90,2)) + sum(~TL(TL(:,1)>=T90,2)) );

%AUC 
[~,~,~,AUC] = perfcurve(TL(:,2),TL(:,1),1);

S_AUC = [S AUC Sens Spec];
end
