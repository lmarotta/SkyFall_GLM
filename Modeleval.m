%Evaluate fall model 

%INPUTS: Features F; True labels L, fvar structure; nz_ind (nonzero
%indices); non-zero model coeffs b;
%Optional Input Args: Thres for computing prediction; Display flag to
%output results

%OUTPUTS: pred: model prediction; conf: prediction confidence 
function [pred,conf,confmat] = Modeleval(F,L,fvar,nz_ind,b,varargin)

if length(varargin) < 2
    Thres =0.5; display =1; %default options
else
    Thres = varargin{1}; 
    display = varargin{2};
end

%%normalize input features
FNZ = F(:,fvar.nzstd);
FN = (FNZ - repmat(fvar.mu(fvar.nzstd),[size(FNZ,1),1])) ./ repmat(fvar.std(fvar.nzstd),[size(FNZ,1),1]); %features normalized
FN = FN(:,nz_ind);
conf= glmval(b, FN, 'logit');
pred= ceil(conf-Thres);

%results
isfall = logical(L);
confmat(:,:)=confusionmat(isfall,pred==1,'order',[false true]);

%plot confusion matrix
if display
    activities = {'Non-Fall','Fall'};
    figure; imagesc(confmat./repmat(sum(confmat,2),[1 2]));
    confmat = confmat./repmat(sum(confmat,2),[1 2]);
    [cmin,cmax] = caxis;
    caxis([0,1]) %set colormap to 0 1
    ax = gca;
    ax.XTick = 1:size(activities,2);
    ax.YTick = 1:size(activities,2);
    set(gca,'XTickLabel',activities,'FontSize',14)
    set(gca,'YTickLabel',activities,'FontSize',14)
    ax.XTickLabelRotation = 45;
    axis square
    title('Fall Detection')
    %add text
    thres = 0.6;
    for i = 1:length(activities)
        for j = 1:length(activities)
            if confmat(i,j) < thres
                col = [1 1 1];
            else
                col = [0 0 0];
            end
            text(i-0.2,j,sprintf('%.3f',confmat(j,i)),'FontSize',12,'FontWeight','bold','Color',col);
        end
    end
%     saveas(gcf,'./Figs/FallDetec.fig')
%     saveas(gcf,'./Figs/FallDetec.jpg')
end    

