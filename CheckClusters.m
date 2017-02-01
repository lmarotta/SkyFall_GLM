fvar.nzstd(fvar.nzstd)=nz_ind;

coeffs=fvar.nzstd;

clust_inds=zeros(length(FeatClusters),1);
halfFeat=sum(cell2mat(FeatClusters(1:end-1,2)));
featCounts=cell2mat(FeatClusters(1:end-1,2));
for i=1:length(FeatClusters)-1
    inds=1+sum(featCounts(1:i-1)):sum(featCounts(1:i-1))+featCounts(i);
    inds=[inds inds+halfFeat];
    clust_inds(i)=sum(coeffs(inds));
    ClustInds{i}=coeffs(inds);
end

clust_inds(i+1)=sum(coeffs(end-6:end));
ClustInds{i}=coeffs(end-6:end);