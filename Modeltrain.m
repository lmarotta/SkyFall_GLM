%Train GLM model
%Input: raw features F; target labels L; Hyperparams: alpha,lambda; 
%no_baro - don't use barometer features if 1
%Output: fvar: structure (normalization coeffs)
%        b: (non zero model coeffs including constant term)
%        nz_ind: non zero coeffs indices 
function [fvar,b,nz_ind]=Modeltrain(F,L,alpha,lambda,no_baro)

epsF = 1e-6; %min std dev for standardizing features
%normalize features
fvar.std= std(F);
fvar.mu= mean(F);
fvar.eps= epsF;
nzstd = (fvar.std./mean(fvar.std) > fvar.eps);
fvar.nzstd= nzstd; %features w nonzero std dev
FNZ = F(:,fvar.nzstd);
FN = zscore(FNZ);

if no_baro %remove barometer features
    pre_baro=sum(nzstd(1:end-nbaro_features));
    fvar.nzstd(end-nbaro_features+1:end)=0;
    FN=FN(:,1:pre_baro);
end

[B,FitInfo] = lassoglm(FN,L,'binomial','Alpha',alpha,'Lambda',lambda); %train the model
sparsity = 100*(length(B)-length(find(B)))/length(B);
b = B(B~=0); %non zero coefficients (B does not include constant term)
b = [FitInfo.Intercept;b];
nz_ind = (B~=0); %non zero indices (skip constant term)

