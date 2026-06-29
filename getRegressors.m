function reg = getRegressors(mask,dataMat,clusterVal,pcaFlag)
% This function gives the regressors that are further used to remove
% nuisance signals in the image. 

if exist('pcaFlag','var') == 0; pcaFlag = 0; end
if exist('clusterVal','var') == 0; clusterVal = 10; end

% Find the pixels in the timecourse that overlaps with the skull/vessel in the skull/vessel map
maskMat = reshape(mask, [size(mask,1)*size(mask,2) 1]);
imageMat = reshape(dataMat,[size(dataMat,1)*size(dataMat,2) size(dataMat,3)]);
imageMat =  imageMat(maskMat,:);
if pcaFlag; [reg.coeff, reg.score, reg.latent, reg.tSq, reg.expVar, reg.mu] = pca(imageMat); end

% kmeans to get the regressors 
reg.idx = kmeans(imageMat,clusterVal,'MaxIter',500);
for iCoeff = 1:clusterVal
    reg.regressor(:,iCoeff) = mean(imageMat((reg.idx == iCoeff),:),1);
end

end

