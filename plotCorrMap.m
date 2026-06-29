function [corrMap] = plotCorrMap(seedMat,allDat,plotFlag)
% Plots and gives out the correlation map between the seed and the rest of
% the image. 

if exist('plotFlag','var')==0 || isempty(plotFlag); plotFlag = 0; end
imSize = size(allDat);
allMat = reshape(allDat,imSize(1)*imSize(2),imSize(3));
cc = corr(seedMat,allMat','rows','complete');
corrMap = reshape(cc,[imSize(1) imSize(2)]);
if plotFlag
    figure; imagesc(corrMap); colormap jet; colorbar;
    clim([ 0 1]);
end
end

