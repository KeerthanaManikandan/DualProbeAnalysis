% Code to check if the FC maps have "effector" and "inter-effector" regions
% Based on Gordon et al 2023 article in Nature 
monkeyName = 'CharlieSheen';

[allDates,datFileNumAll,serverPath,refDate,refDir,refImageName,datFileNameAll,chInCortexProbeA, chInCortexProbeB ,...
    ~,~,probeLabelA,probeLabelB,anesthesiaLevels,heartRate,patchInfo,pairClass,corticalAreaProbeA,...
    corticalAreaProbeB] = getMonkeyParamsDualProbeEphys(monkeyName,commonDir);
hemisphere ='Left';

if exist([refDir refImageName '.bmp'],'file') == 0 % Make sure to get what the reference run is outside the function
    greenMapRef = imread([refDir refImageName '.png']);
else
    greenMapRef = imread([refDir refImageName '.bmp']);
end
greenMapRef  = greenMapRef(:,:,1);

figure; imagesc(greenMapRef);axis image off; colormap gray; hold on;
title('Select the reference site');
refSeed(1,:,:) = ginput(1);

yVals = 0:50:size(greenMapRef,1);
% refSeedT = zeros(20,1,2);
% refSeed(1:20,1,:)= refSeed(1,1,:);
% 
% refSeed(2:20,:,2) = refSeed(1,:,2)-yVals(2:20);
%%
tic;
fileNum = 1; 
for iMap = 1:20
    refSeedT = refSeed;
    refSeedT(:,:,2) = refSeed(:,:,2)-yVals(iMap);
[rsConnMatrix(:,:,iMap),corrMapFinal(:,:,iMap)] = ...
    getRSConnectivityMaps(squeeze(refSeedT(fileNum,:,:))',monkeyName);
end
toc;
save('D:\Data\CharlieSheen_SqM\Left Hemisphere\corrMapsMotorCaudal.mat','rsConnMatrix','corrMapFinal');

%%
figure;
v = VideoWriter('D:\Data\CharlieSheen_SqM\Left Hemisphere\motorMapsCaudalVideo.avi','Uncompressed AVI');  
v.FrameRate = 1;  
open(v);
for iMap = 1:20
    imagesc(ind2rgb(greenMapRef,gray(256))); hold on;
    imagesc(corrMapFinal(:,:,iMap),'alphaData',corrMapFinal(:,:,iMap).*0.9); axis image off;
    colormap jet; clim([0 1]);pause(1);
    frame = getframe(gcf);
    writeVideo(v, frame);
end
close(v);