function [tempHighPass,allDat,nuisanceRemoved,spSmooth,tempBandPass] = getPreProcessedDataRestingState(serverDataPath,dataDir,runName,numFiles,spatialBin,datName)
% Performs downsampling,temporal binning, smoothing, first frame
% subtraction and high pass filtering of data
% tempHighPass - after temporal high pass filtering of data
% allDat - all processing before temporal high pass filtering
% nuisanceRemoved - regresses out the nuisance signals
% spSmooth - spatial smoothing of the image after regression
% tempBandPass - final image after processing

% Check if the downsampled data is already saved
if numFiles == 0
    saveSpatialDownsampledIm(serverDataPath,dataDir,spatialBin);
    numFiles = length(dir([dataDir '\Spatial Downsample SS3' '/*.mat'])); % Reupdate the number of MATLAB files stored
else
    disp('Imaging data has been spatially downsampled...');
end
cd([dataDir '\Spatial Downsample SS3' ]);

% Load the vessel mask
if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
    clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
else
    clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
end
if exist([dataDir '\skullMask.bmp'],'file') == 0
    allCortexMask = imread([dataDir '\skullMask.png']); % Has vessels in this mask
else
    allCortexMask = imread([dataDir '\skullMask.bmp']); % Has vessels in this mask
end

clipMask       = imresize(clipMask,1/3);      % Resize mask
clipMask       = clipMask(:,:,1)>0;           % Converting to 0s and 1s
allCortexMask  = imresize(allCortexMask,1/3); % Resizing cortex mask with vessels
allCortexMask  = allCortexMask(:,:,1)>0;
vesselMask     = ~clipMask & allCortexMask;   % Vessel mask
tempVesselMask = double(vesselMask);

% Get skull mask
if exist([dataDir '\maskSkull0' runName(end) '.bmp'],'file') == 0
    skullMask = imread([dataDir '\maskSkull0' runName(end) '.png']);
else
    skullMask = imread([dataDir '\maskSkull0' runName(end) '.bmp']);
end
skullMask                      = imresize(skullMask,1/3);
skullMask                      = skullMask(:,:,1)>0;
tempSkullMask                  = double(skullMask);

% Append all files if the spatially downsampled data is present
if ~exist('allDat','var')
    allDat = [];
    disp('Appending spatially downsampled data...');
    if numFiles>31; numFiles = 31; end

    parfor iF = 1:numFiles
        tempData{iF,1} = load([datName num2str(iF) '.mat']);
        try tempData{iF,1}.frames_temp;
            tempData{iF,1} = tempData{iF,1}.frames_temp;
        catch
            tempData{iF,1} = tempData{iF,1}.frames_10Hz_temp;
        end
    end

    allDat = single(cat(3,tempData{:})); % Appending data
end

imSize = size(allDat);
allDat = reshape(allDat,[imSize(1)*imSize(2) imSize(3)]);

%%%%%%%%%%%%%%% Step 1: Temporal binning - Downsample from 10 to 2 Hz %%%%%%%%%%%%%%%
binSize      = 5;
disp('Downsampling data from 10Hz to 2Hz...');
for iF = 1:(size(allDat,2)/binSize)
    meanTempBin(:,iF) = mean(allDat(:,((iF-1)*binSize+1):(iF*binSize)),2);
end

%%%%%%%%%%%%%%% Step 2: Temporal smoothing - moving median filter %%%%%%%%%%%%%%%
disp('Temporal smoothing...');
allDat = movmedian(meanTempBin,4,2);

%%%%%%%%%%%%%%% Step 3: First frame subtraction %%%%%%%%%%%%%%%
disp('First frame subtraction...');
firstFrame   = allDat(:,1);
allDat       = ((allDat - firstFrame)./firstFrame)*100; % (delR/R)*100
 
%%%%%%%%%%%%%%% Step 4: Temporal high pass filter (>0.005 Hz) %%%%%%%%%%%%%%%
fs    = 2; % Sampling frequency (after downsampling)
disp('Temporal high pass filter (>0.005 Hz)...');
[b,a] = butter(3,0.005/(fs/2),'high');

parfor iF = 1: size(allDat,1)
    tempHighPass(iF,:) = filtfilt(b,a,double(allDat(iF,:)));
end

tempHighPass     = reshape(tempHighPass,[imSize(1) imSize(2) size(allDat,2)]);

%%%%%%%%%%%%%%% Step 5: Remove nuisance signal along the time course - skull and vessels %%%%%%%%%%%%%%%
% Step 5.1: Find the regressors corresponding to the vessels
disp('Removing nuisance signal along time course...');
regVessel                        = getRegressors(vesselMask,tempHighPass);
tempVesselMask(tempVesselMask>0) = regVessel.idx;

% Step 5.2: Get the regressor for the skull
% Load the skull mask
regSkull                       = getRegressors(skullMask,tempHighPass);
tempSkullMask(tempSkullMask>0) = regSkull.idx;

% Get timecourse of the cortex (not including the vessels)
cortexOnlyMask                 = clipMask & allCortexMask; % No vessels in this mask
timeCourseAll                  = tempHighPass;
timeCourseAll(~cortexOnlyMask) = 0;
timeCourseGlobal               = squeeze(mean(timeCourseAll,[1,2]));

% Regressing out the nuisance signal
X                    = [ones(size(tempHighPass,3),1) regVessel.regressor regSkull.regressor timeCourseGlobal];
tempHighPassReshaped = reshape(tempHighPass,[imSize(1)*imSize(2) size(tempHighPass,3)]);
clear b

parfor iFrame = 1:size(tempHighPassReshaped,1)
    b(:,iFrame)         = regress(tempHighPassReshaped(iFrame,:)',X);
    regFrames(:,iFrame) = tempHighPassReshaped(iFrame,:)' - X*b(:,iFrame);
end

nuisanceRemoved  = reshape(regFrames',[imSize(1) imSize(2) size(regFrames',2)]);

%%%%%%%%%%%%%%% Step 6: Spatial smoothing - Gaussian kernel for 5 pixels (2*ceil(2*sigma)+1 %%%%%%%%%%%%%%%
disp('Spatial smoothing ...');
for iF = 1:size(nuisanceRemoved,3)
   spSmooth(:,:,iF)= imageFilter_LPHP_nsc15(nuisanceRemoved(:,:,iF), 2, 0, []);
end

%%%%%%%%%%%%%%% Step 7: Temporal bandpass filter for 0.01 - 0.1 Hz %%%%%%%%%%%%%%%
clear b a
disp('Temporal band pass filtering from 0.01-0.1 Hz...');
% [b,a]         = butter(3, [0.01/(fs/2) 0.1/(fs/2)],'bandpass');
[z,p,k] = butter(3,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);

spSmoothR     = reshape(spSmooth,[imSize(1)*imSize(2) size(spSmooth,3)]);

parfor iF = 1: size(spSmoothR,1)
    tempBandPassR(iF,:) = filtfilt(sos,g,spSmoothR(iF,:));
end

tempBandPass  = single(reshape(tempBandPassR,[imSize(1) imSize(2) size(spSmooth,3)]));

end

