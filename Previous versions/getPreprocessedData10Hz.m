function tempBandPass10Hz = getPreprocessedData10Hz(serverDataPath,dataDir,runName,numFiles,spatialBin,datName,tempBinFlag)
% Memory saver version of getPreprocessedDataResting state

if ~exist("tempBinFlag",'var'); tempBinFlag = 1; end
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

clipMask                         = imresize(clipMask,1/3); % Resize mask
clipMask                         = clipMask(:,:,1)>0; % Converting to 0s and 1s
allCortexMask                    = imresize(allCortexMask,1/3); % Resizing cortex mask with vessels
allCortexMask                    = allCortexMask(:,:,1)>0;
vesselMask                       = ~clipMask & allCortexMask; % Vessel mask
tempVesselMask                   = double(vesselMask);

% Get skull mask
if exist([dataDir '\maskSkull0' runName(end) '.bmp'],'file') == 0
    skullMask = imread([dataDir '\maskSkull0' runName(end) '.png']);
else
    skullMask = imread([dataDir '\maskSkull0' runName(end) '.bmp']);
end
skullMask                      = imresize(skullMask,1/3);
skullMask                      = skullMask(:,:,1)>0;
tempSkullMask                  = double(skullMask);

% Check if the downsampled data is already saved
if numFiles == 0
    saveSpatialDownsampledIm(serverDataPath,dataDir,spatialBin);
    numFiles = length(dir([dataDir '\Spatial Downsample SS3' '/*.mat'])); % Reupdate the number of MATLAB files stored
else
    disp('Imaging data has been spatially downsampled...');
end
%%
fds = fileDatastore([dataDir '\Spatial Downsample SS3\*.mat' ],'ReadMode','file','ReadFcn',@load);
dataArray = cell(numel(fds.Files),1);
iF = 1; reset(fds); 

while hasdata(fds)
    dataArray{iF} = read(fds);
    dataArray{iF} = dataArray{iF}.frames_temp;
    iF = iF+1; 
end 
allDat = cat(3,dataArray{:}); clear tempData;
imSize = size(allDat);
allDat = single(reshape(allDat,[imSize(1)*imSize(2) imSize(3)]));



% allDat = tall(single(allDat));


    %%%%%%%%%%%%%%% Step 1: Temporal binning - Downsample from 10 to 2 Hz %%%%%%%%%%%%%%%
    binSize      = 5;
    % s = imSize(3) - mod(imSize(3),binSize); % If the values aren't exactly multiples
    if tempBinFlag
        disp('Downsampling data from 10Hz to 2Hz...');
        for iF = 1:(size(allDat,2)/binSize)
            meanTempBin(:,iF) = mean(allDat(:,((iF-1)*binSize+1):(iF*binSize)),2);
        end
    else
        meanTempBin = allDat;
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

    parfor iF = 1: gather(size(allDat,1))
        tempHighPass(iF,:) = filtfilt(b,a,double(gather(allDat(iF,:))));
    end
    % avgDatHP = filtfilt(b,a,double(avgDat));

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
    % spSmooth       = imgaussfilt(nuisanceRemoved,1);
    disp('Spatial smoothing ...');
    for iF = 1:size(nuisanceRemoved,3)
        spSmooth(:,:,iF)= imageFilter_LPHP_nsc15(nuisanceRemoved(:,:,iF), 2, 0, []);
    end


    %%%%%%%%%%%%%%% Step 7: Temporal bandpass filter for 0.01 - 0.1 Hz %%%%%%%%%%%%%%%
    clear b a
    disp('Temporal band pass filtering from 0.01-0.1 Hz...');
    [b,a]         = butter(3, [0.01/(fs/2) 0.1/(fs/2)],'bandpass');
    spSmoothR     = reshape(spSmooth,[imSize(1)*imSize(2) size(spSmooth,3)]);

    parfor iF = 1: size(spSmoothR,1)
        tempBandPassR(iF,:) = filtfilt(b,a,spSmoothR(iF,:));
    end
    tempBandPass  = reshape(tempBandPassR,[imSize(1) imSize(2) size(spSmooth,3)]);

    % Save the processed image file in the folder
    tempBandPass = single(tempBandPass);

    tempBandPass10Hz = [tempBandPass10Hz ; tempBandPass];


end