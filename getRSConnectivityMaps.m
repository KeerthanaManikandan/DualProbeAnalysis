function [rsConnMatrix,avgCorrMapFinalR,seedSig, processedDat,seedLocR] = getRSConnectivityMaps(refSeed,monkeyName)
% This function generates FC map averaged across multiple RS runs for one
% monkey. 

% Get the RS connectivity matrix before  thresholding
% Initializing all the relevant variables
commonDir = 'C:\Users\kem294\Documents\Data';
spatialBin   = 3;
hemisphere   = 'Left';
% saveFlag     = 0; 
% checkBadRuns = 1; 
% monkeyName = 'Bordeaux';
[allDates,allRuns,~,refDir,refImageName, serverPath,greenMapRef] = getMonkeyParamsRS(monkeyName,commonDir,hemisphere);

% Get the necessary variables needed to compute the connectivity maps
clear seedSig
greenMapSize = size(greenMapRef);
greenMap     = double(greenMapRef);
greenMap_RGB = ind2rgb(greenMap,gray(256)); % Convert from grayscale to RGB

% Get the master green mask 
if exist([refDir refImageName '_mask.bmp'],'file') == 0 % Make sure to get what the reference run is outside the function
    greenMapMask = imread([refDir refImageName '_mask.png']);
else
    greenMapMask  = imread([refDir refImageName '_mask.bmp']);
end
greenMapMask = double(greenMapMask(:,:,1) == 255);

% Get the reference seed
refSeed = fliplr(round(refSeed));

% Initializing required variables
numCols = max(cellfun(@(x)size(x,1),allRuns));
tformWithinRuns = cell(size(allDates,1),numCols); 

% Getting maps for all runs individually
for iDate = 1: size(allDates,1)
    for iRun = 1:size(allRuns{iDate,1},1)

        % Get run details
        clear runName expDate serverDataPath dataDir numFiles templateDir datName tempBandPass
        runName       = allRuns{iDate,1}(iRun,:);
        expDate       = allDates(iDate,:);
        if strcmp(expDate,'06_22_2021')
            serverDataPath    = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\302-19_CharlieSheen_SqM\Left Hemisphere\06_22_2021 Craniotomy\' runName];
        else
            serverDataPath = [serverPath expDate '\' runName];
        end

        dataDir  = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        numFiles = length(dir([dataDir '\Spatial Downsample SS3' '/*.mat'])); % loading the downsampled data only
        datName  = 'Data_RS_10Hz_SS3_';
        templateDir = [commonDir '\' monkeyName '_SqM\Left Hemisphere\'];

                
        % Get the processed data or compute the processed data if not saved
        if ~exist([dataDir '\processedFrames.mat'],'file')
            [~,~,~,~,tempBandPass] = getPreProcessedDataRestingState(serverDataPath,dataDir,runName,numFiles,spatialBin,datName);
            save([dataDir '\processedFrames.mat'],'tempBandPass');
            processedDat = tempBandPass; 
        else
            clear processedDatTemp processedDat
            processedDatTemp = load([dataDir '\processedFrames.mat']);
            processedDatTemp = processedDatTemp.tempBandPass;
            processedDat     = processedDatTemp;
        end

        % Load the greens
        if exist([templateDir  expDate '\' allRuns{iDate,1}(iRun,:) '\green0' allRuns{iDate,1}(iRun,5) '_Edited.bmp'],'file') == 0
            greenTemp = imread([dataDir '\green0' allRuns{iDate,1}(iRun,5) '_Edited.png']);
        else
            greenTemp = imread([dataDir '\green0' allRuns{iDate,1}(iRun,5) '_Edited.bmp']);
        end
        greenIm{iDate,iRun}  = greenTemp(:,:,1);

        % Load the clipmask 
         if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
            clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
        else
            clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
         end
         clipMask = imresize(clipMask,1/3); % Resize mask
         clipMask = clipMask(:,:,1)>0; % Converting to 0s and 1s
        
        % Transformation of the seed point to the run-wise maps so as to
        % get individual run-wise correlation maps
        % Linear transformation
        clear tformRunWise
        if exist([dataDir '\tform_RS.mat'],'file') == 0
            [optimizer,metric] = imregconfig('multimodal');
            tformRunWise = imregtform(imresize(greenIm{iDate,iRun},1/spatialBin),imresize(greenIm{iDate,1},1/spatialBin),'affine',optimizer,metric);
            tformWithinRuns{iDate,iRun} = tformRunWise;
            save([dataDir '\tform_RS.mat'],'tformRunWise');
        else
            
            if isempty(tformWithinRuns{iDate,iRun}) 
                load([dataDir '\tform_RS.mat'],'tformRunWise');
                tformWithinRuns{iDate,iRun} = tformRunWise;
            end
        end
        
        % Non-linear transformation 
         if iRun == 1
            if (exist([templateDir expDate '\imageRegisPoints_RS.mat'],'file') == 0) 
                selectPoints = 1; addPoints = 0; checkRegisImage = 1;

                % Get the fixed and moving points and store them in a matfile
                clear temp
                if selectPoints
                    clear movePointsTemp fixedPointsTemp
                    [movPointsTemp,fixedPointsTemp] = cpselect(greenIm{iDate,1},greenMapRef,'Wait',true);
                    movPoints{iDate,1} = movPointsTemp;
                    fixedPoints{iDate,1} = fixedPointsTemp;
                end
                
                if addPoints
                    [movPointsTemp,fixedPointsTemp] = cpselect(greenIm{iDate,1},greenMapRef,movPoints{iDate,1},fixedPoints{iDate,1},'Wait',true);
                    movPoints{iDate,1} = movPointsTemp;
                    fixedPoints{iDate,1} = fixedPointsTemp;
                end

                % In case you want to save again especially when you want to add points or to debug
                if ~(selectPoints || addPoints) 
                    try
                        load([templateDir expDate '\imageRegisPoints_RS.mat']);
                        movPoints{iDate,1} = movPointsTemp;
                        fixedPoints{iDate,1} = fixedPointsTemp;
                    catch
                        disp('Image registration not done, beginning image registration now..');
                        clear movePointsTemp fixedPointsTemp
                        [movPointsTemp,fixedPointsTemp] = cpselect(greenIm{iDate,1},greenMapRef,'Wait',true);
                        movPoints{iDate,1} = movPointsTemp;
                        fixedPoints{iDate,1} = fixedPointsTemp;
                    end
                end

                clear temp OTemp spacingTemp tTemp
                fixedPointsTemp = fixedPoints{iDate,1}; movPointsTemp = movPoints{iDate,1};
                temp = zeros(size(greenMapRef));
                temp(1:size(greenIm{iDate,1},1),1:size(greenIm{iDate,1},2)) = greenIm{iDate,1};
                
                % Perform image registration
                [OTemp,spacingTemp,~] = point_registration(size(greenMapRef),fliplr(fixedPoints{iDate,1}),fliplr(movPoints{iDate,1}));
                [imRTemp{iDate,1},~] = bspline_transform(OTemp,temp,spacingTemp,3); % bicubic spline interpolation
                O{iDate,1} = OTemp;
                spacing{iDate,1} = spacingTemp;

                save([templateDir expDate '\imageRegisPoints_RS.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
                
                if checkRegisImage % Check the green images
                    figure; imshowpair(greenMapRef,imRTemp{iDate,1},'blend'); %colormap gray;
                    title(['Green map from Date: ' allDates(iDate,:)]);
                end
            else % Load the registered image 
                load([templateDir expDate '\imageRegisPoints_RS.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
                movPoints{iDate,1} = movPointsTemp;
                fixedPoints{iDate,1} = fixedPointsTemp;
                O{iDate,1} = OTemp;
                spacing{iDate,1} = spacingTemp;
            end
         end

         % Get connectivity maps
         % Get run-wise correlation maps
         Ox = O{iDate,1}(:,:,1); Oy = O{iDate,1}(:,:,2); dx = spacing{iDate,1}(1); dy = spacing{iDate,1}(2);
         seedLocNonLinear(iDate,:) = bspline_trans_points_double(O{iDate,1},spacing{iDate,1},[refSeed(1) refSeed(2)]); % Non-linear transformation
         [seedLocR(2),seedLocR(1)] = transformPointsInverse(tformWithinRuns{iDate,iRun},seedLocNonLinear(iDate,2),seedLocNonLinear(iDate,1)); % Linear transformation

         if sum(seedLocR<=0)>0 || sum(squeeze(round(seedLocR./spatialBin))<=1)>=1; noSeedSig = false; continue; end % if seed value is negative, incrementing the pixel values
       
         try
             seedSig= calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin), clipMask,fliplr(squeeze(round(seedLocR./spatialBin))),3,processedDat); % Get Gaussian weighted seed signal
             noSeedSig = false;
         catch
             noSeedSig = true;
             continue;
         end

         corrMapTemp(:,:,iDate,iRun) = plotCorrMap(seedSig,processedDat,0); % Get FC map

         % Register the correlation map back to the master green 
         clear corrMapLinear corrMapNonLinear
         scaleRatio = round(greenMapSize(1)/size(processedDat,1),2);
         corrMapLinear = imresize(imwarp(corrMapTemp(:,:,iDate,iRun),tformWithinRuns{iDate,iRun},'OutputView',imref2d(size(corrMapTemp(:,:,iDate,iRun)))),spatialBin,'OutputSize',[greenMapSize(1) greenMapSize(2)]); % Linear transformation
         corrMapNonLinear = bspline_transform_2d_double(Ox,Oy,corrMapLinear,dx,dy,3); % Nonlinear transformation

         % Mutliply the correlation map with master green mask so as to
         % remove the corners outside the cortex 
         corrMapNonLinear = corrMapNonLinear.*greenMapMask;
         corrMap(:,:,iDate,iRun) = imresize(corrMapNonLinear,1/scaleRatio,'OutputSize',[size(processedDat,1) size(processedDat,2)]); % Resize the image
    end
    if noSeedSig; continue; end
end

% Reshape the correlation map
corrMapR = reshape(corrMap,[size(corrMap,1) size(corrMap,2) size(corrMap,3)*size(corrMap,4)]);

% Check if any entity in the matrix is empty
checkMat = squeeze(sum(corrMapR,[1 2])== 0);
corrMapR(:,:,checkMat) = [];

% Average out the correlation maps
avgCorrMap = mean(corrMapR,3,'omitnan');

% Compare the run-wise correlation maps with the average correlation map
% Get the pearson correlation between them
imSize = size(corrMapR);
corrMapMat    = reshape(corrMapR,[imSize(1)*imSize(2) imSize(3)]);
corrMapMatAvg = reshape(avgCorrMap,[imSize(1)*imSize(2) 1]);

for iL = 1:size(corrMapMat,2)
    rho(iL,1) = corr(corrMapMatAvg,corrMapMat(:,iL)); % Get the correlation coefficient
end

% Compute the mean and std dev and reject the runs whose values are greater
% or lesser than 1.5 std dev
rhoM   = mean(rho,'omitnan');
rhoStd = std(rho,'omitnan');
lowThresh = rhoM - (1.5*rhoStd); % rhoM - (1.5*rhoStd);
highThresh = rhoM + (1.5*rhoStd); % rhoM + (1.5*rhoStd);
goodRuns = ((rho >= lowThresh) & (rho <= highThresh));
corrMapMat(:,~goodRuns) = []; % Remove bad runs 

% Get the average of good runs 
corrMapGood = reshape(corrMapMat,[imSize(1) imSize(2) size(corrMapMat,2)]);
rsConnMatrix = mean(corrMapGood,3,'omitnan'); avgCorrMapFinal = rsConnMatrix;
rsConnMatrix = imresize(rsConnMatrix,scaleRatio,'OutputSize',[greenMapSize(1) greenMapSize(2)]);

% Only include pixels which are significantly greater than 0
[~,p]     = ttest(corrMapMat');
reshapedP = reshape(p,[imSize(1) imSize(2)]);
sigPixels = (reshapedP<0.001) & (avgCorrMapFinal>0);
avgCorrMapFinal(~sigPixels) = 0; 

avgCorrMapFinalR = imresize(avgCorrMapFinal,scaleRatio,'OutputSize',[greenMapSize(1) greenMapSize(2)]);
end

function seedSig = calculateSeedSignal(green, clipMask, seed, seedRad,frames)
    % choose seed location (gaussian average of a radius around chosen point)
    % seedRad = 6; % 6 = 250um radius, 12 = 500um radius, 24 = 1000um radius
    clear gaussianDisk green_Seed clipMask_Seed 
    if exist('seedRad','var') == 0; seedRad = 6; end
    seedSig = zeros(size(frames,3),1);
    
    gaussianDisk = fspecial('gaussian',2*seedRad+1,seedRad);
    green_seed = green(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad,:);
    clipMask_seed = clipMask(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad);
    gaussianDisk = gaussianDisk .* double(clipMask_seed);
    gaussianDisk = gaussianDisk ./ sum(gaussianDisk(:));

    for x = 1:2*seedRad+1     % rows are observations, columns are variables
        for y = 1:2*seedRad+1
            seedSig = seedSig + gaussianDisk(x,y) * squeeze(frames(seed(2)+x-seedRad-1,seed(1)+y-seedRad-1,:));

        end
    end
end
