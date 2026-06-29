function allDat = saveSpatialDownsampledIm(serverPath,saveFolder,spatialBin)
% This function saves the spatially downsampled resting state data

if exist('spatialBin','var') == 0; spatialBin = 3; end
if exist('saveFolder','file') == 0; [~,~] = mkdir([saveFolder '\Spatial Downsample SS3']); end % Make the directory to save the mat files

datName  = 'Data_RS_10Hz_SS3_';
numFiles = length(dir([saveFolder '\Spatial Downsample SS3' '/*.mat']));

if numFiles == 0
    % Read .blk files and save them into spatially downsampled data
    disp('Saving spatially downsampled data...');
    cd(serverPath);
    allFiles = dir('Data*.BLK');

    % To make sure that the export of the blk files happen in proper order
    for iF = 1:length(allFiles)
        indB = strfind(allFiles(iF).name,'B'); indB = indB(1);
        indDot = strfind(allFiles(iF).name,'.');
        fnBlk(iF,1) = str2double(allFiles(iF).name(indB+1:indDot-1));
    end

    [~, sortInd] = sortrows(fnBlk);
    allFiles     = allFiles(sortInd);

    if length(allFiles)<181; fileLen = length(allFiles); else fileLen = 181; end

    imageFile = OIReadStim([allFiles(1).folder '\' allFiles(1).name],1,'v'); % Get file size
    imageFile = imresize(imageFile,1/spatialBin);
    imSize = size(imageFile); clear imageFile

    % Read all blk files...
    parfor iFile = 1:fileLen
        % Load the .blk file
        imageFile = OIReadStim([allFiles(iFile).folder '\' allFiles(iFile).name],1,'v');

        % Spatially downsample the image by 1/spatialBin
        imageFile = imresize(imageFile,1/spatialBin);
        imageFileR = reshape(imageFile, [size(imageFile,1)*size(imageFile,2) size(imageFile,3)]);

        % Concatenate the downsampled image
        downSampledIm{iFile,1} = imageFileR;
    end  

    % Save spatially downsampled data in .mat files (300 frames in one file)
    numIter = ceil(length(downSampledIm)/6);
    index = 1;

    for iIter = 1:numIter
        clear frames_temp
        if ((index+6)<= length(downSampledIm))
            frames_temp = reshape(cell2mat(downSampledIm(index:index+5)'),[imSize(1) imSize(2) 300]);
            index = index+6;
        else
            frames_temp = reshape(cell2mat(downSampledIm(index:end)'),imSize(1),imSize(2),[] );
        end
        save([saveFolder '\Spatial Downsample SS3\' datName num2str(iIter) '.mat'],'frames_temp');
    end
    allDat = reshape(cell2mat(downSampledIm(:)'),imSize(1), imSize(2),[]);
end
end

