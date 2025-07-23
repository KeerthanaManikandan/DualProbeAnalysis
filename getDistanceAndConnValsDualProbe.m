function [distSitesAll,connValsAll,refSiteAll, movSiteAll,greenMapRef] = getDistanceAndConnValsDualProbe(monkeyName,allDates,datFileNumAll,refDir,refImageName)
% This function gets the distance between two sites (in mm) and also gets
% the connectivity value of the moving electrode with respect to the
% reference. 

hemisphere ='Left';
if exist([refDir refImageName '.bmp'],'file') == 0 % Make sure to get what the reference run is outside the function
    greenMapRef = imread([refDir refImageName '.png']);
else
    greenMapRef = imread([refDir refImageName '.bmp']);
end
greenMapRef  = greenMapRef(:,:,1);

% Get the distance, connectivity and save it
for iDate = 1:size(allDates,1)
    clc; clear expdate datFileNum distSites connVal refSeed movSeed
    expDate    = allDates(iDate,:);
    datFileNum = datFileNumAll{iDate,1};
    
    if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\distance_conn.mat'],'file')
        for iFile = 1:length(datFileNum)
            fileNum = datFileNum(iFile);
            disp([monkeyName ' ' expDate ' datafile: ' num2str(fileNum)]); 
            if iFile>1
                notSameRef = input('Has the reference seed changed? 1- Yes, 0 - No ');
            else
                notSameRef = 1;
            end
            if notSameRef % New reference site 
                plotCorrMapFlag = 1;
                figure; imagesc(greenMapRef);axis image off; colormap gray; hold on;
                title('Select the reference site');
                refSeed(fileNum,:,:) = ginput(1);
                title('Select the moving site');
                movSeed(fileNum,:,:) = ginput(1);                
                close gcf;
                % Get the average FC map 
                [rsConnMatrix,corrMapFinal] = getRSConnectivityMaps(squeeze(refSeed(fileNum,:,:))',monkeyName);  

                if plotCorrMapFlag % Save the connectivity map for the new site... 
                    figure; imagesc(ind2rgb(greenMapRef,gray(256))); hold on; 
                    imagesc(corrMapFinal,'alphaData',corrMapFinal.*0.9); axis image off; 
                    colormap jet; clim([0 1]);
                    f = gcf; 
                    exportgraphics(f,['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Results\FC_map_' num2str(iFile) '.png'],'Resolution',300);
                    close gcf;
                end 

            else
                figure; imagesc(greenMapRef);axis image off; colormap gray; hold on;
                disp('Select the moving site');
                movSeed(fileNum,:,:) = ginput(1); close gcf;
                refSeed(fileNum,:,:) = refSeed(fileNum-1,:,:);
            end

            % Get the distance (in mm) between the moving seed and the reference seed
            distSites(fileNum,1) = pdist([squeeze(refSeed(fileNum,:,:))' ; squeeze(movSeed(fileNum,:,:))' ],'euclidean')./74; % 74 pixels = 1mm - hence converting to mm 

            % Get the connectivity value of the moving seed with respect to the reference seed
            % 74 pixels = 1000 um, 37 pixels = 500 um, 18 pixels = 250 um 
            % To note that this value holds after the connectivity maps
            % have been brought back to the size of the reference green map
            seedRad = 18; mov = squeeze(movSeed(fileNum,:,:))';
            seedVal = fliplr(round(mov));
            connVal(:,:,fileNum) = (rsConnMatrix(seedVal(1)-seedRad:seedVal(1)+seedRad,seedVal(2)-seedRad:seedVal(2)+seedRad));
        end
        
        if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology'],'dir');[~,~] = mkdir(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology']); end 
        save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\distance_conn.mat'] ,'refSeed','movSeed','distSites','connVal');

        distSitesAll{iDate,1} = distSites;
        connValsAll{iDate,1}  = connVal; 
        refSiteAll{iDate,1}   = refSeed; 
        movSiteAll{iDate,1}   = movSeed;

    else
        allVals =  load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\distance_conn.mat']);
        distSitesAll{iDate,1} = allVals.distSites;
        connValsAll{iDate,1}  = allVals.connVal; 
        refSiteAll{iDate,1}   = allVals.refSeed; 
        movSiteAll{iDate,1}   = allVals.movSeed;
        clear allVals; 
    end
end
end