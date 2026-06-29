function [allDates,allRuns,refDate,refDir,refImageName, serverPath,greenMapRef] = getMonkeyParamsRS(monkeyName,commonDir,hemisphere)
% Initializes the monkey parameters such as experiment dates, runs,
% reference dates, blood vessel maps, folder locations


if strcmp(monkeyName,'CharlieSheen')
    allDates     = ['06_22_2021'; '08_02_2021'; '08_31_2021'; '09_28_2021';'01_11_2022'; '04_04_2022'];
    allRuns      = {['run02';'run03']; ['run00'; 'run01';'run02'];['run00'; 'run01'; 'run02'];['run00'; 'run01' ; 'run02'];...
        ['run00'; 'run01' ; 'run02'];['run00'; 'run01' ; 'run02'];};
    refDate      = '08_31_2021';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    refImageName = 'Charlie Sheen Combined Green 08_31_2021';
    serverPath   = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\302-19_CharlieSheen_SqM\Left Hemisphere\';

elseif strcmp(monkeyName,'Bordeaux')
    allDates     = ['01_06_2020'; '02_10_2020']; % Bordeaux
    allRuns      = {['run01' ;'run03']; ['run00'; 'run01']}; 
    refDate      = '02_10_2020';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    serverPath   = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\Euthanized Animal Data\15-18_Bordeaux_SqM\Left Hemisphere\';
    refImageName = 'Combined Green Edited';

elseif strcmp(monkeyName,'Whiskey')
    allDates     = ['04_25_2022'; '05_09_2022'; '06_28_2022'; '08_08_2022' ];
    allRuns      = {['run00';'run01';'run02']; ['run00'; 'run01';'run02']; ['run00';'run01';'run02']; ['run00';'run01';'run02']};
    refDate      = '05_09_2022';
    refImageName = 'Combined Green';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    serverPath   = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\303-19_Whiskey_SqM\Left Hemisphere\';
end

% Get the master green image 
if exist([refDir refImageName '.bmp'],'file') == 0 % Make sure to get what the reference run is outside the function
    greenMapRef = imread([refDir refImageName '.png']);
else
    greenMapRef = imread([refDir refImageName '.bmp']);
end
greenMapRef  = greenMapRef(:,:,1);
end