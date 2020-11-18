function pupilPipeline(pathParams,videoNameStems, sessionKeyValues)
% Run through transparentTrack stages to perform pupil video pre-processing
%
% Syntax:
%  pupilPipeline(pathParams,videoNameStems, universalKeyValues, sessionKeyValues)
%
% Description:
%   Calls the stages of transparent track to pre-process a set of videos
%   from the MTRP project.
%
% Inputs:
%   pathParams            - Structure. Defines the path to the data and the
%                           processed output. Has the required fields:
%                           dataSourceDirRoot, dataOutputDirRoot,
%                           Approach, Protocol, Subject, Date
%   videoNameStems        - Cell array of char vectors. The names of the
%                           video files to process, omitting the '.mov'
%                           suffix.
%   sessionKeyValues      - Cell array. These key-values are passed to each
%                           of the processing stages.
%
% Outputs:
%   none
%


%% Define universal key values
% These are parameters that are common to all videos processed under the
% MTRP project.
universalKeyValues = {...
    'intrinsicCameraMatrix',[2627.0 0 338.1; 0 2628.1 246.2; 0 0 1],...
    'radialDistortionVector',[-0.3517 3.5353],...
    'eyeLaterality','left',...
    'spectralDomain','nir', ...
    'verbose',true, ...
    'audioTrackSync', true, ...
    'useParallel', true, ...
    'overwriteControlFile', true ...
    };



%% Add the DropBox information to the pathParams
% Get the DropBox base directory
dropboxBaseDir = getpref('mtrpGlarePupil','dropboxBaseDir');

% set common path params
pathParams.dataSourceDirRoot = fullfile(dropboxBaseDir,'MTRP_data', 'pupilGlare_01');
pathParams.dataOutputDirRoot = fullfile(dropboxBaseDir,'MTRP_processing');


%% Create the input and output paths
inputBaseDir = fullfile(pathParams.dataSourceDirRoot, ...
    pathParams.Subject,...
    pathParams.Date,...
    pathParams.Session);

outputBaseDir = fullfile(pathParams.dataOutputDirRoot, ...
    pathParams.Subject,...
    pathParams.Date,...
    pathParams.Session);

% If the outputBaseDir does not exist, make it
if ~exist(outputBaseDir)
    mkdir(outputBaseDir)
end

%% Turn off the warning regarding over-writing the control file
warningState = warning;
warning('off','makeControlFile:overwrittingControlFile');


%% Loop through the videos
for vv = 1:length(videoNameStems)
    
    % Report which video we will now process
    str = sprintf(['%d of %d, ' videoNameStems{vv} '\n\n'],vv,length(videoNameStems));
    fprintf(str);
    
    % The file names
    videoInFileName = fullfile(inputBaseDir,[videoNameStems{vv} '.mp4']);
    convertedVideoName = fullfile(outputBaseDir,[videoNameStems{vv} '_converted.avi']);
    grayVideoName = fullfile(outputBaseDir,[videoNameStems{vv} '_gray.avi']);
    glintFileName = fullfile(outputBaseDir,[videoNameStems{vv} '_glint.mat']);
    perimeterFileName = fullfile(outputBaseDir,[videoNameStems{vv} '_perimeter.mat']);
    correctedPerimeterFileName = fullfile(outputBaseDir,[videoNameStems{vv} '_correctedPerimeter.mat']);
    controlFileName = fullfile(outputBaseDir,[videoNameStems{vv} '_controlFile.csv']);
    pupilFileName = fullfile(outputBaseDir,[videoNameStems{vv} '_pupil.mat']);
    fit3VideoName = fullfile(outputBaseDir,[videoNameStems{vv} '_stage3fit.avi']);
    fit6VideoName = fullfile(outputBaseDir,[videoNameStems{vv} '_stage6fit.avi']);
     
    % Convert
    convertMtrpVideos(videoInFileName, convertedVideoName) 
    
    % Deinterlace
    deinterlaceVideo(convertedVideoName, grayVideoName, ...
        universalKeyValues{:},sessionKeyValues{:});
    
    % Glint
    findGlint(grayVideoName, glintFileName, ...
        universalKeyValues{:},sessionKeyValues{:});
    
    % Perimeter
    findPupilPerimeter(grayVideoName, perimeterFileName, ...
        universalKeyValues{:},sessionKeyValues{:});

    % 3rd Stage Video
    makeFitVideo(grayVideoName, fit3VideoName, ...
        'perimeterFileName',perimeterFileName,...
        'glintFileName',glintFileName,...
        universalKeyValues{:},sessionKeyValues{:});
    
    % Control
    makeControlFile(controlFileName, perimeterFileName, glintFileName, ...
        universalKeyValues{:},sessionKeyValues{:});
    
    % Correct
    applyControlFile(perimeterFileName, controlFileName, correctedPerimeterFileName, ...
        universalKeyValues{:},sessionKeyValues{:});
    
    % Fit
    fitPupilPerimeter(correctedPerimeterFileName, pupilFileName, ...
        universalKeyValues{:},sessionKeyValues{:});
    
    % 6th Stage Video
    makeFitVideo(grayVideoName, fit6VideoName, ...
        'perimeterFileName',correctedPerimeterFileName,...
        'pupilFileName',pupilFileName,...
        'glintFileName',glintFileName,...
        'fitLabel', 'initial', ...
        universalKeyValues{:},sessionKeyValues{:});
    
end

%% Restore the warning state
warning(warningState);


end

