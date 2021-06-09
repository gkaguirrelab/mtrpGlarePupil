%% GLAR_1033_2021-06-09_pupil_session_2
%
% The video analysis pre-processing pipeline for a MTRP session.
%
% To define mask bounds, use:
%{
	glintFrameMask = defineCropMask('pupil_L+S_01.mov','startFrame',10)
	pupilFrameMask = defineCropMask('pupil_L+S_01.mov','startFrame',10)
%}
% For the glint, put a tight box around the glint. For the pupil, define a
% mask area that safely contains the pupil at its most dilated.

%% Session parameters

% Subject and session params.
pathParams.Subject = 'GLAR_1033';
pathParams.Date = '2021-06-09';
pathParams.Session = 'session_2';

%% Analysis Notes

%% Videos
vids = 1:42;

videoNameStems = {};
for ii = vids
    if ii < 10
        vidName = ['trial_0' num2str(ii)];
    else
        vidName = ['trial_' num2str(ii)];
    end
    videoNameStems{end+1} = vidName;
end

% Mask bounds
glintFrameMask = [97   307   332   283];
pupilFrameMask = [178   224   185   334];

% Pupil settings
pupilCircleThreshSet = 0.05;
pupilRangeSets = [40 50];
ellipseEccenLBUB = [0 0.88];
ellipseAreaLB = 0;
ellipseAreaUP = 90000;
pupilGammaCorrection = 0.75;

% Glint settings
glintPatchRadius = 1;
glintThreshold = 0.4;

% Control stage values (after the 3th before the 6th stage)
% Cut settings: 0 for buttom cut, pi/2 for right, pi for top, 3*pi/4 for
% left
candidateThetas = pi;
minRadiusProportion = 0.9;
cutErrorThreshold = 0.25;
%% Loop through video name stems get each video and its corresponding masks
for ii = 1:numel(vids)
    pupilCircleThresh = pupilCircleThreshSet;
    pupilRange = pupilRangeSets;
    videoName = {videoNameStems{ii}};
    % Analysis parameters
    % To adjust these parameters for a given session, use the utility:
    %{
        estimatePipelineParamsGUI('','TOME')
    %}
    % And select one of the raw data .mov files.

    sessionKeyValues = {...
        'pupilGammaCorrection', pupilGammaCorrection, ...
        'startFrame',1, ...
        'nFrames', Inf, ...
        'glintFrameMask',glintFrameMask,...
        'glintGammaCorrection',0.75,...
        'glintThreshold',glintThreshold,...
        'pupilFrameMask',pupilFrameMask,...
        'pupilRange',pupilRange,...
        'pupilCircleThresh',pupilCircleThresh,...
        'glintPatchRadius',glintPatchRadius,...
        'candidateThetas',candidateThetas,...
        'cutErrorThreshold',cutErrorThreshold,...
        'radiusDivisions',50,...
        'ellipseTransparentLB',[0,0,ellipseAreaLB, ellipseEccenLBUB(1), 0],...
        'ellipseTransparentUB',[1280,720,ellipseAreaUP,ellipseEccenLBUB(2), pi],...
        'minRadiusProportion', minRadiusProportion,...
        };

    % Call the pre-processing pipeline
    mtrpPupilPipeline(pathParams,videoName,sessionKeyValues);
    
end

