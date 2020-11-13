%% GLAR_test_2020_11_12_pupil
%
% The video analysis pre-processing pipeline for an LDOG session.
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
pathParams.Subject = 'GLAR_test';
pathParams.Date = '2020-11-12';
pathParams.Session = 'session_1';

% The names of the videos to process
videoNameStems = {...
    'trial_01',... 
    'trial_02',...
    };

% % Stimulus properties
% sets = {[1 7],[2 8],[3 9],[4 10],[5 11], [6, 12]};
% labels = {'pupil_LightFlux_1-6Hz_RightEyeStim','pupil_RodMel_1-6Hz_RightEyeStim',...
%     'pupil_LplusS_1-6Hz_RightEyeStim', 'pupil_LightFlux_1-6Hz_LeftEyeStim',...
%     'pupil_RodMel_1-6Hz_LeftEyeStim', 'pupil_LplusS_1-6Hz_LeftEyeStim'};
% durations = [504,504,504,504,504,504];
% freqs = [1/6,1/6,1/6,1/6,1/6,1/6];

% % There is only one audio TTL pulse 
% checkCountTRs = [112 112 112 112 112 112 112 112 112 112 112 112];

% Mask bounds
glintFrameMaskSet = {...
    [193   356   240   242], ... 
    [193   354   247   244]}; 
pupilFrameMaskSet = {...
    [117   199   111   169], ...
    [106   180   116   169]}; 

pupilCircleThreshSet = [0.057, 0.063];
pupilRangeSets = {[65 80], [71 87]};
candidateThetas = {[pi],[pi]};
ellipseEccenLBUB = {[0.2 0.6],[0.2 0.7]};
glintPatchRadius = [45,45];
minRadiusProportion = [0, 0];
cutErrorThreshold = [2,2];
ellipseAreaLB = [1000, 1000];
ellipseAreaUP = [90000, 90000];
glintThreshold = [0.4, 0.4];
pupilGammaCorrection = [0.75,0.75];

%% Loop through video name stems get each video and its corresponding masks
vids = [1];
for ii = vids
    pupilCircleThresh = pupilCircleThreshSet(ii);
    pupilRange = pupilRangeSets{ii};
    videoName = {videoNameStems{ii}};
    glintFrameMask = glintFrameMaskSet{ii};
    pupilFrameMask = pupilFrameMaskSet{ii};
    % Analysis parameters
    % To adjust these parameters for a given session, use the utility:
    %{
        estimatePipelineParamsGUI('','TOME')
    %}
    % And select one of the raw data .mov files.

    sessionKeyValues = {...
        'pupilGammaCorrection', pupilGammaCorrection(ii), ...
        'startFrame',1, ...
        'nFrames', Inf, ...
        'glintFrameMask',glintFrameMask,...
        'glintGammaCorrection',0.75,...
        'glintThreshold',glintThreshold(ii),...
        'pupilFrameMask',pupilFrameMask,...
        'pupilRange',pupilRange,...
        'pupilCircleThresh',pupilCircleThresh,...
        'glintPatchRadius',glintPatchRadius(ii),...
        'candidateThetas',candidateThetas{ii},...
        'cutErrorThreshold',cutErrorThreshold(ii),...
        'radiusDivisions',50,...
        'ellipseTransparentLB',[0,0,ellipseAreaLB(ii), ellipseEccenLBUB{ii}(1), 0],...
        'ellipseTransparentUB',[1280,720,ellipseAreaUP(ii),ellipseEccenLBUB{ii}(2), pi],...
        'minRadiusProportion', minRadiusProportion(ii),...
        };

    % Call the pre-processing pipeline
    mtrpPupilPipeline(pathParams,videoName,sessionKeyValues);
    
end

% %% Call the frequency fitting pipeline
% fourierFitPipeline(pathParams,videoNameStems,sets,labels,durations,freqs);

