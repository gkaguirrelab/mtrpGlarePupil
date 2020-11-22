function [ output ] = averageAcrossTrials( observerID, dateID, sessionName, varargin )
% Load and average pupil responses across trials
%
% Syntax:
%   output = myFunc(input)
%
% Description:
%   Lorem ipsum dolor sit amet, consectetur adipiscing elit. Aenean euismod
%   nulla a tempor scelerisque. Maecenas et lobortis est. Donec et turpis
%   sem. Sed fringilla in metus ut malesuada. Pellentesque nec eros
%   efficitur, pellentesque nisl vel, dapibus felis. Morbi eu gravida enim.
%   Sed sodales ipsum eget finibus dapibus. Fusce sagittis felis id orci
%   egestas, non convallis neque porttitor. Proin ut mi augue. Cras posuere
%   diam at purus dignissim, vel vestibulum tellus ultrices
%
% Inputs:
%   observerID            - Char vector.
%   dateID                - Char vector.
%   sessionName           - Char vector.
%
% Optional key/value pairs:
%  'experimentName'       - Char vector.
%  'dropBoxBaseDir'       - Char vector. Default value is taken from the
%                           localHook preferences
%  'dataDir'              - Char vector.
%  'processingDir'        - Char vector.
%  'rmseThresh'           - Scalar. Frames with an ellipse fit that have an
%                           RMSE to the pupil perimeter that is higher than
%                           this value are discarded.
%  'blinkFrameBuffer'     - Scalar. This many frames on either side of
%                           a frame above the rmseThresh and discarded.
%
% Outputs:
%   baz                   - Cell. Baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz
%
% Examples:
%{
    observerID = 'GLAR_test';
    dateID = '2020-11-17';
    sessionName = 'session_1'
    experimentName = '';
    averageAcrossTrials( observerID, dateID, sessionName, 'experimentName', experimentName )
%}



%% input parser
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('observerID',@ischar);
p.addRequired('dateID',@ischar);
p.addRequired('sessionName',@ischar);

% Optional
p.addParameter('experimentName','pupilGlare_01',@ischar);
p.addParameter('dropBoxBaseDir',getpref('mtrpGlarePupil','dropboxBaseDir'),@ischar);
p.addParameter('dataDir',fullfile('MTRP_data'),@ischar);
p.addParameter('processingDir',fullfile('MTRP_processing'),@ischar);
p.addParameter('rmseThresh',0.5,@isscalar);
p.addParameter('blinkFrameBuffer',2,@isscalar);

% parse
p.parse(observerID, dateID, sessionName, varargin{:})


%% Load the trial order
% Write code here to load and process the output data file from metropsis.
trials = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
trialTypes = unique(trials);


%% Average pupil responses by trial type
for tt = 1:length(trialTypes)
    idx = find(trials == trialTypes(tt));
    
    data = [];
    
    for ii = 1:length(idx)
        
        % Load this puilData file
        trialName = sprintf('trial_%02d_pupil.mat',idx(ii));
        pupilDataFile = fullfile(...
            p.Results.dropBoxBaseDir,...
            p.Results.processingDir,...
            p.Results.experimentName,...
            observerID,dateID,sessionName,...
            trialName );
        load(pupilDataFile,'pupilData');
        
        % Obtain the pupil area vector, filter blinks, convert to % change
        area = pupilData.initial.ellipses.values(:,3);
        rmse = pupilData.initial.ellipses.RMSE;
        badIdx = rmse > p.Results.rmseThresh;   
        for rr = -p.Results.blinkFrameBuffer:p.Results.blinkFrameBuffer
            area( circshift(badIdx,rr) ) = nan;
        end
        meanArea = nanmean(area);
        area = (area - meanArea)/meanArea;

        % Add this trial data to the data array
        data(ii,:) = area;

    end

end
