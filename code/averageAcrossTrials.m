function [ data ] = averageAcrossTrials( observerID, dateID, sessionName, trials, varargin )
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
    observerID = 'GLAR_briana';
    dateID = '2020-12-01';
    sessionName = 'session_1'
    experimentName = 'pupilGlare_01';
    T = readtable('BRIANA HAGGERTY_11.txt');
    trials(strcmp(T.Condition,'Glow'))=1;
    trials(strcmp(T.Condition,'Halo'))=2;
    trials(strcmp(T.Condition,'Uniform'))=3;
    [ data ] = averageAcrossTrials( observerID, dateID, sessionName, trials, 'experimentName', experimentName );
%}



%% input parser
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('observerID',@ischar);
p.addRequired('dateID',@ischar);
p.addRequired('sessionName',@ischar);
p.addRequired('trials',@isvector);

% Optional
p.addParameter('experimentName','pupilGlare_01',@ischar);
p.addParameter('dropBoxBaseDir',getpref('mtrpGlarePupil','dropboxBaseDir'),@ischar);
p.addParameter('dataDir',fullfile('MTRP_data'),@ischar);
p.addParameter('processingDir',fullfile('MTRP_processing'),@ischar);
p.addParameter('rmseThresh',0.5,@isscalar);
p.addParameter('blinkFrameBuffer',2,@isscalar);
p.addParameter('plotColors',{'r',[0.5 0.5 0.5],'k'},@iscell);
p.addParameter('plotLabels',{'glow','halo','uniform'},@iscell);

% parse
p.parse(observerID, dateID, sessionName, trials, varargin{:})


%% Load the trial order
% Write code here to load and process the output data file from metropsis.
trialTypes = unique(trials);


%% Prepare the variables and figures
data = [];
figHandle = figure();

%% Average pupil responses by trial type
for tt = 1:length(trialTypes)

    idx = find(trials == trialTypes(tt));
        
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
        data(tt,ii,:) = area;
    end
    
    % Get the mean and SEM of the response for this trial type
    yMean = nanmean(squeeze(data(tt,:,:)));
    ySEM = nanstd(squeeze(data(tt,:,:))) / sqrt(length(idx));

    plotHandles(tt) = plot(yMean,'-','Color',p.Results.plotColors{tt});
    hold on
    plot(yMean+ySEM,':','Color',p.Results.plotColors{tt});
    plot(yMean-ySEM,':','Color',p.Results.plotColors{tt});
        
end

legend(plotHandles,p.Results.plotLabels);

end