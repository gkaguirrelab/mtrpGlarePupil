function [ data, figHandle ] = averageAcrossTrials( observerID, dateID, sessionName, trials, varargin )
% Load and average pupil responses across trials
%
% Syntax:
%   [ data, figHandle ] = averageAcrossTrials( observerID, dateID, sessionName, trials, varargin )
%
% Description:
%   Loads, cleans, and averages pupil responses obtained as part of the
%   Metropsis Glare experiment
%
% Inputs:
%   observerID            - Char vector.
%   dateID                - Char vector.
%   sessionName           - Char vector.
%   trials                - Numeric vector. Identifies the stimulus type
%                           for each trial. Make sure that these identities
%                           match the labels.
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
%  'createPlot'           - Logical. Do we create a plot of the data?
%  'plotColors','plotLabels' - Cell array. Colors and labels for the trial
%                           types that are shown in the plot.
%
% Outputs:
%   data                  - Matrix of dimensions [a,b,c], where a is the
%                           number of trial types, b is the number of
%                           trials of each type, an c is the number of time
%                           samples per trial. The matrix will contain nans
%                           for censored time points with bad ellipse fits,
%                           or if there are fewer trials of one type or
%                           another.
%   figHandle             - Handle to the plot.
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
%{
    observerID = 'GLAR_briana';
    dateID = '2020-12-01';
    sessionName = 'session_2'
    experimentName = 'pupilGlare_01';
    T = readtable('BRIANA HAGGERTY_12.txt');
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
p.addParameter('createPlot',true,@islogical);
p.addParameter('plotColors',{'r',[0.5 0.5 0.5],'k'},@iscell);
p.addParameter('plotLabels',{'glow','halo','uniform'},@iscell);

% parse
p.parse(observerID, dateID, sessionName, trials, varargin{:})


%% Load the trial order
% Write code here to load and process the output data file from metropsis.
trialTypes = unique(trials);


%% Prepare the variables and figures
data = [];
figHandle = [];

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
    
    
    %% Create a plot if requested
    if p.Results.createPlot
        
        % Create the figure
        if tt == 1
            figHandle = figure();
        end
        
        % Get the mean and SEM of the response for this trial type
        yMean = nanmean(squeeze(data(tt,:,:)));
        ySEM = nanstd(squeeze(data(tt,:,:))) / sqrt(length(idx));
        
        % Set the first second of the response to a value of zero
        yMean = yMean - nanmean(yMean(1:60));
        
        % Set the x temporal support
        xVals = (1:length(yMean))/60;
        
        % Plot
        plotHandles(tt) = plot(xVals,yMean,'-','Color',p.Results.plotColors{tt});
        hold on
        plot(xVals,yMean+ySEM,':','Color',p.Results.plotColors{tt});
        plot(xVals,yMean-ySEM,':','Color',p.Results.plotColors{tt});
        
        if tt == length(trialTypes)
            legend(plotHandles,p.Results.plotLabels);
            title([observerID ' - ' dateID ' - ' sessionName],'interpreter', 'none');
            ylabel('Proportion change pupil area');
            xlabel('Time [secs]');
        end
        
    end
    
end

end