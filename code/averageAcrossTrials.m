function [ data, figHandle ] = averageAcrossTrials( observerID, dateID, sessionName, trials, varargin )
% Load and average pupil responses across trials
%
% Syntax:
%   [ data, figHandle ] = averageAcrossTrials( observerID, dateID, sessionName, trials, varargin )
%
% Description:
%   Loads, cleans, and averages pupil responses obtained as part of the
%   Metropsis Glare experiment.
%
%   If there is a trial from a given session that is judged to be "bad" for
%   some reason, it may be excluded from the analysis by changing the name
%   of the pupil.mat data file for that trial. For example:
%
%       trial_41_pupil.mat  -->  trial_41_pupil_BAD.mat
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
%  'nFramesBaseline'      - Scalar. The number of frames at the start of a
%                           the trial to be used as the baseline prior to
%                           the start of the stimulus.
%  'createPlot'           - Logical. Do we create a plot of the data?
%  'plotColors','plotLabels' - Cell array. Colors and labels for the trial
%                           types that are shown in the plot.
%
% Outputs:
%   data                  - Cell array. Each entry in the cell is a matrix
%                           of dimensions [a,b], where a is the number of
%                           trials, and b is the number of time samples per
%                           trial. The matrix will contain nans for
%                           censored time points with bad ellipse fits, or
%                           if there are fewer trials of one type or
%                           another. There is one cell entry in data for
%                           each cell entry passed in the trials variable.
%   figHandle             - Handle to the plot.
%
% Examples:
%{
    observerID = 'BRIANA HAGGERTY';
    dateID = '2020-12-14';
    sessionName = {'session_1','session_2','session_3','session_4'};
    experimentName = 'pupilGlare_01';

    for ss=1:4
        T = readtable(sprintf(['/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/MTRP_data/Exp_002GN/Subject_' observerID '/' observerID '_%d.txt'],ss));
        trials{ss}(strcmp(T.Condition,'Glow'))=1;
        trials{ss}(strcmp(T.Condition,'Halo'))=2;
        trials{ss}(strcmp(T.Condition,'Uniform'))=3;
    end

    [ data ] = averageAcrossTrials( observerID, dateID, sessionName, trials, 'experimentName', experimentName );
%}
%{
    observerID = 'GLAR_01';
    dateID = '2020-12-22';
    sessionName = {'session_1','session_2','session_3','session_4'};
    experimentName = 'pupilGlare_01';

    for ss=1:4
        T = readtable(sprintf(['/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/MTRP_data/Exp_002GN/Subject_' observerID '/' observerID '_%d.txt'],ss));
        trials{ss}(strcmp(T.Condition,'Glow'))=1;
        trials{ss}(strcmp(T.Condition,'Halo'))=2;
        trials{ss}(strcmp(T.Condition,'Uniform'))=3;
    end

    [ data ] = averageAcrossTrials( observerID, dateID, sessionName, trials, 'experimentName', experimentName, 'rmseThresh',1, 'blinkFrameBuffer',4 );
%}


%% input parser
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('observerID',@ischar);
p.addRequired('dateID',@ischar);
p.addRequired('sessionName',@(x)(isvector(x) | iscell(x)));
p.addRequired('trials',@(x)(isvector(x) | iscell(x)));

% Optional
p.addParameter('experimentName','pupilGlare_01',@ischar);
p.addParameter('dropBoxBaseDir',getpref('mtrpGlarePupil','dropboxBaseDir'),@ischar);
p.addParameter('dataDir',fullfile('MTRP_data'),@ischar);
p.addParameter('processingDir',fullfile('MTRP_processing'),@ischar);
p.addParameter('rmseThresh',1.5,@isscalar);
p.addParameter('blinkFrameBuffer',4,@isscalar);
p.addParameter('nFramesBaseline',45,@isscalar);
p.addParameter('createPlot',true,@islogical);
p.addParameter('verbose',true,@islogical);
p.addParameter('plotColors',{'r',[0.5 0.5 0.5],'k'},@iscell);
p.addParameter('plotLabels',{'glow','halo','uniform'},@iscell);

% parse
p.parse(observerID, dateID, sessionName, trials, varargin{:})


%% Determine the number of sessions to process
% Also adjust passed variable types
if iscell(sessionName)
    nSessions = length(sessionName);
else
    nSessions = 1;
    sessionName = {sessionName};
    trials = {trials};
end

%% Prepare the variables and figures
figHandle = [];

%% Get the trialTypes from the first session
trialTypes = unique(trials{1});

%% Create a variable to hold the data
data = {};

%% Average pupil responses by trial type
for tt = 1:length(trialTypes)
    
    trialData = [];
    
    %% Loop over sessions
    for ss = 1:nSessions
        
        % Get the trials for this session
        idx = find(trials{ss} == trialTypes(tt));
        
        for ii = 1:length(idx)
            
            % Identify this pupilData file
            trialName = sprintf('trial_%02d_pupil.mat',idx(ii));
            pupilDataFile = fullfile(...
                p.Results.dropBoxBaseDir,...
                p.Results.processingDir,...
                p.Results.experimentName,...
                observerID,dateID,sessionName{ss},...
                trialName );
            
            % Check if the file exists
            if isfile(pupilDataFile)
                
                % Load it and add it to the data array
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
                trialData(end+1,:) = area;
                
            else
                
                % The file does not exist. Report this
                if p.Results.verbose
                    str = ['Pupil file does not exist; skipping: ' pupilDataFile '\n'];
                    fprintf(str);
                end
                
            end
        end
        
    end
    
    % Set the initial baseline of the response to a value of zero
    yMean = nanmean(trialData);
    baseAdjust = nanmean(yMean(1:p.Results.nFramesBaseline));
    trialData = trialData - baseAdjust;
    
    %% Create a plot if requested
    if p.Results.createPlot
        
        % Create the figure
        if tt == 1
            figHandle = figure();
        end
        
        % Get the mean and SEM of the response for this trial type
        yMean = nanmean(trialData);
        ySamples = sum(~isnan(trialData));
        ySEM = nanstd(trialData) ./ sqrt(ySamples);
        
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
    
    %% Add the trialData to the data variable
    data{tt} = trialData;
    
end

end