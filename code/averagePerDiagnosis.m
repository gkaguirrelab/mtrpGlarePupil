function [ data, figHandle ] = averagePerDiagnosis( diagnosis, varargin )
% Load and average pupil responses across trials for MwA, MwoA, or HaF
% subjects
%
% Syntax:
%   [ data, figHandle ] = averagePerDiagnosis( diagnosis, varargin )
%
% Description:
%   Loads, cleans, and averages pupil responses across subjects with the
%   same POEM diagnosis obtained as part of the Metropsis Glare experiment
%
% Inputs:
%   diagnosis            - Char vector. Must be one of the following: 'mwa'
%                          for migraine with aura subjects, 'mwoa' for
%                          migraine without aura subjects, or 'haf' for
%                          control subjects.
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
%   data                  - Matrix of dimensions [a,b,c], where a is the
%                           number of trial types, b is the number of
%                           trials of each type, an c is the number of time
%                           samples per trial. The matrix will contain nans
%                           for censored time points with bad ellipse fits,
%                           or if there are fewer trials of one type or
%                           another.
%   figHandle             - Handle to the plot.
%

%% set up subject session info and determine which subjects to process

% subject session information
GLAR_1001 = {'GLAR_1001', '2021-02-17', {'session_1','session_2','session_3','session_4'}};
GLAR_1002 = {'GLAR_1002', '2021-02-17', {'session_1','session_2','session_3','session_4'}};
GLAR_1003 = {'GLAR_1003', '2021-02-23', {'session_1','session_2','session_3','session_4'}};
GLAR_1004 = {'GLAR_1004', '2021-02-25', {'session_1','session_2','session_3','session_4'}};
GLAR_1005 = {'GLAR_1005', '2021-02-25', {'session_1','session_2','session_3','session_4'}};
GLAR_1006 = {'GLAR_1006', '2021-02-25', {'session_1','session_2','session_3','session_4'}};
GLAR_1007 = {'GLAR_1007', '2021-03-01', {'session_1','session_2','session_3','session_4'}};
GLAR_1008 = {'GLAR_1008', '2021-03-04', {'session_1','session_2','session_3','session_4'}};
GLAR_1009 = {'GLAR_1009', '2021-03-23', {'session_1','session_2','session_3','session_4'}};
GLAR_1010 = {'GLAR_1010', '2021-03-26', {'session_1','session_2','session_3','session_4'}};
GLAR_1011 = {'GLAR_1011', '2021-03-30', {'session_1','session_2','session_3','session_4'}};
GLAR_1012 = {'GLAR_1012', '2021-04-01', {'session_1','session_2','session_3','session_4'}};
GLAR_1013 = {'GLAR_1013', '2021-04-02', {'session_1','session_2','session_3','session_4'}};
GLAR_1014 = {'GLAR_1014', '2021-04-05', {'session_1','session_2','session_3','session_4'}};
GLAR_1015 = {'GLAR_1015', '2021-04-06', {'session_1','session_2','session_3','session_4'}};
GLAR_1016 = {'GLAR_1016', '2021-04-08', {'session_1','session_2','session_3','session_4'}};
GLAR_1017 = {'GLAR_1017', '2021-04-12', {'session_1','session_2','session_3','session_4'}};
GLAR_1020 = {'GLAR_1020', '2021-05-06', {'session_1','session_2','session_3','session_5'}};
GLAR_1021 = {'GLAR_1021', '2021-05-06', {'session_1','session_2','session_3','session_4'}};
GLAR_1022 = {'GLAR_1022', '2021-05-10', {'session_1','session_2','session_3','session_4'}};
GLAR_1023 = {'GLAR_1023', '2021-05-13', {'session_1','session_2','session_3','session_4'}};
GLAR_1025 = {'GLAR_1025', '2021-05-18', {'session_1','session_2','session_3','session_4'}};
GLAR_1026 = {'GLAR_1026', '2021-05-19', {'session_1','session_2','session_3','session_5'}};
GLAR_1028 = {'GLAR_1028', '2021-05-24', {'session_1','session_2','session_3','session_4'}};

% subject POEM category
mwaSubjects = {GLAR_1010, GLAR_1011, GLAR_1012, GLAR_1013, GLAR_1014,...
    GLAR_1015, GLAR_1016, GLAR_1017, GLAR_1020, GLAR_1021};
mwoaSubjects = {GLAR_1023, GLAR_1025, GLAR_1026, GLAR_1028};
hafSubjects = {GLAR_1001, GLAR_1002, GLAR_1003, GLAR_1004, GLAR_1005,...
    GLAR_1006, GLAR_1007, GLAR_1008, GLAR_1009, GLAR_1022};

%% input parser
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('diagnosis',@ischar);

% Optional
p.addParameter('dropBoxBaseDir',getpref('mtrpGlarePupil','dropboxBaseDir'),@ischar);
p.addParameter('createPlot',true,@islogical);
p.addParameter('plotColors',{'r',[0.5 0.5 0.5],'k'},@iscell);
p.addParameter('plotLabels',{'glow','halo','uniform'},@iscell);

% parse
p.parse(diagnosis, varargin{:})

% select POEM category
if strcmp(diagnosis,'mwa')
    pC = mwaSubjects;
    groupName = 'Average across migraine with aura subjects';
elseif strcmp(diagnosis,'mwoa')
    pC = mwoaSubjects;
    groupName = 'Average across migraine without aura subjects';
elseif strcmp(diagnosis,'haf')
    pC = hafSubjects;
    groupName = 'Average across headache-free control subjects';
else
    error('Diagnosis string input must be either mwa, mwoa, or haf.');
end

%% Prepare the variables and figures
figHandle = [];
glow = [];
halo = [];
uniform = [];

%% process trials for each subject

% Turn off a table loading warning
warnState = warning();
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

% Hide the figures as we make them
%set(0,'DefaultFigureVisible','off');

% Loop through the list of subjects
for ii = 1:length(pC)
    % setup subject session information
    subject = pC{ii};
    observerID = subject{1};
    dateID = subject{2};
    sessionName = subject{3};
    
    for ss=1:4
        fileName = fullfile(p.Results.dropBoxBaseDir,sprintf(['MTRP_data/Exp_002GN/Subject_' observerID '/' observerID '_%d.txt'],ss));
        T = readtable(fileName);
        trials{ss}(strcmp(T.Condition,'Glow'))=1;
        trials{ss}(strcmp(T.Condition,'Halo'))=2;
        trials{ss}(strcmp(T.Condition,'Uniform'))=3;
    end
    
    % get session data
    [subData] = averageAcrossTrials(observerID, dateID, sessionName, trials, 'createPlot', false);
    
    % Get the mean of the response for each trial type
    glowMean = nanmean(subData{1});
    haloMean = nanmean(subData{2});
    uniformMean = nanmean(subData{3});
    
    % add subject averages to respective POEM category matrices
    glow = [glow; glowMean];
    halo = [halo; haloMean];
    uniform = [uniform; uniformMean];
end


% Restore the warning state
warning(warnState);

%% create plot if requested

if p.Results.createPlot
    
    figHandle = figure();
    trialTypes = {glow, halo, uniform};
    
    for tt = 1:length(trialTypes)
        
        typeData = trialTypes{tt};
        
        % Get the mean and SEM of the response for each trial type
        yMean = nanmean(typeData);
        ySamples = sum(~isnan(typeData));
        ySEM = nanstd(typeData) ./ sqrt(ySamples);
        
        % Set the x temporal support
        xVals = (1:length(yMean))/60;
        % Plot
        plotHandles(tt) = plot(xVals,yMean,'-','Color',p.Results.plotColors{tt});
        hold on
        plot(xVals,yMean+ySEM,':','Color',p.Results.plotColors{tt});
        plot(xVals,yMean-ySEM,':','Color',p.Results.plotColors{tt});
        
        if tt == length(trialTypes)
            legend(plotHandles,p.Results.plotLabels);
            title(groupName,'interpreter', 'none');
            ylabel('Proportion change pupil area');
            xlabel('Time [secs]');
            ylim([-0.15 0.05])
        end
        
    end
end

end