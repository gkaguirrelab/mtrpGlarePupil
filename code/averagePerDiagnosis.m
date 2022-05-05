function [ data, figHandles ] = averagePerDiagnosis( varargin )
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
GLAR_1029 = {'GLAR_1029', '2021-05-26', {'session_1','session_2','session_3','session_4'}};
GLAR_1030 = {'GLAR_1030', '2021-05-28', {'session_1','session_2','session_3','session_4'}};
GLAR_1031 = {'GLAR_1031', '2021-06-01', {'session_1','session_2','session_3','session_4'}};
GLAR_1032 = {'GLAR_1032', '2021-06-04', {'session_1','session_2','session_4','session_5'}};
GLAR_1033 = {'GLAR_1033', '2021-06-09', {'session_1','session_2','session_3','session_4'}};
GLAR_1034 = {'GLAR_1034', '2021-06-22', {'session_1','session_2','session_3','session_4'}};
GLAR_1035 = {'GLAR_1035', '2021-06-25', {'session_1','session_2','session_3','session_4'}};
GLAR_1036 = {'GLAR_1036', '2021-07-01', {'session_1','session_2','session_3','session_4'}};
GLAR_1037 = {'GLAR_1037', '2021-07-07', {'session_1','session_2','session_3','session_4'}};
GLAR_1038 = {'GLAR_1038', '2021-08-17', {'session_1','session_2','session_3','session_4'}};
GLAR_1040 = {'GLAR_1040', '2021-10-07', {'session_1','session_2','session_3','session_4'}};
GLAR_1041 = {'GLAR_1042', '2021-10-13', {'session_1','session_2','session_3','session_4'}};
GLAR_1042 = {'GLAR_1042', '2021-10-13', {'session_1','session_2','session_3','session_4'}};
GLAR_1043 = {'GLAR_1043', '2021-11-02', {'session_1','session_2','session_3','session_4'}};
GLAR_1044 = {'GLAR_1044', '2021-11-03', {'session_1','session_2','session_3','session_4'}};
GLAR_1045 = {'GLAR_1045', '2021-11-16', {'session_1','session_2','session_3','session_4'}};
GLAR_1047 = {'GLAR_1047', '2021-12-09', {'session_1','session_2','session_3','session_4'}};
GLAR_1048 = {'GLAR_1048', '2021-12-14', {'session_1','session_2','session_3','session_4'}};
GLAR_1049 = {'GLAR_1049', '2022-02-01', {'session_1','session_2','session_3','session_5'}};
GLAR_1050 = {'GLAR_1050', '2022-02-03', {'session_1','session_2','session_3','session_4'}};
GLAR_1051 = {'GLAR_1051', '2022-04-23', {'session_1','session_2','session_3','session_4'}};


% subject POEM category
% Migraine with Aura
subjectSets{1} = {GLAR_1010, GLAR_1011, GLAR_1012, GLAR_1013, GLAR_1014,...
    GLAR_1015, GLAR_1016, GLAR_1017, GLAR_1020, GLAR_1021, GLAR_1034,...
    GLAR_1035, GLAR_1038, GLAR_1048, GLAR_1051};
% Migraine without aura
subjectSets{2} = {GLAR_1023, GLAR_1025, GLAR_1026, GLAR_1028, GLAR_1030,...
    GLAR_1033, GLAR_1040, GLAR_1041, GLAR_1042, GLAR_1043, GLAR_1044,...
    GLAR_1045, GLAR_1047, GLAR_1049, GLAR_1050};
% Headache free
subjectSets{3} = {GLAR_1001, GLAR_1002, GLAR_1003, GLAR_1004, GLAR_1005,...
    GLAR_1006, GLAR_1007, GLAR_1008, GLAR_1009, GLAR_1022, GLAR_1029,...
    GLAR_1031, GLAR_1032, GLAR_1036, GLAR_1037};

%% input parser
p = inputParser; p.KeepUnmatched = false;

% Optional
p.addParameter('dropBoxBaseDir',getpref('mtrpGlarePupil','dropboxBaseDir'),@ischar);
p.addParameter('createPlot',true,@islogical);
p.addParameter('plotColors',{'r',[0 0 1],'k'},@iscell);
p.addParameter('plotLabels',{'Glow','Halo','Uniform'},@iscell);
p.addParameter('diagnosis',{'Migraine with aura','Migraine without aura','Headache free'},@iscell);

% parse
p.parse( varargin{:})


%% Prepare the variables and figures
figHandle = [];
glow = [];
halo = [];
uniform = [];

%% Load the data

% Turn off a table loading warning
warnState = warning();
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

% Loop through the list of subjects
for dd = 1:length(p.Results.diagnosis)
    pC = subjectSets{dd};
    
    glow = [];
    halo = [];
    uniform = [];
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
        [subData] = averageAcrossTrials(observerID, dateID, sessionName, trials, 'createPlot', false, 'verbose', false);
        
        % Get the mean of the response for each trial type
        glowMean = nanmean(subData{1});
        haloMean = nanmean(subData{2});
        uniformMean = nanmean(subData{3});
        
        % add subject averages to respective POEM category matrices
        glow = [glow; glowMean];
        halo = [halo; haloMean];
        uniform = [uniform; uniformMean];
    end
    
    data{dd} = {glow, halo, uniform};
end

% Restore the warning state
warning(warnState);


%% Obtain the data means
meanData = cellfun(@(x) [nanmean(x{1},2), nanmean(x{2},2), nanmean(x{3},2)],data,'UniformOutput',false);


%% Provide a table of the mean responses
T = array2table(cell2mat(cellfun(@(x) 100.*mean(x)',meanData,'UniformOutput',false)));
T.Properties.RowNames=p.Results.plotLabels;
T.Properties.VariableNames=p.Results.diagnosis;

fprintf('Table of mean percent change pupil response by group and stimulus:\n')
T
fprintf('\n')

%% Report t-test for all subject effects
[~,pVal,ci,stats] = ttest(...
    -100.*cell2mat(cellfun(@(x) x(:,1)',meanData,'UniformOutput',false)),...
    -100.*cell2mat(cellfun(@(x) x(:,3)',meanData,'UniformOutput',false))...
    );
fprintf('All subjects, glow-uniform: t(df) = %2.2f (%d); p = %2.3f; mean (95CI) = %2.1f (%2.1f,%2.1f)  \n',stats.tstat,stats.df,pVal,mean(ci),ci(1),ci(2));
[~,pVal,ci,stats] = ttest(...
    -100.*cell2mat(cellfun(@(x) x(:,2)',meanData,'UniformOutput',false)),...
    -100.*cell2mat(cellfun(@(x) x(:,3)',meanData,'UniformOutput',false))...
    );
fprintf('All subjects, halo-uniform: t(df) = %2.2f (%d); p = %2.3f; mean (95CI) = %2.1f (%2.1f,%2.1f)  \n',stats.tstat,stats.df,pVal,mean(ci),ci(1),ci(2));
[~,pVal,ci,stats] = ttest(...
    -100.*cell2mat(cellfun(@(x) x(:,1)',meanData,'UniformOutput',false)),...
    -100.*cell2mat(cellfun(@(x) x(:,2)',meanData,'UniformOutput',false))...
    );
fprintf('All subjects, glow-halo: t(df) = %2.2f (%d); p = %2.3f; mean (95CI) = %2.1f (%2.1f,%2.1f)  \n',stats.tstat,stats.df,pVal,mean(ci),ci(1),ci(2));


%% Report glow-uniform effects
if p.Results.createPlot
    figure
end
divs=15;
for dd = 1:length(p.Results.diagnosis)
    yVals{dd} = -100.*(meanData{dd}(:,1)-meanData{dd}(:,3));
    xVals = repmat(dd,1,length(yVals{dd}));
    [N,~,bin] = histcounts(yVals{dd},'BinWidth',0.5);
    for xx=1:length(N)
        count = N(xx);
        idx=find(bin==xx);
        xVals(idx)=xVals(idx)+linspace(-(count-1)/divs,(count-1)/divs,count);
    end
    if p.Results.createPlot
        scatter(xVals,yVals{dd},200,...
            'MarkerEdgeColor','none',...
            'MarkerFaceColor','k',...
            'MarkerFaceAlpha',0.5);
        hold on
        plot([dd-0.25,dd+0.25],[mean(yVals{dd}),mean(yVals{dd})],'-r')
    end
end
if p.Results.createPlot
    plot([0.5 3.5],[0 0],':k');
    xlim([0.5 3.5]);
    xticks([1 2 3]);
    xticklabels(p.Results.diagnosis)
    ylabel('Pupil response glow - uniform [%∆]','FontSize',16);
end

%  t-tests
[~,pVal,ci,stats] = ttest(cell2mat(yVals'));
fprintf('All subjects, glow-uniform: t(df) = %2.2f (%d); p = %2.3f; mean (95CI) = %2.1f (%2.1f,%2.1f)  \n',stats.tstat,stats.df,pVal,mean(ci),ci(1),ci(2));
[~,pVal,ci,stats] = ttest2(yVals{1},yVals{3});
fprintf('MwA vs HaF, glow-uniform: t(df) = %2.2f (%d); p = %2.3f; mean (95CI) = %2.1f (%2.1f,%2.1f)  \n',stats.tstat,stats.df,pVal,mean(ci),ci(1),ci(2));
[~,pVal,ci,stats] = ttest2(yVals{1},yVals{2});
fprintf('MwA vs MwoA, glow-uniform: t(df) = %2.2f (%d); p = %2.3f; mean (95CI) = %2.1f (%2.1f,%2.1f)  \n',stats.tstat,stats.df,pVal,mean(ci),ci(1),ci(2));


%% Report [illusory brightness]-uniform effects
if p.Results.createPlot
    figure
end
divs=15;
for dd = 1:length(p.Results.diagnosis)
    yVals{dd} = -100.*(0.5*(meanData{dd}(:,1) + meanData{dd}(:,2))-meanData{dd}(:,3));
    xVals = repmat(dd,1,length(yVals{dd}));
    [N,~,bin] = histcounts(yVals{dd},'BinWidth',0.5);
    for xx=1:length(N)
        count = N(xx);
        idx=find(bin==xx);
        xVals(idx)=xVals(idx)+linspace(-(count-1)/divs,(count-1)/divs,count);
    end
    if p.Results.createPlot
        scatter(xVals,yVals{dd},200,...
            'MarkerEdgeColor','none',...
            'MarkerFaceColor','k',...
            'MarkerFaceAlpha',0.5);
        hold on
        plot([dd-0.25,dd+0.25],[mean(yVals{dd}),mean(yVals{dd})],'-r')
    end
end
if p.Results.createPlot
    plot([0.5 3.5],[0 0],':k');
    xlim([0.5 3.5]);
    xticks([1 2 3]);
    xticklabels(p.Results.diagnosis)
    ylabel('Pupil response illusory bright - uniform [%∆]','FontSize',16);
end

% Report t-tests
[~,pVal,ci,stats] = ttest(cell2mat(yVals'));
fprintf('All subjects, bright-uniform: t(df) = %2.2f (%d); p = %2.3f; mean (95CI) = %2.1f (%2.1f,%2.1f)  \n',stats.tstat,stats.df,pVal,mean(ci),ci(1),ci(2));
[~,pVal,ci,stats] = ttest2(yVals{1},yVals{3});
fprintf('MwA vs HaF, bright-uniform: t(df) = %2.2f (%d); p = %2.3f; mean (95CI) = %2.1f (%2.1f,%2.1f)  \n',stats.tstat,stats.df,pVal,mean(ci),ci(1),ci(2));
[~,pVal,ci,stats] = ttest2(yVals{1},yVals{2});
fprintf('MwA vs MwoA, bright-uniform: t(df) = %2.2f (%d); p = %2.3f; mean (95CI) = %2.1f (%2.1f,%2.1f)  \n',stats.tstat,stats.df,pVal,mean(ci),ci(1),ci(2));



%% create time-series plot if requested

if p.Results.createPlot
    
    for dd = 1:length(p.Results.diagnosis)
        figHandles{dd} = figure();
                
        trialTypes = data{dd};
        plotHandles = [];

        for tt = 1:length(trialTypes)
            
            typeData = trialTypes{tt};
            
            % Get the mean and SEM of the response for each trial type
            yMean = 100.*nanmean(typeData);
            ySamples = sum(~isnan(typeData));
            ySEM = 100.*nanstd(typeData) ./ sqrt(ySamples);
            
            % Set the x temporal support
            xVals = (1:length(yMean))/60;
            
            % Plot
            pl = subplot(1,1,1);
            xx = 0:.025:xVals(end);
            yy = csaps(xVals,yMean,0.9999,xx);
            yySEM = csaps(xVals,ySEM,0.9999,xx);
            plotHandles(tt) = plot(xx,yy,'-','Color',p.Results.plotColors{tt},'LineWidth',2);
            pl.Box = 'off';
            hold on
            time = [xx, fliplr(xx)];
            inBetween = [yy+yySEM, fliplr(yy-yySEM)];
            patch(time,inBetween,p.Results.plotColors{tt},'EdgeColor','none','FaceAlpha',0.08);
            patch(time,inBetween,p.Results.plotColors{tt},'EdgeColor','none','FaceAlpha',0.08);
            
            if tt == length(trialTypes)
                lgd = legend(plotHandles,p.Results.plotLabels,'Location', 'southeast');
                lgd.FontSize = 16;
                title([p.Results.diagnosis{dd}], 'FontSize', 16);
                ylabel('Percent change pupil area', 'FontSize', 16);
                xlabel('Time [secs]', 'FontSize', 16);
                ylim([-15 5])
            end
            
        end
    end
end

end