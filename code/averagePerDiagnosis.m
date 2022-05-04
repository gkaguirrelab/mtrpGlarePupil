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

% Hide the figures as we make them
%set(0,'DefaultFigureVisible','off');

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


%% calculate the bootstrap auc
myFun = @(x) sum(nanmean(x));
nBoots = 1000;

%% Obtain the data means
meanData = cellfun(@(x) [nanmean(x{1},2), nanmean(x{2},2), nanmean(x{3},2)],data,'UniformOutput',false);

% Make a plot of mean effects
figure
divs=15;
binWidth=1;
for dd = 1:length(p.Results.diagnosis)
    yVals = -100.*(meanData{dd}(:,1)-meanData{dd}(:,3));
    xVals = repmat(dd,1,length(yVals));
    [N,~,bin] = histcounts(yVals,'BinWidth',0.5);
    for xx=1:length(N)
        count = N(xx);
        idx=find(bin==xx);
        xVals(idx)=xVals(idx)+linspace(-(count-1)/divs,(count-1)/divs,count);
    end
    scatter(xVals,yVals,200,...
        'MarkerEdgeColor','none',...
        'MarkerFaceColor','k',...
        'MarkerFaceAlpha',0.5);
    hold on
    plot([dd-0.25,dd+0.25],[mean(yVals),mean(yVals)],'-r')
end
plot([0.5 3.5],[0 0],':k');
xlim([0.5 3.5]);
xticks([1 2 3]);
xticklabels(p.Results.diagnosis)
ylabel('Pupil response glow - uniform [%∆]','FontSize',16);

%% Report means
bootData = [data{1}{1}; data{2}{1}; data{3}{1}];
bootstat = sort(bootstrp(nBoots,myFun,bootData));
fprintf('Mean glow = %2.2f \n',mean(bootstat));
bootData = [data{1}{2}; data{2}{2}; data{3}{2}];
bootstat = sort(bootstrp(nBoots,myFun,bootData));
fprintf('Mean halo = %2.2f \n',mean(bootstat));
bootData = [data{1}{3}; data{2}{3}; data{3}{3}];
bootstat = sort(bootstrp(nBoots,myFun,bootData));
fprintf('Mean uniform = %2.2f \n',mean(bootstat));



% All subjects, glow vs. uniform
bootData = [data{1}{1}-data{1}{3}; data{2}{1}-data{2}{3}; data{3}{1}-data{3}{3}];
bootstat = sort(bootstrp(nBoots,myFun,bootData));
n = size(bootData,1);
t = mean(bootstat) / std(bootstat);
df = n-1;
pVal = tpdf(t,df);
fprintf('glow vs. uniform, all subjects: mean (95 CI) = %2.2f [%2.2f to %2.2f], t(df) = %2.1f (%d), p=%2.3f.\n',mean(bootstat),bootstat(nBoots*0.025),bootstat(nBoots*0.975),t,df,pVal);

% All subjects, glow vs. halo
bootData = [data{1}{1}-data{1}{2}; data{2}{1}-data{2}{2}; data{3}{1}-data{3}{2}];
bootstat = sort(bootstrp(nBoots,myFun,bootData));
n = size(bootData,1);
t = mean(bootstat) / std(bootstat);
df = n-1;
pVal = tpdf(t,df);
fprintf('glow vs. halo, all subjects: mean (95 CI) = %2.2f [%2.2f to %2.2f], t(df) = %2.1f (%d), p=%2.3f.\n',mean(bootstat),bootstat(nBoots*0.025),bootstat(nBoots*0.975),t,df,pVal);

% All subjects, halo vs uniform
bootData = [data{1}{2}-data{1}{3}; data{2}{2}-data{2}{3}; data{3}{2}-data{3}{3}];
bootstat = sort(bootstrp(nBoots,myFun,bootData));
n = size(bootData,1);
t = mean(bootstat) / std(bootstat);
df = n-1;
pVal = tpdf(t,df);
fprintf('halo vs. uniform, all subjects: mean (95 CI) = %2.2f [%2.2f to %2.2f], t(df) = %2.1f (%d), p=%2.3f.\n',mean(bootstat),bootstat(nBoots*0.025),bootstat(nBoots*0.975),t,df,pVal);


% MwA, glow vs uniform
bootData = [data{1}{1}-data{1}{3}];
bootstat = sort(bootstrp(nBoots,myFun,bootData));
n = size(bootData,1);
t = mean(bootstat) / std(bootstat);
df = n-1;
pVal = tpdf(t,df);
fprintf('glow vs. uniform, MwA: mean (95 CI) = %2.2f [%2.2f to %2.2f], t(df) = %2.1f (%d), p=%2.3f.\n',mean(bootstat),bootstat(nBoots*0.025),bootstat(nBoots*0.975),t,df,pVal);

% HaF, glow vs uniform
bootData = [data{2}{1}-data{2}{3}];
bootstat = sort(bootstrp(nBoots,myFun,bootData));
n = size(bootData,1);
t = mean(bootstat) / std(bootstat);
df = n-1;
pVal = tpdf(t,df);
fprintf('glow vs. uniform, MwoA: mean (95 CI) = %2.2f [%2.2f to %2.2f], t(df) = %2.1f (%d), p=%2.3f.\n',mean(bootstat),bootstat(nBoots*0.025),bootstat(nBoots*0.975),t,df,pVal);

% HaF, glow vs uniform
bootData = [data{3}{1}-data{3}{3}];
bootstat = sort(bootstrp(nBoots,myFun,bootData));
n = size(bootData,1);
t = mean(bootstat) / std(bootstat);
df = n-1;
pVal = tpdf(t,df);
fprintf('glow vs. uniform, HaF: mean (95 CI) = %2.2f [%2.2f to %2.2f], t(df) = %2.1f (%d), p=%2.3f.\n',mean(bootstat),bootstat(nBoots*0.025),bootstat(nBoots*0.975),t,df,pVal);


% MwA, structured vs uniform
bootData = [(data{1}{1}+data{1}{2})./2-data{1}{3}];
bootstat = sort(bootstrp(nBoots,myFun,bootData));
n = size(bootData,1);
t = mean(bootstat) / std(bootstat);
df = n-1;
pVal = tpdf(t,df);
fprintf('structured vs. uniform, MwA: mean (95 CI) = %2.2f [%2.2f to %2.2f], t(df) = %2.1f (%d), p=%2.3f.\n',mean(bootstat),bootstat(nBoots*0.025),bootstat(nBoots*0.975),t,df,pVal);

% Non-MwA, structured vs uniform
bootData = [(data{2}{1}+data{2}{2})./2-data{2}{3}; (data{3}{1}+data{3}{2})./2-data{3}{3}];
bootstat = sort(bootstrp(nBoots,myFun,bootData));
n = size(bootData,1);
t = mean(bootstat) / std(bootstat);
df = n-1;
pVal = tpdf(t,df);
fprintf('structured vs. uniform, non-MwA: mean (95 CI) = %2.2f [%2.2f to %2.2f], t(df) = %2.1f (%d), p=%2.3f.\n',mean(bootstat),bootstat(nBoots*0.025),bootstat(nBoots*0.975),t,df,pVal);


% Interaction, migraine w aura vs. control with glow vs. uniform
mwaData = data{1}{1}-data{1}{3};
mwaN = size(mwaData,1);
mwaBoot = sort(bootstrp(nBoots,myFun,mwaData));
mwaSD = std(mwaBoot)/sqrt(1/mwaN);
hafData = data{3}{1}-data{3}{3};
hafN = size(hafData,1);
hafBoot = sort(bootstrp(nBoots,myFun,hafData));
hafSD = std(hafBoot)/sqrt(1/hafN);
tStat = (mean(mwaBoot)-mean(hafBoot)) / sqrt( mwaSD^2/mwaN + hafSD^2/hafN   );
df = mwaN+hafN-2;
pStat = tpdf(tStat,df);
fprintf('MwA vs. HaF x glow vs. uniform:  t-value (df), p = %2.2f (%d), %2.2f \n',tStat,df,pStat);

% Interaction, migraine w aura vs. control with structured vs. uniform
mwaData = (data{1}{1}+data{1}{2})./2-data{1}{3};
mwaN = size(mwaData,1);
mwaBoot = sort(bootstrp(nBoots,myFun,mwaData));
mwaSD = std(mwaBoot)/sqrt(1/mwaN);
hafData = (data{3}{1}+data{3}{2})./2-data{3}{3};
hafN = size(hafData,1);
hafBoot = sort(bootstrp(nBoots,myFun,hafData));
hafSD = std(hafBoot)/sqrt(1/hafN);
tStat = (mean(mwaBoot)-mean(hafBoot)) / sqrt( mwaSD^2/mwaN + hafSD^2/hafN   );
df = mwaN+hafN-2;
pStat = tpdf(tStat,df);
fprintf('POST-HOC MwA vs. HaF x structured vs. uniform:  t-value (df), p = %2.2f (%d), %2.2f \n',tStat,df,pStat);


%% Fit the response
%{
trialTypes = data{1};

temporalFit = tfeTPUP('verbosity','none');

% Set up some default params
[initialParams,vlbParams,vubParams] = temporalFit.defaultParams;
defaultParamsInfo.nInstances = 1;

vlbParams.paramMainMatrix(1) = -500;
vubParams.paramMainMatrix = [1000 500 30 0 0 0];

for tt = 1:length(trialTypes)
    
    typeData = trialTypes{tt};
    
    % create the packet for fitting
    thePacket.response.values = nanmean(typeData)*100;
    thePacket.response.timebase = (0:1/60:(length(thePacket.response.values)-1)/60)*1000;
    thePacket.stimulus.timebase = thePacket.response.timebase;
    thePacket.stimulus.values = zeros(size(thePacket.response.timebase));
    thePacket.stimulus.values(30:90)=1;
    thePacket.kernel = [];
    thePacket.metaData = [];
    
    [paramsFit,fVal,modelResponseStruct{tt}] = ...
        temporalFit.fitResponse(thePacket,...
        'defaultParamsInfo', defaultParamsInfo,...
        'initialParams',initialParams,...
        'vlbParams',vlbParams,...
        'vubParams',vubParams);
    
end
%}


%% create plot if requested

if p.Results.createPlot
    
    for dd = 1:length(p.Results.diagnosis)
        figHandles{dd} = figure();
                
        trialTypes = data{dd};
        plotHandles = [];

        for tt = 1:length(trialTypes)
            
            typeData = trialTypes{tt};
            
            % Get the mean and SEM of the response for each trial type
            yMean = nanmean(typeData);
            ySamples = sum(~isnan(typeData));
            ySEM = nanstd(typeData) ./ sqrt(ySamples);
            
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
%             plot(xVals,yMean+ySEM,':','Color',p.Results.plotColors{tt});
%             plot(xVals,yMean-ySEM,':','Color',p.Results.plotColors{tt});
            %        plot(xVals,modelResponseStruct{tt}.values/100,'-b');            
            
            if tt == length(trialTypes)
                lgd = legend(plotHandles,p.Results.plotLabels,'Location', 'southeast');
                lgd.FontSize = 16;
                title([p.Results.diagnosis{dd}], 'FontSize', 16);
                ylabel('Proportion change pupil area', 'FontSize', 16);
                xlabel('Time [secs]', 'FontSize', 16);
                ylim([-0.15 0.05])
            end
            
        end
    end
end

end