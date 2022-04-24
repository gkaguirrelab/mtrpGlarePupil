%% figures
% This script produces figures for the paper.

%% temporal structure

names = {'Onset tone' 'Delay (100-500msec)' 'Stimuli' 'Recording' 'Process time'};
startTimes={[0],[100],[1600],[600],[4600]};
endTimes={[100],[600],[2600],[4600],[6450]};
timeline(names, startTimes, endTimes);