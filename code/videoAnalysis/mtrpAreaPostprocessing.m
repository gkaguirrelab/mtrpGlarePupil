function mtrpAreaPostprocessing(pupilFileName, outputFileName, varargin)
% Threshold pupil area with RMSE values
%
% Syntax:
%  mtrpAreaPostprocessing(pupilFileName, outputFileName, varargin)
%
% Description:
%   This script loads pupil.mat files, extracts the area vector and
%   thresholds it using the RMSE values in order to get rid of blinks in
%   the values.
%
% Required inputs:
%   pupilFileName         - String. Path to the pupil.mat file
%   outputFileName        - String. Output file path and name
%   
% Optional inputs:
%   threshold             - Number. Threshold to use to remove blinks.
%                           Timepoints that have larger RMSE than this
%                           threshold will be removed from the area vector.
%                           Default = 0.9
%                           
% Outputs:
%   none
%

%% parse input and define variables
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;

% Required
p.addRequired('pupilFileName',@isstr);
p.addRequired('outputFileName',@isstr);

% Optional
p.addParameter('threshold',0.9,@isnumeric);

% parse
p.parse(pupilFileName, outputFileName, varargin{:})

%% Process area 

% Load the pupil file 
load(pupilFileName)

% Extract the area and the RMSE vectors
area = pupilData.initial.ellipses.values(:,3);
rmse = pupilData.initial.ellipses.RMSE;

% Threshold RMSE values and get the indices of values larger than the
% threshold
largeRmseIndices = find(rmse >= p.Results.threshold);

% Remove these from the area vector
area(largeRmseIndices) = NaN;

% Calculate percentage change 
areaPercentageChange = diff(area)./area(1:end-1);

% Plot area
figure('visible', 'off');
plot(area)
plotSavePath = strrep(outputFileName,'_area','_areaPlot');
plotSavePath = strrep(plotSavePath,'.mat','.png');
saveas(gcf, plotSavePath)

% Save the new area and percentage change file 
save(outputFileName, 'area', 'areaPercentageChange')
