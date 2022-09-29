clear all
clc

% We need Kendrick's knkutils toolbox. Check if it exists. Pull if it
% doesn't.
currentScriptPath = matlab.desktop.editor.getActiveFilename;
matlabPath = currentScriptPath(1:end-50);
knkutilsFolder = fullfile(matlabPath, 'toolboxes', 'knkutils');
if ~isfolder(knkutilsFolder)
    fprintf('Pulling knkutils \n')
    system(['git clone https://github.com/cvnlab/knkutils.git ' knkutilsFolder]);
end
addpath(genpath(knkutilsFolder));

% Load glare, halo, stripe, uniform stim and concatanate across time
stimulusPath = fullfile(matlabPath, 'projects', 'mtrpGlarePupil', 'code', 'stimulus');
stimulusFiles = dir(stimulusPath);
stimulusFiles(1:2) = [];
stimulus = [];
for ii = 1:length(stimulusFiles) 
    ss = importdata(fullfile(stimulusFiles(ii).folder, stimulusFiles(ii).name));
    stimulus(:,:,ii) = ss;
end

% Specify typical values for V1, V2, V3, and hv4 areas
c = [0.9276, 0.9928, 0.9941, 0.9472]; % Strength of second order contrast 
n = [0.1814, 0.1285, 0.1195, 0.1152]; % Exponent of power-law linearity
sd = [0.9308, 1.0738, 1.4671, 2.1242]; % SD of 2D gaussian applied to image
spacing = [1/sd(1), 1/sd(2), 1/sd(3), 1/sd(4)]; % Number of sd that separate 2D gaussians 

% Specify save names and struct containers
responseSaveNames = {'response1', 'response2', 'response3', 'response4'};
cacheSaveNames = {'cache1', 'cache2', 'cache3', 'cache4'};
responses = [];
caches = [];
gammaResponsesMean = [];
%% Run models 
for ii = 1:length(responseSaveNames)
    % Initiate caches
    caches.(cacheSaveNames{ii}) = [];
    
    % Run SOC Model which gives us the BOLD.
    [responses.(responseSaveNames{ii}),caches.(cacheSaveNames{ii})] = socmodel(stimulus,150,[0 1 0.5],1.2,{37.5 -1 1 8 2 .01 2 0}, ...
                                                                               1,.5,sd(ii),spacing(ii),n(ii),c(ii),caches.(cacheSaveNames{ii})); 
                                                                           
    % RUN THE SOV model now.                                                                      
    % Get stimulus2 stage from cache. This contains gabor fitted images. We
    % get these and calculate the complex-cell energy here again seprately, 
    % because we will deviate from this point for the SOV model calculation.
    gaborFittedStim = caches.(cacheSaveNames{ii}).stimulus2;
    complexCellStimulus = sqrt(blob(gaborFittedStim.^2,2,2));
    
    % Set number of gabor orientation and set resolution.
    nrGaborOrientations = 8;
    res = sqrt(size(complexCellStimulus,2)/nrGaborOrientations);

    % Make gaussian filters, make full. Not sure if this is correct. 
    [d,d,d,d,filters] = applymultiscalegaussianfilters(randn(1,res^2),[sd(ii)],[1/sd(ii)],.01,2);
    filters = full(filters);
%     nz = filters ~= 0;
    
    % Do the model calculation.
    % The original model commented out below. 
    % modelfun = @(pp,dd) pp(4) * var(reshape(reshape(dd,[],res*res) * gaufun1(pp),[],nrOrientations),[],2).^pp(5);
  
    % I hardcode the prf gain or pp(4) to 1. 
    % pp(5) is power law power, so pass the n set. 
    modelfit = 1 * var(reshape(reshape(complexCellStimulus,[],res*res) * filters,[],nrGaborOrientations),[],2).^n(ii);

    % Separate the matrix into gamma values
    combined = [];
    combined = [combined modelfit(1:4:length(modelfit)) modelfit(2:4:length(modelfit)) modelfit(3:4:length(modelfit)) modelfit(4:4:length(modelfit))];
    
    % Get the mean across PRFs
    gammaResponsesMean(:,ii) = mean(combined)';
%     F = sqrt(size(gausGamma',2));
%     response = reshape(gausGamma',F,F,[]);
end                

% visualize the results
f = figure; setfigurepos([100 100 700 300]);
subplotCounter = 0;
stimNames = {'glare', 'halo', 'stripe', 'uniform'};
fprintf('Plotting SOC model - BOLD.... \n')
for p=1:4
  subplot(4,5,subplotCounter + 1);
  imagesc(stimulus(:,:,p),[0 1]); axis image tight; colormap(gray); colorbar; title(stimNames{p});
  subplot(4,5,subplotCounter + 2);
  imagesc(responses.response1(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V1');
  subplot(4,5,subplotCounter + 3);
  imagesc(responses.response2(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V2');
  subplot(4,5,subplotCounter + 4);
  imagesc(responses.response3(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V3');
  subplot(4,5,subplotCounter + 5);
  imagesc(responses.response4(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('hv4');
  subplotCounter = subplotCounter + 5;
end
t = table(['stim1'; 'stim2'; 'stim3'; 'stim4'], gammaResponsesMean(:,1), gammaResponsesMean(:,2), gammaResponsesMean(:,3), gammaResponsesMean(:,4));
t.Properties.VariableNames = ["Stimulus SOV - gamma", "V1 mean", "V2 mean", "V3 mean", "hV4 mean"];
t