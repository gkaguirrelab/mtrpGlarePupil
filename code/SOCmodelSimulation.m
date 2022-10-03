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
    %% SOC model, Kay - BOLD
    % Initiate caches
    caches.(cacheSaveNames{ii}) = [];
    
    % Run SOC Model which gives us the BOLD.
    [responses.BOLD.(responseSaveNames{ii}),caches.(cacheSaveNames{ii})] = socmodel(stimulus,150,[0 1 0.5],1.2,{37.5 -1 1 8 2 .01 2 0}, ...
                                                                               1,.5,sd(ii),spacing(ii),n(ii),c(ii),caches.(cacheSaveNames{ii})); 
                                                                           
    %% SOV model, Hermes - Gamma                                                          
    % Get stimulus2 stage from SOC image cache. This contains gabor fitted
    % images. We get these and re-calculate complex-cell energy model
    % again. This is where the SOV model deviates from the SOC model
    gaborFittedStim = caches.(cacheSaveNames{ii}).stimulus2;
    complexCellStimulus = sqrt(blob(gaborFittedStim.^2,2,2));
    
    % Set number of gabor orientation and get the resolution.
    nrGaborOrientations = 8;
    res = sqrt(size(complexCellStimulus,2)/nrGaborOrientations);
    
    % Commented out below is the original model for Gamma calculation. 
    % The formula is response = gain * var(gaborFilteredStim * 2DpRF).^powerLawPower
    % modelfun = @(pp,dd) pp(4) * var(reshape(reshape(dd,[],res*res) * gaufun1(pp),[],nrOrientations),[],2).^pp(5);
    
    % Create a set of gaussian filters, make full. Check this with Geoff.
    [d,d,d,d,filters] = applymultiscalegaussianfilters(randn(1,res^2),[sd(ii)],[1/sd(ii)],.01,2);
    filters = full(filters);
%     nz = filters ~= 0;
  
    % Loop through the pRF filters and calculate the model for each pRF.
    % We hardcode the prf gain or to 1. And pass the same power-law values 
    % that we used for the SOC model. We add the results into a one big
    % vector.
    modelfit = [];
    for tt = 1:length(filters)
        fit = 1 * var(reshape(reshape(complexCellStimulus,[],res*res) * filters(:,tt),[],nrGaborOrientations),[],2).^n(ii);
        modelfit = [modelfit; fit];
    end
    
    % Separate the vector into gamma values for each stimulus.
    combined = [];
    combined = [combined modelfit(1:4:length(modelfit)) modelfit(2:4:length(modelfit)) modelfit(3:4:length(modelfit)) modelfit(4:4:length(modelfit))];
    
%     % Get the mean across PRFs
%     gammaResponsesMean(:,ii) = mean(combined)';

    % Reshape the vectors into images and add it into the gamma
    F = sqrt(size(combined',2));
    responses.gamma.(responseSaveNames{ii}) = reshape(combined,F,F,[]);
end                

% visualize the results for BOLD
f = figure; setfigurepos([100 100 700 300]);
subplotCounter = 0;
stimNames = {'glare', 'halo', 'stripe', 'uniform'};
fprintf('Plotting SOC model - BOLD.... \n')
for p=1:4
  subplot(4,5,subplotCounter + 1);
  imagesc(stimulus(:,:,p),[0 1]); axis image tight; colormap(gray); colorbar; title(stimNames{p});
  subplot(4,5,subplotCounter + 2);
  imagesc(responses.BOLD.response1(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V1');
  subplot(4,5,subplotCounter + 3);
  imagesc(responses.BOLD.response2(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V2');
  subplot(4,5,subplotCounter + 4);
  imagesc(responses.BOLD.response3(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V3');
  subplot(4,5,subplotCounter + 5);
  imagesc(responses.BOLD.response4(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('hv4');
  subplotCounter = subplotCounter + 5;
end
sgtitle('SOC model - BOLD, Kay') 

% Do the same for gamma
% visualize the results for BOLD
f2 = figure; setfigurepos([100 100 700 300]);
subplotCounter = 0;
stimNames = {'glare', 'halo', 'stripe', 'uniform'};
fprintf('Plotting SOC model - gamma.... \n')
for p=1:4
  subplot(4,5,subplotCounter + 1);
  imagesc(stimulus(:,:,p),[0 1]); axis image tight; colormap(gray); colorbar; title(stimNames{p});
  subplot(4,5,subplotCounter + 2);
  imagesc(responses.gamma.response1(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V1');
  subplot(4,5,subplotCounter + 3);
  imagesc(responses.gamma.response2(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V2');
  subplot(4,5,subplotCounter + 4);
  imagesc(responses.gamma.response3(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V3');
  subplot(4,5,subplotCounter + 5);
  imagesc(responses.gamma.response4(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('hv4');
  subplotCounter = subplotCounter + 5;
end
sgtitle('SOV model - gamma, Hermes') 