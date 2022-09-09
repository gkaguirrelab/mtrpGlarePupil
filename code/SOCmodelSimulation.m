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
c = [0.9276, 0.9928, 0.9941, 0.9472];
n = [0.1814, 0.1285, 0.1195, 0.1152];
sd = [0.9308, 1.0738, 1.4671, 2.1242];

% Specify save names and struct containers
responseSaveNames = {'response1', 'response2', 'response3', 'response4'};
cacheSaveNames = {'cache1', 'cache2', 'cache3', 'cache4'};
responses = [];
caches = [];
%% Run models 
for ii = 1:length(responseSaveNames)
    % Initiate caches
    caches.(cacheSaveNames{ii}) = [];
    % Run SOC Model
    [responses.(responseSaveNames{ii}),caches.(cacheSaveNames{ii})] = socmodel(stimulus,150,[0 1 0.5],1.2,{37.5 -1 1 8 2 .01 2 0}, ...
                                                                               1,.5,sd(ii),1/sd(ii),n(ii),c(ii),caches.(cacheSaveNames{ii}));   
%     % Get stimulus2 stage from cache. This contains gabor fitted images. We
%     % get these and calculate the complex-cell energy here again, because
%     % we will deviate here for the SOV model
%     gaborFittedStim = caches.(cacheSaveNames{ii}).stimulus2;
%     complexCellStimulus = sqrt(blob(gaborFittedStim.^2,2,2));
%     
%     % Get resolution
%     nrGaborOrientations = 8;
%     res = sqrt(size(stimulus,2)/nrGaborOrientations);
%     [~,xx,yy] = makegaussian2d(res,2,2,2,2);
%     
%     % Given a set of parameters this outputs a 2D gaussian
%     gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));
% 
%     % SOV model
%     modelfun = @(pp,dd) pp(4) * var(reshape(reshape(dd,[],res*res) * gaufun1(pp),[],nrOrientations),[],2).^pp(5);
%     modelfit = modelfun(params, stimulus);
end                

% visualize the results
figure; setfigurepos([100 100 700 300]);
subplotCounter = 0;
stimNames = {'glare', 'halo', 'stripe', 'uniform'};
fprintf('Plotting... \n')
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