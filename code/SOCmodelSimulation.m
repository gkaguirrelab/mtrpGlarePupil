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

% initialize cache for model runs
cache1 = [];
cache2 = [];
cache3 = [];
cache4 = [];

% Run the model with typical V1 response values
c = 0.9276;
n = 0.1814;
sd = 0.9308;
[response,cache1] = socmodel(stimulus,150,[0 1 0.5],1.2,{37.5 -1 1 8 2 .01 2 0}, ...
                            1,.5,sd,1/sd,n,c,cache1);

% Run it for typical V2 vals
c = 0.9928;
n = 0.1285;
sd = 1.0738;
[response2,cache2] = socmodel(stimulus,150,[0 1 0.5],1.2,{37.5 -1 1 8 2 .01 2 0}, ...
                            1,.5,sd,1/sd,n,c,cache2);

% Run it for typical V3 vals
c =  0.9941;
n = 0.1195;
sd = 1.4671;
[response3,cache3] = socmodel(stimulus,150,[0 1 0.5],1.2,{37.5 -1 1 8 2 .01 2 0}, ...
                            1,.5,sd,1/sd,n,c,cache3);   
                        
% Run it for typical hV4 vals
c = 0.9472;
n = 0.1152;
sd = 2.1242;
[response4,cache4] = socmodel(stimulus,150,[0 1 0.5],1.2,{37.5 -1 1 8 2 .01 2 0}, ...
                            1,.5,sd,1/sd,n,c,cache4);                         
                        

% visualize the results
figure; setfigurepos([100 100 700 300]);
subplotCounter = 0;
stimNames = {'glare', 'halo', 'stripe', 'uniform'};
fprintf('Plotting... \n')
for p=1:4
  subplot(4,5,subplotCounter + 1);
  imagesc(stimulus(:,:,p),[0 1]); axis image tight; colormap(gray); colorbar; title(stimNames{p});
  subplot(4,5,subplotCounter + 2);
  imagesc(response(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V1');
  subplot(4,5,subplotCounter + 3);
  imagesc(response2(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V2');
  subplot(4,5,subplotCounter + 4);
  imagesc(response3(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('V3');
  subplot(4,5,subplotCounter + 5);
  imagesc(response4(:,:,p),[0 1.5]); axis image tight; colormap(gray); colorbar; title('hv4');
  subplotCounter = subplotCounter + 5;
end