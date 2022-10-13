clear all; clc

% Load Hermes and Aguirre stimulus
load('C:\Users\ozenc\Documents\MATLAB\projects\Paper_Hermes_2019_eLife\data\stimuli\task-soc_stimuli.mat', 'stimuli');
hermesStimulus = stimuli; 

% Load Aguirre stimulus
currentScriptPath = strsplit(mfilename('fullpath'),filesep);
currentScriptPath = currentScriptPath(1:end-1);
stimulusPath = fullfile(currentScriptPath{:}, 'stimulus');
stimulusFiles = dir(stimulusPath);
stimulusFiles(1:2) = [];
aguirreStimulus = [];
for ii = 1:length(stimulusFiles) 
    ss = importdata(fullfile(stimulusFiles(ii).folder, stimulusFiles(ii).name));
    aguirreStimulus(:,:,ii) = ss;
end

% Combine stimuli into a single cell 
stimuli = {}; 
stimuli{1} = single(hermesStimulus);
stimuli{2} = single(aguirreStimulus);

% Set parameters gamma and soc
gaussianStd = 40;
powerLawExponent = 0.1814;
secondOrderStrength = 0.9276;

% Helper function
socfun = @(dd,wts,c) bsxfun(@minus,dd,c*(dd*wts)).^2 * wts;

% Loop through both stimulus set
for ii = 1:length(stimuli)
    % If Hermes stimuli just downsample to 240x240
    if isequal(ii,1)
        stimuli{ii} = processmulti(@imresize,stimuli{ii},[240 240],'cubic');
    % If Aguirre stimuli, first downsample to 240x240, crop to make the 
    % visual field similar to Hermes' stimuli, and resample to 240x240 
    elseif isequal(ii,2)
        temp = zeros(240,240,size(stimuli{ii},3),'single');
        for im = 1:size(stimuli{ii},3)
            procIm = imresize(imcrop(imresize(single(stimuli{ii}(:,:,im)),[240 240],'cubic'), [10, 10, 240 - 2*10, 240 - 2*10]), [240 240]);
            temp(:,:,im) = procIm;
        end
        stimuli{ii} = temp;
    end

    % Rescale to 0-254, values 0-1, substract background lum 0.5
    if isequal(ii, 1)
        stimuli{ii}(stimuli{ii} < 0) = 0;
        stimuli{ii}(stimuli{ii} > 254) = 254;
        stimuli{ii} = stimuli{ii}/254 - 0.5;
    else
        rng = [min(aguirreStimulus(:)) max(aguirreStimulus(:)) 0.5];
        stimuli{ii} = normalizerange(stimuli{ii},0,1,rng(1),rng(2)) - normalizerange(rng(3),0,1,rng(1),rng(2));
    end
end
clear temp;

% Apply Gabor filters to the stimuli.  filters occur at different positions,
% orientations, and phases.  there are several parameters that govern the
% design of the filters:
filt_prop.cycles = 60*(270/240);    %   the number of cycles per image is 60*(270/240)
filt_prop.bandwidth = -1;           %   the spatial frequency bandwidth of the filters is 1 octave
filt_prop.spacings=1;               %   the separation of adjacent filters is 1 std dev of the Gaussian envelopes
                                    %     (this results in a 135 x 135 grid of positions)
filt_prop.orientations=8;           %   filters occur at 8 orientations
filt_prop.phases=2;                 %   filters occur at 2 phases (between 0 and pi)
filt_prop.thres=0.01;               %   the Gaussian envelopes are thresholded at .01
filt_prop.scaling=2;                %   filters are scaled to have an equivalent Michelson contrast of 1
filt_prop.mode=0;                   %   the dot-product between each filter and each stimulus is computed

% Initialize gamma output 
gamma = {};
bold = {};

for ii = 1:2
    % Pad
    stimuli{ii} = placematrix(zeros(270,270,size(stimuli{ii},3),'single'),stimuli{ii});

    % after this step, stimulus is images x phases*orientations*positions.
    stimuli{ii} = applymultiscalegaborfilters(reshape(stimuli{ii},270*270,[])', ...
      filt_prop.cycles,filt_prop.bandwidth,filt_prop.spacings,filt_prop.orientations,...
      filt_prop.phases,filt_prop.thres,filt_prop.scaling,filt_prop.mode);

    % Complex cell energy. 
    stimuli{ii} = sqrt(blob(stimuli{ii}.^2,2,2));

    % Get resolution create a gaussian. Gaussian values are hard coded.
    % Peak is in the middle of the image, the size was preset to a large 2D. 
    res = sqrt(size(stimuli{ii},2)/filt_prop.orientations);
    pRF = vflatten(makegaussian2d(res,[],[],gaussianStd,gaussianStd)/(2*pi*gaussianStd^2));
    gamma{ii} = 1 * var(reshape(reshape(stimuli{ii},[],res*res) * pRF,[],filt_prop.orientations),[],2).^powerLawExponent;
    
    % Population term in divisive normalization and repeat it for each
    % orientation
    stimulusPOP = blob(stimuli{ii},2,8)/8;
    stimulusPOP = upsamplematrix(stimulusPOP,8,2,[],'nearest');

    % apply divisive normalization to the complex-cell outputs.  there are two parameters
    % that influence this operation: an exponent term (r) and a semi-saturation term (s).
    % the parameter values specified here were determined through a separate fitting
    % procedure (see paper for details).  for the purposes of this script, we will
    % simply hard-code the parameter values here and not worry about attempting to fit
    % the parameters.
    r = 1;
    s = 0.5;
    stimuli{ii} = stimuli{ii}.^r ./ (s.^r + stimulusPOP.^r);
    clear stimulusPOP;

    % sum across orientation.  after this step, stimulus is images x positions.
    imEnergyMean = blob(stimuli{ii},2,8);
    
    % Calculate BOLD (SCO model) 
    
    bold{ii} = 1*socfun(imEnergyMean, pRF, restrictrange(secondOrderStrength,0,1)).^powerLawExponent;
end

figure 
subplot(2,1,1)
bar(bold{1})
ylim([0 1])
title('Broadband - BOLD')
subplot(2,1,2)
bar(gamma{1})
ylim([0 1])
title('Narrowband')

figure
subplot(1,2,1)
bar(bold{2})
set(gca,'XTickLabel',{'glare', 'halo', 'stripe', 'uniform'});
xtickangle(45)
ylim([0 0.6])
title('Broadband - BOLD')
subplot(1,2,2)
bar(gamma{2})
set(gca,'XTickLabel',{'glare', 'halo', 'stripe', 'uniform'});
xtickangle(45)
ylim([0 0.6])
title('Narrowband')