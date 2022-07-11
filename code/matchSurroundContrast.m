%% Calculate integrated contrast for Glare experiment
%
% The Glare experiment makes use of radial stimuli that have a center and a
% surround, placed on a background. The surround can contain a gradient
% that changes in contrast linearly from the inner ring of the surround
% annulus to the outer ring. For any given contrast value selected for the
% center and background, this routine provides the average contrast across
% the surround.
%
% Contrast is specified in the range of -1 to 1, with the half-on state of
% the display at zero, and the max and min primary settings corresponding
% to -1 and 1, respectively. This odd contrast convention is what is used
% by the Metropsis / Psykinematix system.
%
% The values defined for the inner and outer radii of the surround for the
% three stimuli create equal integrated contrast for all three stimulus
% types.
%


%% Housekeeping
%close all

%% Define properties of the stimulus

% The radius (in degrees visual angle) of the the surround
radius_g = 9;

% The radius (in degrees visual angle) of the center
radius_c = 3.4;

% Contrast of the background
contrast_background = 0.0;

% Contrast of the center
contrast_center = 1.0;

% The names of the stimuli
stimulus_labels = {'glare','halo','uniform','stripe'};

% Contrast of the inner edge of the surround for the glare, halo, and
% uniform
contrast_inner = [1, -0.5, 0.137, 0.137 ];

% Contrast at the outer edge of the surround for glare, halo, and uniform
contrast_outer = [-0.5, 0.607, 0.137, 0.607 ];

% Number of sinusoidal cycles to include in the "stripe" stimulus
numCycles = 3;

% Find the phase shift that places the start of the cosine at the same
% point that the uniform occupies adjacent to the pedastal.
phaseShift = pi;

%% Create visualization of the stimuli
% Loop through the sets

for ss = 1:length(contrast_inner)
    
    % Define the X, Y positions of the stimulus image over a domain of -10 to
    % 10 degrees of visual angle
    [X, Y] = meshgrid(-10 : 0.01 : 10, -10 : 0.01 : 10);
    
    % R is the distance of each point in the grid from the image center
    R = sqrt(X .^ 2 + Y .^ 2);
    
    % Initialize the stimulus image with the background
    stimImage = zeros(size(X)) + contrast_background;

    % Render the center
    stimImage(R < radius_c) = contrast_center;
    
    % Render the stimus surround
    if ss <=3
        gradient = @(r) (contrast_inner(ss) - (contrast_inner(ss) - contrast_outer(ss)) .* ((r - radius_c) ./ (radius_g - radius_c)));
    else
        gradient = @(r) contrast_inner(ss) + contrast_outer(ss) .* sin(phaseShift + ((r - radius_c) ./ (radius_g - radius_c)).*2.*pi.*numCycles);
    end

    for ii = 1:size(stimImage, 1)
        for jj = 1:size(stimImage, 2)
            if R(ii, jj) > radius_c && R(ii, jj) <= radius_g
                stimImage(ii, jj) = gradient(R(ii, jj));
            end
        end
    end
    
    % Throw a -1 in the corner to render the full image range.
    stimImage(1,1)=-1;
    
    % Display the stimulus
    visStim(stimImage);
    
    % Add a plot title
    title(stimulus_labels{ss})
    
    % Display a cross-section through the stimulus
    subplot(1, 2, 2);
    xRange = -10 : 0.01 : 10;
    plot(xRange, stimImage(round(length(xRange)/2),:), 'k', 'LineWidth', 2);
    ylim([-1 1]);
    ylabel('Contrast');
    xlabel('Position [degs]');
    axis square
    
    
    %% Calculate the integrated contrast of the surround
    % This is a unit test of the polar integration
    %{
	circle = @(r) ones(size(r));
	polarFun = @(theta, r) circle(r) .* r;
	q = integral2(polarFun, -pi, pi, 0, 1);
    %}
    
    % Define the polar functions
    polarFun  = @(theta, r) gradient(r) .* r ;
    radiusFun = @(theta, r) r;
    
    % Integrate
    val = integral2(polarFun, -pi, pi, radius_c, radius_g) ./ integral2(radiusFun, -pi, pi, radius_c, radius_g);
    
    % Report the contrast of the surround
    fprintf(['The integrated contrast of the ' stimulus_labels{ss} ' surround is %.3f \n'], val);
    
end




%% Helper function
function visStim(stimulus)
figure();
subplot(1, 2, 1);
stimulus = stimulus - min(stimulus(:));
stimulus = stimulus ./ max(stimulus(:));
imshow(stimulus);
end

