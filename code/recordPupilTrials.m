% recordPupilTrials
%
% This script records pupil videos via an ffmpeg command in response to "t"
% keypress events transmitted from a Metropsis experiment (via a TTL->USB
% converter).
%
% To determine which device to record from, issue thecommand
%{
	ffmpeg -f avfoundation -list_devices true -i ""
%}
% in the terminal. Identify which device number we want, and place that in
% the quotes after the -i in the command stem below.


% Set the parameters of the experiment
protocolParams.trialDurationSecs = 4;
protocolParams.recordCommandStem = 'ffmpeg -hide_banner -video_size 640x480 -framerate 60.500094 -f avfoundation -i "1" -t trialDurationSecs "videoFileOut.mp4"';
protocolParams.experimentName = 'pupilGlare_01';

% Query the user for the subject ID, date, and session
prompt = {'Subject ID:','yyyy-mm-dd','Session name:'};
dlgtitle = 'Subject and session info';
dims = [1 35];
definput = {'GLAR_01',datestr(now, 'yyyy-mm-dd'),'session_1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

% Store the responses in the protocol params
protocolParams.observerID = answer{1};
protocolParams.todayDate = answer{2};
protocolParams.sessionName = answer{3};

% Set the save location for the data
protocolParams.dropBoxBaseDir = getpref('mtrpGlarePupil','dropboxBaseDir');

% Create a directory for the data (if it does not already exist)
protocolParams.dataOutDir = fullfile(...
    protocolParams.dropBoxBaseDir,...
    'MTRP_data',...
    protocolParams.experimentName,...
    protocolParams.observerID,...
    protocolParams.todayDate,...
    protocolParams.sessionName);

if ~exist(protocolParams.dataOutDir,'dir')
    mkdir(protocolParams.dataOutDir)
end

% Put the instructions on the console
fprintf(['Preparing to record video trials for ' protocolParams.observerID '.\n']);
fprintf('Adjust camera.\n');

% Enter the camera adjustment loop
figHandle = figure();
ButtonHandle = uicontrol('Style', 'PushButton', ...
    'String', 'Stop loop', ...
    'Callback', 'delete(gcbf)');

stillRecording = true;
while stillRecording
    
    % Record a video snippet
    tmpVid = [tempname '.mp4'];
    vidCommand = protocolParams.recordCommandStem;
    vidCommand = strrep(vidCommand,'trialDurationSecs','0.33');
    vidCommand = strrep(vidCommand,'videoFileOut.mp4',tmpVid);
    [status,cmdout] = system(vidCommand);
    
    % Extact the mid time point of that snippet
    tmpIm = [tempname '.jpg'];
    extractCommand = ['ffmpeg -ss 00:00:00.17 -i ' tmpVid ' -vframes 1 -q:v 2 ' tmpIm];
    [status,cmdout] = system(extractCommand);
    
    % Load the image
    imagesc(imread(tmpIm));
    axis off
    hold on
    drawnow
    
    % Delete the temp files
    delete(tmpVid);
    delete(tmpIm);
    
    % Check if we have hit the stop button
    if ~ishandle(ButtonHandle)
        stillRecording = false;
    end
    
end

% Start the experiment
protocolParams.experimentStartTime = datestr(now);
fprintf('Starting the experiment. Press "q" when done.\n');
trialIdx = 1;
trialTimes = [];
stillRecording = true;

while stillRecording
    
    % Announce we are ready for a trial
    fprintf('Trial %d. Waiting for t...',trialIdx);
    
    % Assemble to video recording command
    vidOutFile = fullfile(protocolParams.dataOutDir,sprintf('trial_%02d.mp4',trialIdx));
    vidCommand = protocolParams.recordCommandStem;
    vidCommand = strrep(vidCommand,'trialDurationSecs',num2str(protocolParams.trialDurationSecs));
    vidCommand = strrep(vidCommand,'videoFileOut.mp4',vidOutFile);
    
    % Wait for a keypress
    keyPress = getkey(1, 'non-ascii');

    % Act upon the keyPress
    switch keyPress
        case 'q'
            stillRecording = false;
            fprintf('exiting.\n');
        case 't'
            trialTimes(trialIdx) = now;
            fprintf('recording...');
            [status,cmdout] = system(vidCommand);
            fprintf('done.\n');
            trialIdx = trialIdx+1;
    end
    
end

% Save the protocol params
protocolParams.trialTimes = trialTimes;
paramFileOut = fullfile(protocolParams.dataOutDir,'protocolParams.mat');
save(paramFileOut,'protocolParams');
