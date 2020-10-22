

p.Results.channelIDs = [0];
p.Results.frequencyInHz = 2000;
p.Results.recordingDurationSecs = 10/1000;
p.Results.verbose = false;


%% Instantiate a LabJack object
labjackOBJ = LabJackU6('verbosity', double(false));

%% Listen to the LabJack in 10 msec increments until a TTL pulse is seen

% Configure analog input sampling
labjackOBJ.configureAnalogDataStream(p.Results.channelIDs, p.Results.frequencyInHz);

waiting = true;
waitCount = 1;
while waiting
    
    labjackOBJ.startDataStreamingForSpecifiedDuration(p.Results.recordingDurationSecs);

    
    if max(labjackOBJ.data) > 1
        waiting = false;
    else
        waitCount = waitCount+1;
    end
end

foo=1;


% Close-up shop
%    labjackOBJ.shutdown();
 
    

recordCommandStem = 'ffmpeg -hide_banner -video_size 1280x720 -framerate 60.000240 -f avfoundation -i "1" -t 10 "foo.mp4"';
