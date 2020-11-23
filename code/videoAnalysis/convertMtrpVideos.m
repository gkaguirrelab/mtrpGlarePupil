function convertMtrpVideos(videoInFileName, videoOutFileName, varargin) 
% Converts mtrp videos to mpeg4 format
%
% Syntax:
%  convertMtrpVideos(videoInFileName, varargin)  
%
% Description:
%  This script re-encodes metropsis pupil videos to mpeg4 format to make
%  them work with MATLAB on multiple OS platforms. 
%
% Inputs:
%   videoInFileName       - String. Path and name to the original video.
%   videoOutFileName      - String. Save path and name
%
% Optional parameters
%   ffmpegPath            - String. path to ffmpeg executable. Use this if
%                           you have not set a path to ffmpeg already and 
%                           cannot call it from the terminal without
%                           including the path.
%
% Outputs:
%   none
%

%% parse input and define variables
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;

% Required
p.addRequired('videoInFileName',@isstr);
p.addRequired('outputPath',@isstr);

% Optional
p.addParameter('ffmpegPath','',@isstr);

% parse
p.parse(videoInFileName, videoOutFileName, varargin{:})

%% Convert the video with ffmpeg
ffmpegCommand = [p.Results.ffmpegPath 'ffmpeg -i' ' ' '"' videoInFileName '"' ' ' '-c:v mpeg4 -q:v 0' ' ' '"' videoOutFileName '"'];
system(ffmpegCommand)

end
