%% Tracking using particle filters
%
% This file is used as an inital test script for implementing 
% "3D Trajectory Prediction of Basketball Shot Using Filtering Techniques
% and Computer Vision" project. The project is a self picked topic for implementation
% in the Appied Estimation course at% the KTH Royal Institute of Technology in 2021.
%
% Authors : 
% Matthew William Lock (mwlock@kth.se)
% Miguel Garcia Naude (magn2@kth.se)

% The code in this file is done in sections with each section being
% indepenent of each other. The idea of tracking an object in a video
% stream is incrementely built in each section in a "chronological" order.

%% Loading a playing a video stream

% Clean environment and load video
close all; clc; clear;
videos = ["shot_1.mp4"];
video = videos(1);

% Get number of frames
v = VideoReader(video); numFrames = 0;
while hasFrame(v)
    readFrame(v);
    numFrames = numFrames + 1;
end

% Load frames
for i = 1:numFrames
    frames(:,:,:,i) = read(v,i);
end

% Loop through and display frames
for i = 1:numFrames
        imshow(frames(:,:,:,i));
        hold on;
        pause(1/v.FrameRate);    
end

%% 1. Convert frames to HSV space

% Clean environment and load video
close all; clc; clear;
videos = ["shot_1.mp4"];
video = videos(1);

% Get number of frames
v = VideoReader(video); numFrames = 0;
while hasFrame(v)
    readFrame(v);
    numFrames = numFrames + 1;
end

% Load frames
for i = 1:numFrames
    frame = read(v,i);
    frames(:,:,:,i) = frame;
    frames_hsv(:,:,:,i) = rgb2hsv(frame);
end
