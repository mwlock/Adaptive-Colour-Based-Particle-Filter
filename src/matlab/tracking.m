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
videos = ["shot_1_vid_high_res.mp4"];
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

% Load reference frame
reference_frame = frames(:,:,:,55);
figure;
subplot(1,2,1);
imshow(reference_frame);
title('Reference frame');

% Load binary image
file = matfile('shot_1_vid_high_res_binary_frame_55.mat');
binaryImage = file.binaryImage;
subplot(1,2,2);
imshow(binaryImage);
title('Binary mask');

% Get masked region
reference_frame = reshape(reference_frame, [], 3);
mask_reshaped = reshape(binaryImage, [], 1);
masked_pixels = reference_frame(mask_reshaped,:);
masked_pixels = reshape(masked_pixels, [], 1, 3);

% convert RGB to HSV
masked_pixels_HSV = rgb2hsv(masked_pixels);

% Get Colours histogram
figure;
subplot(1,3,1);
[h_counts,binLocations_h] = imhist(masked_pixels_HSV(:,:,1),8);
title('H');
subplot(1,3,2);
[s_counts,binLocations_s] = imhist(masked_pixels_HSV(:,:,2),8);
title('S');
subplot(1,3,3);
[v_counts,binLocations_v] = imhist(masked_pixels_HSV(:,:,3),4);
title('V');



%% Get target object

% Clean environment and load video
close all; clc; clear;
videos = ["shot_1.mp4"];
video = videos(1);

% Load a frame
v = VideoReader(video); 
frame = read(v,40);
% imshow(frame);

% convert to greyscale image
frame = rgb2gray(frame);

% Adaptive thresholding
T = adaptthresh(frame, 0.95);
BW = imbinarize(frame,T);

imshowpair(frame, BW, 'montage')

%% 1. Convert frames to HSV space

% Clean environment and load video
close all; clc; clear;
videos = ["shot_1_vid_high_res.mp4"];
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
