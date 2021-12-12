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

% The code in this file is done in sections with each section meant 
% to run one after the next. The idea of tracking an object in a video
% stream is incrementely built in each section in a "chronological" order.

%% Get colour distriubtion of the target

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

% Get image parameters
image_height = size(reference_frame,1);
image_width = size(reference_frame,2);

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
h_counts = h_counts';
imhist(masked_pixels_HSV(:,:,1),8);
title('H');
subplot(1,3,2);
[s_counts,binLocations_s] = imhist(masked_pixels_HSV(:,:,2),8);
s_counts = s_counts';
imhist(masked_pixels_HSV(:,:,2),8);
title('S');
subplot(1,3,3);
[v_counts,binLocations_v] = imhist(masked_pixels_HSV(:,:,3),4);
v_counts = v_counts';
v_counts(8)=0;
imhist(masked_pixels_HSV(:,:,3),4);
title('V');

% Generate target histogram
target_histogram = [h_counts;s_counts;v_counts];

%% Initiate PF

% Generate initial particles set
M = 200;                % number of particles
S = zeros(3,M);         % set of particles   

S(1,:) = rand(1,M)*(image_height-1)+1;       % y     
S(2,:) = rand(1,M)*(image_width-1)+1;        % x

%% Predict (with rectangles drawn around particles)

% Look at how particles move over time
R = diag([10 10]);                                  % process noise 

% Size of rectangles to draw
rect_width = 30;
rect_height = 30; 

% Predict 100 times (merely an example)
for i = 1:100
    S = predict_noise(S,R,M);
    plot(S(2,:),S(1,:),'.');
    xlim([0 image_width]);
    ylim([0 image_height]);

    % Draw rectangles
    for r =1:M
        xLeft = S(2,r) - rect_width/2;
        yBottom = S(1,r) - rect_height/2;
        rectangle('Position',[xLeft,yBottom,rect_width,rect_height],'EdgeColor','b','LineWidth',1);
    end

    pause(1/60);    
end

%% Calculate particle weights (for single frame)

% Loop over all histograms
for hist_index = 1:M
    
    % Create two logical matrices
    logical_image_1 = false(image_height,image_width);
    logical_image_2 = false(image_height,image_width);
    
    % Create regions of interest
    region_x = round(S(2,hist_index)-rect_width/2)+(0:rect_width-1);
    region_y = round(S(1,hist_index)-rect_height/2)+(0:rect_height-1);

    % Check if regions of interest are out of bounds and correct
    region_x = region_x ((image_width+1 > region_x) & (region_x > 0));
    region_y = region_y ((image_height+1 > region_y) & (region_y > 0));


    % Bound regions of interest
    logical_image_1(:,region_x) = true;
    logical_image_2(region_y,:) = true;
    logical_image = and(logical_image_1,logical_image_2);

end
