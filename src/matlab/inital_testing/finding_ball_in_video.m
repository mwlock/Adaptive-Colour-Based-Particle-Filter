%% Finding ball in video
%
% This file is used as an inital test script for implementing 
% "3D Trajectory Prediction of Basketball Shot Using Filtering Techniques
% and Computer Vision" project. The project is a self picked topic for implementation
% in the Appied Estimation course at% the KTH Royal Institute of Technology in 2021.
%
% Authors : 
% Matthew William Lock (mwlock@kth.se)
% Miguel Garcia Naude (magn2@kth.se)

%% Initalise environment

close all; clc; clear;
videos = ["shot_1.mp4"];
video = videos(1);

%% Load video
v = VideoReader(video);

% Get number of frames
numFrames = 0;
while hasFrame(v)
    readFrame(v);
    numFrames = numFrames + 1;
end

%% Find ball in specific frame

frame = read(v,55);
image = frame;

ball_radius_range = [2 100]; % note that this metric is in pixels        
[centers, radii, metric] = imfindcircles(image,ball_radius_range,'ObjectPolarity','bright');

% Retain the (n) strongest circles according to the metric values.
N = 10;
n= N;
if(size(centers,1) < N)
    n = size(centers,1);
end

centersStrongn = centers(1:n,:); 
radiiStrongn = radii(1:n);
metricStrongn = metric(1:n);

% Show image
imshow(image);

% Draw the (n) strongest circle perimeters over the original image.
viscircles(centersStrongn, radiiStrongn,'EdgeColor','g');

% Plot Center points
axis on;
hold on;
if size(centersStrongn,2)>0
    plot(centersStrongn(:,1),centersStrongn(:,2),'g+','MarkerSize', 15)
end
title('No pre-processing')

%% Load video frames

% Parameters
ball_radius_range = [5 100]; % note that this metric is in pixels  
N = 15;
n= N;

% Loop through frames
for i = 1:numFrames

    frame = read(v,i);
    image = frame;

          
    [centers, radii, metric] = imfindcircles(image,ball_radius_range,'sensitivity',0.99);
    
    % Retain the (n) strongest circles according to the metric values.

    if(size(centers,1) < N)
        n = size(centers,1);
    end
    
    centersStrongn = centers(1:n,:); 
    radiiStrongn = radii(1:n);
    metricStrongn = metric(1:n);
    
    % Show image
    imshow(image,[]);
    drawnow;
    
    % Draw the (n) strongest circle perimeters over the original image.
    viscircles(centersStrongn, radiiStrongn,'EdgeColor','g');
    
    % Plot Center points
    axis on;
    hold on;
    if size(centersStrongn,2)>0
        plot(centersStrongn(:,1),centersStrongn(:,2),'g+','MarkerSize', 15)
    end
    title('No pre-processing');
    drawnow;
    hold off;

end
