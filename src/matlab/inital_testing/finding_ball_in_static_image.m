%% Finding ball in static image
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
images = ["freeshot.jpg","non_ideal_shot_1.png","non_ideal_shot_2.png"];
image = imread(images(1));


%% Find ball in static frame (https://se.mathworks.com/help/images/ref/imfindcircles.html)

ball_radius_range = [15 100]; % note that this metric is in pixels        
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
plot(centersStrongn(:,1),centersStrongn(:,2),'g+','MarkerSize', 15)
title('No pre-processing')

%% Find ball in static frame with pre-processing

k = 6;                % number of clusters used
L = 8;                  % number of iterations
seed = 28;              % seed used for random initialization
scale_factor = 0.6;     % image downscale factor
image_sigma = 1.0;      % image preblurring scale

[ segm, centers ] = kmeans_segm(image, k, L, seed);

Inew = mean_segments(image, segm);
I = overlay_bounds(image, segm);

image = Inew;

ball_radius_range = [15 100]; % note that this metric is in pixels        
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
plot(centersStrongn(:,1),centersStrongn(:,2),'g+','MarkerSize', 15)
title('K-means pre-processing')



