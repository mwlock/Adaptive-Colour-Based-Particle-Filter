%% Background subtraction to find moving objects
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
reference_frame = read(v,1);

%% Load all frames
for i = 1:numFrames
    frames(:,:,:,i) = read(v,i);
end

%% Compute subtrackted frame
subtrackted_frames = zeros([size(reference_frame) numFrames]);

for i = 1:numFrames
    subtrackted_frames(:,:,:,i) = abs(im2double(frames(:,:,:,i)) - im2double(reference_frame));
    reference_frame = frames(:,:,:,i);
end

%% Pure Background scubtraction

% Loop through frames
for i = 1:numFrames

        imshow(subtrackted_frames(:,:,:,i));

        drawnow;
        axis on;
        hold on;
        pause(1/30);
    
end

%% Binary Background scubtraction

% Loop through frames
for i = 1:numFrames
        
        mask = im2bw(subtrackted_frames(:,:,:,i),0.3);
        rgbImage = frames(:,:,:,i);
        maskedRgbImage = bsxfun(@times, rgbImage, cast(mask,class(rgbImage)));
        imshow(maskedRgbImage);

        drawnow;
        axis on;
        hold on;
        pause(1/30);
    
end

%% TRying to find circles with background subtraction

ball_radius_range = [2 100]; % note that this metric is in pixels  
N = 15;
n= N;

% Loop through frames
for i = 1:numFrames
    
    image = subtrackted_frames(:,:,:,i);
          
    [centers, radii, metric] = imfindcircles(image,ball_radius_range,'sensitivity',0.2);
    
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
    drawnow;

end



