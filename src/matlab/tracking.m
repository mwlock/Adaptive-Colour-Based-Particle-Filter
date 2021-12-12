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

%% Get colour distribution of the target

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
target_histogram = target_histogram / sum(target_histogram(1,:),2);

% Reshape reference back to original
reference_frame = reshape(reference_frame,image_height,image_width,3);

%% Initiate PF

% Generate initial particles set
M = 200;                % number of particles
S = zeros(3,M);         % set of particles   

S(1,:) = rand(1,M)*(image_height-1)+1;       % y     
S(2,:) = rand(1,M)*(image_width-1)+1;        % x

%% Rectangle dimensions

% Size of rectangles to draw
rect_width = 32;
rect_height = 32; 

%% Predict (with rectangles drawn around particles)

% Look at how particles move over time
R = diag([10 10]);                                  % process noise  

% Predict 100 times (merely an example - shows moving particles)

for i = 1:100
    
    hold off
    imshow(reference_frame);
    hold on;
    
    % Predict motion of particles
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

%% Test with particle placed directly on the ball

s_test = [148;706];

% Shown particle on image
subplot(1,2,1);
imshow(reference_frame)
hold on;
plot(s_test(2,:),s_test(1,:),'.');
xlim([0 image_width]);
ylim([0 image_height]);

% Draw rectable around image
xLeft = s_test(2,:) - rect_width/2;
yBottom =s_test(1,:) - rect_height/2;
rectangle('Position',[xLeft,yBottom,rect_width,rect_height],'EdgeColor','b','LineWidth',2);

logical_image=get_rectangle_mask_from_sample(s_test,image_height,image_width,rect_height,rect_width);
maskedRgbImage = bsxfun(@times, reference_frame, cast(logical_image, 'like', reference_frame));

subplot(1,2,2);
imshow(maskedRgbImage);

% Get histogram
histogram = get_histogram(reference_frame,logical_image);

% Get distance
dist_intermediate = sum((target_histogram - histogram).^2,2);
distance = sqrt(sum(dist_intermediate.^2));

fprintf('Distance: %d\n',distance);


%% Calculate particle weights (for single frame)

% Measurement noise
sigma = 0.5;

% distances and particle weights
distances = Inf(1,M);
particle_weights = zeros(1,M);

tic()
% Loop over all histograms
for hist_index = 1:M
    
    % Get image mask
    [logical_image, out_of_image] = get_rectangle_mask_from_sample(S(:,hist_index),image_height,image_width,rect_height,rect_width);
    
    % Check if particle is out of image
    if out_of_image
        distances(hist_index) = distances(hist_index);
        continue;
    end

    % Get distance
    histogram = get_histogram(reference_frame,logical_image);
    dist_intermediate = sum((target_histogram - histogram).^2,2);
    distance = sqrt(sum(dist_intermediate.^2));
    distances(hist_index) = distance;
end

% Normalise distances
% distances = distances/max(distances);

weights = exp(-distances.^2/(2*sigma^2))/(2*pi*sigma);
weights = weights/sum(weights);

time = toc();
sprintf('Time taken to compute weigths %.2f',time)

%% Track static ball
R = diag([10 10]);                                  % process noise 

% Generate initial particles set
M = 200;                % number of particles
S = zeros(3,M);         % set of particles  

% distances and particle weights
weights = zeros(1,M);
distances = zeros(1,M);
particle_weights = zeros(1,M);

S(1,:) = rand(1,M)*(image_height-1)+1;       % y     
S(2,:) = rand(1,M)*(image_width-1)+1;        % x

% Size of rectangles to draw
rect_width = 32;
rect_height = 32; 

for i = 1:100
    
    hold off
    imshow(reference_frame);
    hold on;
    
    % Predict motion of particles
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

    % Update weights
    for hist_index = 1:M
        [logical_image, out_of_image] = get_rectangle_mask_from_sample(S(:,hist_index),image_height,image_width,rect_height,rect_width);
        
        % Check if particle is out of image
        if out_of_image
            distances(hist_index) = distances(hist_index);
            continue;
        end
    
        % Get distance
        histogram = get_histogram(reference_frame,logical_image);
        dist_intermediate = sum((target_histogram - histogram).^2,2);
        distance = sqrt(sum(dist_intermediate.^2));
        distances(hist_index) = distance;
    end
    
    % update weights
    weights = exp(-distances.^2/(2*sigma^2))/(2*pi*sigma);
    S(3,:) = weights/sum(weights);

    % Perform resampling every 5 steps
    S = pf_systematic_resample(S,M);


    pause(1/60);  
    
end

