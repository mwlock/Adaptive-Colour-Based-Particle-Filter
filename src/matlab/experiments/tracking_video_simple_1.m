%% Track static ball

% Clean environment and load video
clc; clear;
videos = ["shot_1_vid_low_res.mp4"];
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

reference_frame = read(v,55);

% Get image parameters
image_height = size(reference_frame,1);
image_width = size(reference_frame,2);

% Get target distribution
% Get Colours histogram

% Load binary image
file = matfile('shot_1_vid_high_res_binary_frame_55.mat');
binaryImage = file.binaryImage;

% Get masked region
reference_frame = reshape(reference_frame, [], 3);
mask_reshaped = reshape(binaryImage, [], 1);
masked_pixels = reference_frame(mask_reshaped,:);
masked_pixels = reshape(masked_pixels, [], 1, 3);

% convert RGB to HSV
masked_pixels_HSV = rgb2hsv(masked_pixels);

% Generate target histogram
[h_counts,binLocations_h] = imhist(masked_pixels_HSV(:,:,1),8);
h_counts = h_counts';
imhist(masked_pixels_HSV(:,:,1),8);
[s_counts,binLocations_s] = imhist(masked_pixels_HSV(:,:,2),8);
s_counts = s_counts';
[v_counts,binLocations_v] = imhist(masked_pixels_HSV(:,:,3),4);
v_counts = v_counts';
v_counts(8)=0;
v_counts=zeros(1,8);
target_histogram = [h_counts;s_counts;v_counts];
target_histogram = target_histogram / sum(target_histogram(1,:),2);

% Reshape reference back to original
reference_frame = reshape(reference_frame,image_height,image_width,3);

R = diag([100 100]);                                  % process noise 

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
rect_width = 30;
rect_height = 30; 

% Measurement noise
sigma = 0.05;

% Keep track of whether the target has converged
tracking = false;

% Enable or disable initialisaiton
reinit_particles = false;

for i = 1:numFrames

    image = frames(:,:,:,i);
    
    hold off
    imshow(image);
    hold on;
    
    % Predict motion of particles
    S = predict_noise(S,R,M);
    
    plot(S(2,:),S(1,:),'.');
    xlim([0 image_width]);
    ylim([0 image_height]);

    % Update weights
    for hist_index = 1:M
        [logical_image, out_of_image] = get_rectangle_mask_from_sample(S(:,hist_index),image_height,image_width,rect_height,rect_width);
        
        % Check if particle is out of image
        if out_of_image
            distances(hist_index) = distances(hist_index);
            continue;
        end
    
        % Get distance
        histogram = get_histogram(image,logical_image);
        dist_intermediate = sum((target_histogram - histogram).^2,2);
        distance = sqrt(sum(dist_intermediate.^2));
        distances(hist_index) = distance;
    end
    
    % update weights
    weights = exp(-distances.^2/(2*sigma^2))/(2*pi*sigma);
    S(3,:) = weights/sum(weights);

    % Perform resampling every 5 steps
    if min(distances)<0.15
        S = pf_systematic_resample(S,M);
    end

    % plot estimated position 
    mean_x = sum(S(2,:).*S(3,:));
    mean_y = sum(S(1,:).*S(3,:));

    xLeft = mean_x - rect_width/2;
    yBottom = mean_y - rect_height/2;
    rectangle('Position',[xLeft,yBottom,rect_width,rect_height],'EdgeColor','g','LineWidth',3);

    std_x = std(S(2,:));
    std_y = std(S(1,:));

    % Check if tracking
    if std_x < 10 &&  std_y < 10 && ~tracking
        tracking = true;
    else
        tracking = false;
    end

    % Check for standard variation of particles
    if std_x > 15 &&  std_y > 15 && tracking && reinit_particles
        S(1,:) = rand(1,M)*(image_height-1)+1;       % y     
        S(2,:) = rand(1,M)*(image_width-1)+1;        % x
        S(3,:) = ones(1,M)*1/M;
        tracking = false;
    end

    fprintf('Tracking %d\nstd x : %0.3f\nstd y: %0.3f\n',tracking,std_x,std_y);

    pause(1/60);  
    
end
