%% Track ball in video
% DO NOT CHANGE PARAMTERS

% Clean environment and load video
clc; clear;
video = "horse_v1.mp4";

% Get number of frames
v = VideoReader(video); numFrames = 0;
while hasFrame(v)
    readFrame(v);
    numFrames = numFrames + 1;
    fprintf('Number of frames : %d \n',numFrames);
end

% Set downsampling size
scale = 0.3;

% Get deminsions of frame
ref_frame = read(v,1);
ref_frame = imresize(ref_frame,scale);
dimensions = size(ref_frame);

% Preallocate matrices
frames = uint8(zeros([dimensions numFrames]));
hsv_frames = zeros([dimensions numFrames]);

% Load frames
for i = 1:numFrames
    frames(:,:,:,i) = uint8(imresize(read(v,i),scale));
    hsv_frames(:,:,:,i) = rgb2hsv(frames(:,:,:,i));
    fprintf('Progress \t%d/%d \t(%.2f%%) \n',i,numFrames,i/numFrames*100);
end

%% Get colour params
reference_frame= frames(:,:,:,15);
reference_frame_hsv = hsv_frames(:,:,:,15);

% Get image parameters
image_height = size(reference_frame_hsv,1);
image_width = size(reference_frame_hsv,2);

% Load binary image
file = matfile('bee_in_horsev1_binary_frame_15.mat');
binaryImage = imresize(file.binaryImage,scale);

% Get reference histogram
target_histogram = get_histogram(reference_frame_hsv,binaryImage,false);
target_histogram = reshape(target_histogram',1,[]);
original_target_histogram = target_histogram;

% Show pixels of interest
imshow(bsxfun(@times, reference_frame, cast(binaryImage, 'like', reference_frame)));

%% Track with dynamic target distribution

R = diag([100 100])*scale;                                  % process noise 

M = 200;                % number of particles
S = zeros(3,M);         % set of particles  

% distances and particle weights
weights = zeros(1,M);
distances = ones(1,M);
 
S(1,:) = rand(1,M)*(image_height-1)+1;       % y     
S(2,:) = rand(1,M)*(image_width-1)+1;        % x

% Size of rectangles to draw
rect_width = 67*scale;
rect_height = 42*scale; 

% Measurement noise
sigma = 0.2;

% Resampling threshold
resampling_thresh = 0.2;

% Keep track of whether the target has converged
tracking = false;

% Enable or disable initialisaiton
reinit_particles = false;

% Track mean state observation probability
% mean_state_observation_probabilities = zeros(1,numFrames);
clear mean_state_observation_probabilities;

% Specify contribution of mean state distribution and probability threshold
mean_state_observation_prob_max=1;          % used for graphing
mean_state_observation_prob_thresh = 0.7;
alpha = 0;

target_histogram = original_target_histogram;

% Decide mode for mean hist calculation
% mean_hist_calc_mode = 0;        % mean histogram
mean_hist_calc_mode = 1;        % hist_at_mean

for i = 1:numFrames
    
    % Get image
    image = frames(:,:,:,i);
    hsv_image = hsv_frames(:,:,:,i);
    
    % Plot image
    % subplot(1,2,1);
    hold off
    imshow(image);
    hold on;
    
    % Predict motion of particles
    S = predict_noise(S,R,M);
    
    plot(S(2,:),S(1,:),'.');
    xlim([0 image_width]);
    ylim([0 image_height]);

    % Store histograms
    histograms = zeros(M,size(target_histogram,2));

    % Update weights
    for hist_index = 1:M
        [logical_image, out_of_image] = get_rectangle_mask_from_sample(S(:,hist_index),image_height,image_width,rect_height,rect_width);
        
        % Check if particle is out of image
        if out_of_image
            distances(hist_index) = distances(hist_index);
            continue;
        end
    
        % Get histogram
        histogram = get_histogram(hsv_image,logical_image,false);
        histogram = reshape(histogram',1,[]);
        histograms(hist_index,:) = histogram;

        % Get distance   
        distance = bhattacharyya_distance(target_histogram,histogram);
        distances(hist_index) = distance;
    end

    fprintf('Min distance %.2f \n\n', min(distances));
    
    % update weights
    weights = 1/(sqrt(2*pi)*sigma)*exp(-distances.^2/(2*sigma^2));
    S(3,:) = weights/sum(weights);

    % Get mean position 
    mean_x = sum(S(2,:).*S(3,:));
    mean_y = sum(S(1,:).*S(3,:));

    % Plot mean positons
    xLeft = mean_x - rect_width/2;
    yBottom = mean_y - rect_height/2;
    rectangle('Position',[xLeft,yBottom,rect_width,rect_height],'EdgeColor','g','LineWidth',1);

    % Calculate mean state histogram
    if mean_hist_calc_mode == 0
        mean_state_histogram = sum(bsxfun(@times, histograms, S(3,:)'),1);
    elseif mean_hist_calc_mode == 1
        [logical_image_mean_hist, out_of_image] = get_rectangle_mask_from_sample([mean_y;mean_x],image_height,image_width,rect_height,rect_width);
        % imshow(bsxfun(@times, image, cast(logical_image_mean_hist, 'like', image)));    

        % Check if particle is out of image
        if out_of_image
            continue;
        end
    
        % Get histogram
        mean_state_histogram = get_histogram(hsv_image,logical_image_mean_hist,false);
        mean_state_histogram = reshape(mean_state_histogram',1,[]);
    end

    % Calculate distance to mean distribution
    mean_state_hist_dist = bhattacharyya_distance(target_histogram,mean_state_histogram);

    % Calculate mean state observation probability
    mean_state_observation_prob = exp(-mean_state_hist_dist.^2/(2*sigma^2));
    mean_state_observation_probabilities(i)=mean_state_observation_prob; %#ok<SAGROW> 

    % Apply histogram update
    if mean_state_observation_prob > mean_state_observation_prob_thresh
        rectangle('Position',[xLeft,yBottom,rect_width,rect_height],'EdgeColor','r','LineWidth',1);
        target_histogram = (1-alpha) * target_histogram + alpha * mean_state_histogram;
    end

    % Perform resampling 
    if min(distances)< resampling_thresh || true
        S = pf_systematic_resample(S,M);
    end

    std_x = std(S(2,:));
    std_y = std(S(1,:));

    % Check if tracking
    if std_x < 10 &&  std_y < 10 && ~tracking
        tracking = true;
    end 

    if std_x >= 10 &&  std_y >= 10 && tracking
        tracking = false;
    end

    % Check for standard variation of particles
    if std_x > 15 &&  std_y > 15 && tracking && reinit_particles
        S(1,:) = rand(1,M)*(image_height-1)+1;       % y     
        S(2,:) = rand(1,M)*(image_width-1)+1;        % x
        S(3,:) = ones(1,M)*1/M;
        tracking = false;
    end

    fprintf('Frame %d\nTracking %d\nstd x : %0.3f\nstd y: %0.3f\n',i,tracking,std_x,std_y);

    % plot observation prob    
    %     subplot(1,2,2);
    %     hold off;
    %     plot(mean_state_observation_probabilities);
    %     mean_state_observation_prob_max = max(mean_state_observation_prob_max,mean_state_observation_prob);
    %     axis([0 numFrames 0 mean_state_observation_prob_max]);
    %     drawnow


    pause(1/160);  
    
end

figure;
plot(mean_state_observation_probabilities);
