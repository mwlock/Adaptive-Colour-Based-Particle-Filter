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

frames = uint8(zeros([dimensions numFrames]));
hsv_frames = zeros([dimensions numFrames]);

blur_std = 0.5;

% Load frames
for i = 1:numFrames
    frames(:,:,:,i) = uint8(imresize(read(v,i),scale));
    % hsv_frames(:,:,:,i) = rgb2hsv(imgaussfilt(frames(:,:,:,i),blur_std));
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

R = diag([15 15 5 5])*scale;                                  % process noise 

M = 200;                % number of particles
S = zeros(5,M);         % set of particles  

% distances and particle weights
weights = zeros(1,M);
distances = ones(1,M);

% Size of rectangles to draw
rect_width = 67*scale;
rect_height = 42*scale; 

% Set region of interest (ROI) limits
roi_min = 10;
roi_max = image_width;

% Init particles
S(1,:) = rand(1,M)*(image_height-1)+1;      % y     
S(2,:) = rand(1,M)*(image_width-1)+1;       % x
S(3,:) = ones(1,M)*rect_height;             % roc height
S(4,:) = ones(1,M)*rect_width;              % roc width     
S(5,:) = ones(1,M)/M; 

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
mean_state_observation_prob_thresh = 0.5;
alpha = 0.05;

% Retain original target distribution
target_histogram = original_target_histogram;

% Decide mode for mean hist calculation
% mean_hist_calc_mode = 0;        % mean histogram
mean_hist_calc_mode = 1;        % hist_at_mean

% Store histograms
histograms = zeros(M,size(target_histogram,2));

for i = 40:numFrames
    
    % Get image
    image = frames(:,:,:,i);
    hsv_image = hsv_frames(:,:,:,i);
    
    % Plot image
    % subplot(1,2,1);
    hold off;
    imshow(image);
    hold on;
    
    % Predict motion of particles
    S = pf_predict_noise_and_region_size(S,R,M,roi_min,roi_max);
    
    % Plot particles
    plot(S(2,:),S(1,:),'.');
    xlim([0 image_width]);
    ylim([0 image_height]);

    % Plot particle regions
    %     for m = 1:M
    %         % Get region size
    %         rect_width = S(4,m);
    %         rect_height = S(3,m);
    %         xLeft = S(2,m) - rect_width/2;
    %         yBottom = S(1,m) - rect_height/2;
    %         rectangle('Position',[xLeft,yBottom,rect_width,rect_height],'EdgeColor','b','LineWidth',1);
    %     end   

    % Update weights
    for hist_index = 1:M

        % Get region size
        rect_width = S(4,hist_index);
        rect_height = S(3,hist_index);

        % Region coordinates
        coordinates = [S(1,hist_index),S(2,hist_index)];

        % Get mask
        [logical_image, out_of_image] = get_rectangle_mask_from_sample(S(:,hist_index),image_height,image_width,rect_height,rect_width);
        
        % Plot particle regions
        %         xLeft = round(S(2,hist_index) - rect_width/2);
        %         yBottom = round(S(1,hist_index) - rect_height/2);
        %         rectangle('Position',[xLeft,yBottom,rect_width,rect_height],'EdgeColor','b','LineWidth',1);
        
        % Check if particle is out of image
        if out_of_image
            distances(hist_index) = distances(hist_index);
            continue;
        end
    
        % Get histogram
        histogram = get_weighted_histogram(hsv_image,logical_image,coordinates,5);
        histogram = reshape(histogram',1,[]);
        histograms(hist_index,:) = histogram;

        % Get distance   
        distance = bhattacharyya_distance(target_histogram,histogram);
        distances(hist_index) = distance;
    end

    fprintf('Min distance %.2f \n\n', min(distances));
    
    % update weights
    weights = 1/(sqrt(2*pi)*sigma)*exp(-distances.^2/(2*sigma^2));
    weights = weights/sum(weights);
    S(5,:) = weights;

    % Get mean position 
    mean_x = round(sum(S(2,:).*S(5,:)));
    mean_y = round(sum(S(1,:).*S(5,:)));
    rect_height = round(sum(S(3,:).*S(5,:)));
    rect_width = round(sum(S(4,:).*S(5,:)));

    % Plot mean positons
    xLeft = round(mean_x - rect_width/2);
    yBottom = round(mean_y - rect_height/2);
    rectangle('Position',[xLeft,yBottom,rect_width,rect_height],'EdgeColor','g','LineWidth',1);

    % Calculate mean state histogram
    if mean_hist_calc_mode == 0
        mean_state_histogram = sum(bsxfun(@times, histograms, S(5,:)'),1);
    elseif mean_hist_calc_mode == 1
        [logical_image, out_of_image] = get_rectangle_mask_from_sample([mean_y;mean_x],image_height,image_width,rect_height,rect_width);
        % imshow(bsxfun(@times, image, cast(logical_image, 'like', image)));    

        % Check if particle is out of image
        if out_of_image
            continue;
        end
    
        % Get histogram
        mean_state_histogram = get_histogram(hsv_image,logical_image,false);
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
        S(5,:) = ones(1,M)*1/M;
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
