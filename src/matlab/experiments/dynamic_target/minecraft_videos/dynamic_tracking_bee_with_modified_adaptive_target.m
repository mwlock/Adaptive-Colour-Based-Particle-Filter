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

% Load frames
for i = 1:numFrames
    frames(:,:,:,i) = uint8(imresize(read(v,i),scale));
    hsv_frames(:,:,:,i) = rgb2hsv(frames(:,:,:,i));
    fprintf('Progress \t%d/%d \t(%.2f%%) \n',i,numFrames,i/numFrames*100);
end

%% Get colour params
reference_frame = frames(:,:,:,15);

% Get image parameters
image_height = size(reference_frame,1);
image_width = size(reference_frame,2);

% Load binary image
file = matfile('bee_in_horsev1_binary_frame_15.mat');
binaryImage = imresize(file.binaryImage,scale);

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
[v_counts,binLocations_v] = imhist(masked_pixels_HSV(:,:,3),1);
v_counts = v_counts';
v_counts(8)=0;

target_histogram = [h_counts;s_counts;v_counts];
target_histogram = target_histogram / sum(target_histogram(1,:),2);
target_histogram = reshape(target_histogram',1,[]);
original_target_histogram = target_histogram;

% Reshape reference back to original
reference_frame = reshape(reference_frame,image_height,image_width,3);
imshow(bsxfun(@times, reference_frame, cast(binaryImage, 'like', reference_frame)));

%% Track with dynamic target distribution

R = diag([50 50 0 0 0 0])*scale;                                  % process noise 

M = 200;                % number of particles
S = zeros(7,M);         % set of particles  

% distances and particle weights
weights = ones(1,M)*1/M;
distances = ones(1,M);

% Specify contribution of mean state distribution and probability threshold
mean_state_observation_prob_max=1;          % used for graphing
mean_state_observation_prob_thresh = 0.7;

% Alpha limits
alpha_min = 0.05;
alpha_max = 0.2;
beta = (alpha_max - alpha_min);

% Init random particles
S(1,:) = rand(1,M)*(image_height-1)+1;                          % y     
S(2,:) = rand(1,M)*(image_width-1)+1;                           % x
S(3,:) = (rand(1,M)*image_height-image_height/2)/20*0;          % y'
S(4,:) = (rand(1,M)*image_width-image_width/2)/20*0;            % x'
S(5,:) = alpha_min + (alpha_min_max-alpha_min) .* rand(1,M);    % alpha
S(6,:) = (rand(1,M)*beta-beta/2)/10;                            % alpha'
S(7,:) = weights;

% Size of rectangles to draw
rect_width = 67*scale;
rect_height = 42*scale; 

% Measurement noise
sigma = 0.1;

% Resampling threshold
resampling_thresh = 0.2;

% Keep track of whether the target has converged
tracking = false;

% Enable or disable initialisaiton
reinit_particles = false;

% Track mean state observation probability
% mean_state_observation_probabilities = zeros(1,numFrames);
clear mean_state_observation_probabilities;

target_histogram = original_target_histogram;

% Decide mode for mean hist calculation
% mean_hist_calc_mode = 0;        % mean histogram
mean_hist_calc_mode = 1;        % hist_at_mean

% Store histograms
histograms = zeros(M,8*3);
hypothesis_histograms = zeros(M,8*3);

for i = 1:numFrames
    
    % Get image
    image = frames(:,:,:,i);
    hsv_image = hsv_frames(:,:,:,i);
    
    % Plot image
    %     subplot(1,2,1);
    hold off
    imshow(image);
    hold on;
    
    % Predict motion of particles
    S = pf_predict_motion_and_colour(S,R,M,alpha_min,alpha_max);
    
    % Plot particles
    plot(S(2,:),S(1,:),'.');
    xlim([0 image_width]);
    ylim([0 image_height]);

    % Get mean position 
    if i == 1
        mean_x = sum(S(2,:).*S(7,:));
        mean_y = sum(S(1,:).*S(7,:));

        % Calculate mean state histogram
        if mean_hist_calc_mode == 0
            mean_state_histogram = sum(bsxfun(@times, histograms, S(7,:)'),1);
        elseif mean_hist_calc_mode == 1
            [logical_image, out_of_image] = get_rectangle_mask_from_sample([mean_y;mean_x],image_height,image_width,rect_height,rect_width);
            % imshow(bsxfun(@times, reference_frame, cast(logical_image, 'like', reference_frame)));    
    
            % Check if particle is out of image
            if out_of_image
                continue;
            end
        
            % Get histogram
            mean_state_histogram = get_histogram(hsv_image,logical_image,false);
            mean_state_histogram = reshape(histogram',1,[]);
        end

    end

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

        % Calculate hypothesis histogram
        hypothesis_target = (1-S(5,hist_index))*target_histogram + S(5,hist_index)*mean_state_histogram;
        hypothesis_histograms(hist_index,:) = hypothesis_target;

        % Get distance   
        distance = bhattacharyya_distance(hypothesis_target,histogram);
        distances(hist_index) = distance;
        
        if ~isreal(distance)
            disp('kill me');
        end

    end    

    fprintf('Min distance %.2f \n\n', min(distances));
    
    % update weights
    %     weights = 1/(sqrt(2*pi)*sigma)*exp(-distances.^2/(2*sigma^2));
    weights = exp(-distances.^2/(2*sigma^2));
    S(7,:) = weights/sum(weights);

    % Get mean position 
    mean_x = sum(S(2,:).*S(7,:));
    mean_y = sum(S(1,:).*S(7,:));

    % Plot mean positons
    xLeft = mean_x - rect_width/2;
    yBottom = mean_y - rect_height/2;
    rectangle('Position',[xLeft,yBottom,rect_width,rect_height],'EdgeColor','g','LineWidth',1);
    
    % Update target histogram?
    target_histogram = (1 - sum(S(5,:).*S(7,:)))*target_histogram + sum(S(5,:).*S(7,:))* mean_state_histogram ;

    % Calculate mean state histogram
    if mean_hist_calc_mode == 0
        mean_state_histogram = sum(bsxfun(@times, histograms, S(7,:)'),1);
    elseif mean_hist_calc_mode == 1
        [logical_image, out_of_image] = get_rectangle_mask_from_sample([mean_y;mean_x],image_height,image_width,rect_height,rect_width);
        % imshow(bsxfun(@times, reference_frame, cast(logical_image, 'like', reference_frame)));    

        % Check if particle is out of image
        if out_of_image
            continue;
        end
    
        % Get histogram
        mean_state_histogram = get_histogram(hsv_image,logical_image,false);
        mean_state_histogram = reshape(histogram',1,[]);
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
        S(7,:) = ones(1,M)*1/M;
        tracking = false;
    end

    fprintf('Tracking %d\nstd x : %0.3f\nstd y: %0.3f\n',tracking,std_x,std_y);

    % plot observation prob    
    %     subplot(1,2,2);
    %     hold off;
    %     plot(mean_state_observation_probabilities);
    %     mean_state_observation_prob_max = max(mean_state_observation_prob_max,mean_state_observation_prob);
    %     axis([0 numFrames 0 mean_state_observation_prob_max]);
    %     drawnow


    pause(1/160);  
    
end