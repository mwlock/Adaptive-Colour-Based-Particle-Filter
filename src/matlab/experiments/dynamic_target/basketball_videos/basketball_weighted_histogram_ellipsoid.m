%% Track ball in video
% DO NOT CHANGE PARAMTERS

% Clean environment and load video
clc; clear;
video = "shot_1_vid_high_res.mp4";

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

% Get image parameters
image_height = dimensions(1);
image_width = dimensions(2);

frames = uint8(zeros([dimensions numFrames]));
hsv_frames = zeros([dimensions numFrames]);

% Apply blur to image
blur_std = 0.5;

% Load frames
for i = 1:numFrames
    frames(:,:,:,i) = uint8(imresize(read(v,i),scale));
    % hsv_frames(:,:,:,i) = rgb2hsv(imgaussfilt(frames(:,:,:,i),blur_std));
    hsv_frames(:,:,:,i) = rgb2hsv(frames(:,:,:,i));
    fprintf('Progress \t%d/%d \t(%.2f%%) \n',i,numFrames,i/numFrames*100);
end

% Create meshgrid
x_max = size(ref_frame,2);
y_max = size(ref_frame,1);
x = 1:x_max;
y = 1:y_max;
[X,Y] = meshgrid(x,y);

%% Tune initial histogram
reference_frame= frames(:,:,:,55);
reference_frame_hsv = hsv_frames(:,:,:,55);

% Parameters
center_x_coordinate = 212;
center_y_coordinate = 45;
Hx_initial = 5;
Hy_initial = 5; 

% Show tuning
imshow(reference_frame);
hold on;
ellipse(Hx_initial,Hy_initial,0,center_x_coordinate,center_y_coordinate,'g');

%% Show mask
logical_image = get_ellipse_mask(center_x_coordinate,center_y_coordinate,Hx_initial,Hy_initial,X,Y);
imshow(bsxfun(@times, reference_frame, cast(logical_image, 'like', reference_frame)));

%% Get reference histogram
a = sqrt((Hx_initial/2)^2 + (Hy_initial/2)^2);
tic;
target_histogram = get_weighted_histogram(reference_frame_hsv,logical_image,[center_y_coordinate,center_x_coordinate],a,X,Y);
toc;
target_histogram = reshape(target_histogram',1,[]);
original_target_histogram = target_histogram;
%% Track with dynamic target distribution

R = diag([20 20 5 5])*scale; % process noise 

M = 400;                % number of particles
S = zeros(5,M);         % set of particles  

% distances and particle weights
weights = zeros(1,M);
distances = ones(1,M);

% Size of rectangles to draw
% rect_width = 67*scale;
% rect_height = 42*scale; 
Hx = Hx_initial;
Hy = Hy_initial; 

% Set region of interest (ROI) limits
roi_min = 0;
roi_max = image_width;

% Init particles
S(1,:) = rand(1,M)*(image_height-1)+1;      % y     
S(2,:) = rand(1,M)*(image_width-1)+1;       % x
S(3,:) = ones(1,M)*Hy;                      % roi height
S(4,:) = ones(1,M)*Hx;                      % roi width     
S(5,:) = ones(1,M)/M;                       % weights

% Measurement noise
sigma = 0.2;

% Resampling threshold
resampling_thresh = 0.4;

% Keep track of whether the target has converged
tracking = false;

% Enable or disable initialisaiton
reinit_particles = false;

% Track mean state observation probability
% mean_state_observation_probabilities = zeros(1,numFrames);
clear mean_state_observation_probabilities;

% Specify contribution of mean state distribution and probability threshold
mean_state_observation_prob_max=1;          % used for graphing
mean_state_observation_prob_thresh = 0.45;
alpha = 0.05;

% Retain original target distribution
target_histogram = original_target_histogram;

% Decide mode for mean hist calculation
% mean_hist_calc_mode = 0;        % mean histogram
mean_hist_calc_mode = 1;        % hist_at_mean

% Store histograms
histograms = zeros(M,size(target_histogram,2));

% Save the mean state
mean_state = zeros(5,M);

% save all samples
all_samples = zeros(5,M,numFrames);

for i = 20:numFrames

    % Get image
    image = frames(:,:,:,i);
    hsv_image = hsv_frames(:,:,:,i);
    
    % Plot image
    subplot(1,2,1);
    hold off;
    imshow(image);
    hold on;
    
    % Predict motion of particles
    S = pf_predict_noise_and_region_size(S,R,M,roi_min,roi_max);
    
    % Plot particles
    % plot(S(2,:),S(1,:),'.');    
    xlim([0 image_width]);
    ylim([0 image_height]);

    % Plot particle regions
%     for m = 1:M
%         % Get region size
%         Hx = S(4,m);
%         Hy = S(3,m);
%         x = S(2,m);
%         y = S(1,m);
%         ellipse(Hx,Hy,0,x,y,'g');
%     end   

    % Update weights
    
    tic;
    for hist_index = 1:M

        % Get region size
        Hx = S(4,hist_index);
        Hy = S(3,hist_index);
        a = sqrt((Hx/2)^2 + (Hy/2)^2);        

        % Region coordinates
        x = S(2,hist_index);
        y = S(1,hist_index);
        coordinates = [y,x];

        % Get mask      
        logical_image = get_ellipse_mask(coordinates(2),coordinates(1),Hx,Hy,X,Y);
        % imshow(bsxfun(@times, image, cast(logical_image, 'like', image)))

        % Plot particle regions
        % ellipse(Hx,Hy,0,x,y,'g');
    
        % Get histogram          
        histogram = get_weighted_histogram(hsv_image,logical_image,coordinates,a,X,Y);

        if sum(histogram ==1)
            distances(hist_index) = 1;
        else
            % Reshape histograme
            histogram = reshape(histogram',1,[]);
            histograms(hist_index,:) = histogram;
    
            % Get distance   
            distance = bhattacharyya_distance(target_histogram,histogram);
            distances(hist_index) = distance;
        end

    end
    toc;

    fprintf('Min distance %.2f \n\n', min(distances));
    
    % update weights
    weights = 1/(sqrt(2*pi)*sigma)*exp(-distances.^2/(2*sigma^2));
    weights = weights/sum(weights);
    S(5,:) = weights;

    % Get mean position 
    mean_x = round(sum(S(2,:).*S(5,:)));
    mean_y = round(sum(S(1,:).*S(5,:)));
    Hy_mean = round(sum(S(3,:).*S(5,:)));
    Hx_mean = round(sum(S(4,:).*S(5,:)));
    a = sqrt((Hy_mean/2)^2 + (Hx_mean/2)^2); 
    coordinates = [mean_y,mean_x];

    % Plot mean positons
    ellipse(Hx_mean,Hy_mean,0,mean_x,mean_y,'g');

    % Calculate mean state histogram
    if mean_hist_calc_mode == 0
        mean_state_histogram = sum(bsxfun(@times, histograms, S(5,:)'),1);
    elseif mean_hist_calc_mode == 1
        logical_image = get_ellipse_mask(coordinates(2),coordinates(1),Hx_mean,Hy_mean,X,Y);
        % imshow(bsxfun(@times, image, cast(logical_image, 'like', image)));   

        % Get histogram
        mean_state_histogram = get_weighted_histogram(hsv_image,logical_image,coordinates,a,X,Y);
        mean_state_histogram = reshape(mean_state_histogram',1,[]);
    end

    % Calculate distance to mean distribution
    mean_state_hist_dist = bhattacharyya_distance(target_histogram,mean_state_histogram);

    % Calculate mean state observation probability
    mean_state_observation_prob = exp(-mean_state_hist_dist.^2/(2*sigma^2));
    mean_state_observation_probabilities(i)=mean_state_observation_prob; %#ok<SAGROW> 

    % Apply histogram update
    if mean_state_observation_prob > mean_state_observation_prob_thresh
        ellipse(Hx_mean,Hy_mean,0,mean_x,mean_y,'r');
        target_histogram = (1-alpha) * target_histogram + alpha * mean_state_histogram;
    end

    % Save the mean state
    mean_state(:,i) = [mean_y;mean_x;Hy_mean;Hx_mean;mean_state_observation_prob];

    % Perform resampling 
    if min(distances)< resampling_thresh || true
        S = pf_systematic_resample(S,M);
    end

    % Save partilces
    all_samples(:,:,i)=S;

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
    subplot(1,2,2);
    hold off;
    plot(mean_state_observation_probabilities);
    mean_state_observation_prob_max = max(mean_state_observation_prob_max,mean_state_observation_prob);
    axis([0 numFrames 0 mean_state_observation_prob_max]);
    drawnow


%     pause(1/160);  
    
end

%% Result plotter

for i = 20:numFrames

    image = frames(:,:,:,i);    
    hold off;
    imshow(image);
    hold on;
    
    % plot samples
    plot(all_samples(2,:,i),all_samples(1,:,i),'.');

    % Draw ellipses around partilces
    %     for m = 1:M
    %         % Get region size
    %         Hx = all_samples(4,m,i);
    %         Hy = all_samples(3,m,i);
    %         x = all_samples(2,m,i);
    %         y = all_samples(1,m,i);
    %         ellipse(Hx,Hy,0,x,y,'g');
    %     end   

    y = mean_state(1,i);
    x = mean_state(2,i);
    Hy = mean_state(3,i);
    Hx = mean_state(4,i);
    observation_prob = mean_state(5,i);
    c = 'g';
    if observation_prob > mean_state_observation_prob_thresh
        c = 'r';
    end

    ellipse(Hx,Hy,0,x,y,c);

    if i == 180
        disp(i);
    end


    pause(1/60);


end

%% plot observation graph

plot(mean_state(5,20:end))
xlim([0 size(mean_state(5,20:end),2)]);
% ylim([0 image_height]);

%% other stuff
xlim([230 576])
ylim([50 300])

%% more stuff
xlim([0 213])
t = xlabel('Frame');
t = ylabel('Mean state bservation probability \pi_{E(S)}');
t = yline(0.45,'r','Observation probability threshold \pi_T');
set(gca,'FontSize',32);
t.FontSize = 32;
