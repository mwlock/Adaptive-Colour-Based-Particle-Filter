% This file is to return colour histogram of region in an image
%
% "3D Trajectory Prediction of Basketball Shot Using Filtering Techniques
% and Computer Vision" project. The project is a self picked topic for implementation
% in the Appied Estimation course at% the KTH Royal Institute of Technology in 2021.
%
% Authors : 
% Matthew William Lock (mwlock@kth.se)
% Miguel Garcia Naude (magn2@kth.se)

% observation likelihood
%           image                
%           binaryImage
% Outputs: 
%           histogram            3X8    

function histogram = get_histogram(image, binaryImage,convert_to_hsv)

    % Check i fconversion to HSV is needed
    if nargin == 2
        convert_to_hsv = true;
    end
    
    % Get masked region
    reference_frame = reshape(image, [], 3);
    mask_reshaped = reshape(binaryImage, [], 1);
    masked_pixels = reference_frame(mask_reshaped,:);
    masked_pixels = reshape(masked_pixels, [], 1, 3);
    
    % convert RGB to HSV
    if convert_to_hsv
        masked_pixels = rgb2hsv(masked_pixels);
    end
    
    % Get Colours histogram

    [h_counts,~] = imhist(masked_pixels(:,:,1),8);
    h_counts = h_counts';
    h_counts(8)=0;

    [s_counts,~] = imhist(masked_pixels(:,:,2),8);
    s_counts = s_counts';
    s_counts(8)=0;

    [v_counts,~] = imhist(masked_pixels(:,:,3),1);
    v_counts = v_counts';
    v_counts(8)=0;
    % v_counts = zeros(1,8);
    
    % Generate target histogram
    histogram = [h_counts;s_counts;v_counts];
    histogram = histogram / sum(histogram(1,:),2);  % divide by total number of pixels