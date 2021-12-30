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
    
    % Parameters
    h_bins = 8;
    s_bins = 8;
    v_bins = 1;
    
    % Matrix sizing
    max_bins = max([h_bins,s_bins,v_bins]);
    bins = zeros(1,max_bins);
    
    % Get Colours histogram
    [h_counts,~] = imhist(masked_pixels(:,:,1),h_bins);
    h_counts = h_counts';
    h = bins;
    h(1:size(h_counts,2)) = h_counts;

    [s_counts,~] = imhist(masked_pixels(:,:,2),s_bins);
    s_counts = s_counts';
    s = bins;
    s(1:size(h_counts,2)) = s_counts;

    [v_counts,~] = imhist(masked_pixels(:,:,3),v_bins);
    v_counts = v_counts';
    v = bins;
    v(1:size(v_counts,2)) = v_counts;
    
    % Generate target histogram
    histogram = [h;s;v];
    histogram = histogram / sum(histogram(1,:),2);  % divide by total number of pixels