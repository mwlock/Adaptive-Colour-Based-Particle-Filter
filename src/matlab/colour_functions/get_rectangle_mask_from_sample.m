% This file is used to generate image masks from samples
% Mask returned is a rectangle of specified domensions, with the center
% determined by the sample position
%
% "3D Trajectory Prediction of Basketball Shot Using Filtering Techniques
% and Computer Vision" project. The project is a self picked topic for implementation
% in the Appied Estimation course at% the KTH Royal Institute of Technology in 2021.
%
% Authors : 
% Matthew William Lock (mwlock@kth.se)
% Miguel Garcia Naude (magn2@kth.se)

function [mask,out_of_image] = get_rectangle_mask_from_sample(sample,image_height,image_width,rect_height,rect_width)
    
    mask = false(image_height,image_width);

    out_of_image = false;   % Boolean that keeps track if particle is inside image
    logical_image_1 = false(image_height,image_width);  
    logical_image_2 = false(image_height,image_width);
    
    % Create regions of interest
    region_x = round(sample(2,:)-rect_width/2)+(0:rect_width-1);
    region_y = round(sample(1,:)-rect_height/2)+(0:rect_height-1);
    
    % Check if regions of interest are out of bounds and correct
    region_x = region_x ((image_width+1 > region_x) & (region_x > 0));
    region_y = region_y ((image_height+1 > region_y) & (region_y > 0));

    % Check if sample is outside of image
    if isempty(region_x) || isempty(region_y)
        out_of_image = true;
        return;
    end
    
    % Create Mask using logical values of region_x & region_y
    logical_image_1(:,region_x) = true;
    logical_image_2(region_y,:) = true;
    mask = and(logical_image_1,logical_image_2);

end