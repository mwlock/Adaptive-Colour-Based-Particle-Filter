% Finding ball in video
%
% This file is used for the prediction step of the particle filter, where
% the the motion model uses constant velocity in the x and y directions and
% adds noise.
%
% "3D Trajectory Prediction of Basketball Shot Using Filtering Techniques
% and Computer Vision" project. The project is a self picked topic for implementation
% in the Appied Estimation course at% the KTH Royal Institute of Technology in 2021.
%
% Authors : 
% Matthew William Lock (mwlock@kth.se)
% Miguel Garcia Naude (magn2@kth.se)

% This function performs the prediction step.
% Inputs:
%           S(t-1)            5XN           input particle set
%           R                 4X4           process model covariance matrix          
% Outputs:   
%           S_bar(t)          5XN           output particle set
function [S_bar] = pf_predict_noise_and_region_size(S,R,M,roi_min,roi_max)

    % Noise (4xM) (to generate noise for each state hypothesis/particle)
    noise = [sqrt(R) * randn(4,M); zeros(1,M)];

    % Predicted Particle Set
    S_bar = S + noise;
    
    % Check for y limits
    smaller_y = S_bar(3,:) < roi_min;
    bigger_y = S_bar(3,:) > roi_max;
    S_bar(3,smaller_y) = roi_min;
    S_bar(3,bigger_y) = roi_max;

    % Check for x limits
    smaller_x = S_bar(4,:) < roi_min;
    bigger_x = S_bar(4,:) > roi_max;
    S_bar(4,smaller_x) = roi_min;
    S_bar(4,bigger_x) = roi_max;
    
end