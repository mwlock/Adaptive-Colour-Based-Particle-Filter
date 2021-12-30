% Finding ball in video
%
% This file is used for the prediction step of the particle filter, where
% the the motion model using a random walk model. The size of the region is
% adjusted through the same model.
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
%           S(t-1)            7XN           input particle set
%           R                 6X4           process model covariance matrix          
% Outputs:   
%           S_bar(t)          7XN           output particle set
function [S_bar] = pf_predict_motion_and_colour(S,R,M,alpha_min,alpha_max)

    % Noise (4xM) (to generate noise for each state hypothesis/particle)
    noise = [sqrt(R) * randn(6,M); zeros(1,M)];
    
    % Propogate particles though motion model
    S(1,:) = round(S(1,:) + S(3,:));
    S(2,:) = round(S(2,:) + S(4,:));
    S(5,:) = S(5,:) + S(6,:);
    
    % Predicted Particle Set
    S_bar = S + noise;

    % Check for limits
    smaller = S_bar(5,:) < alpha_min;
    bigger = S_bar(5,:) > alpha_max;
    S_bar(5,smaller) = alpha_min;
    S_bar(5,bigger) = alpha_max;

end