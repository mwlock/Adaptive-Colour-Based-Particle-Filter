% Finding ball in video
%
% This file is used for the prediction step of the particle filter, where
% the the motion model used the identity matrix for noise and simply adds
% noise to the previous state.
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
%           S(t-1)            3XN           input particle set
%           R                 2X2           process model covariance matrix          
% Outputs:   
%           S_bar(t)          3XN           output particle set
function [S_bar] = predict_noise(S,R,M)

    dzero = zeros(1,M);

    % Noise (4xM) (to generate noise for each state hypothesis/particle)
    noise = [sqrt(R) * randn(2,M); zeros(1,M)];
    
    % Predicted Particle Set
    S_bar = round(S + noise);
    
end