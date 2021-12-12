% Finding ball in video
%
% This file is used for calculating the weights of particles in a particle
% filter.
%
% "3D Trajectory Prediction of Basketball Shot Using Filtering Techniques
% and Computer Vision" project. The project is a self picked topic for implementation
% in the Appied Estimation course at% the KTH Royal Institute of Technology in 2021.
%
% Authors : 
% Matthew William Lock (mwlock@kth.se)
% Miguel Garcia Naude (magn2@kth.se)

% This function calcultes the weights for each particle based on the
% observation likelihood
%           S_bar(t)            3XM
% Outputs: 
%           S_bar(t)            4XM
function S_bar = pf_weights(S_bar, outlier)

    
    

    % Psi => for each observation, this is the most likely set of particles
    Psi= Psi(1,find(~outlier),:);
    
    product =prod(Psi,2);
    weights=product/sum(product);

    % Note: these are the predicted weights
    S_bar(4,:)=weights;

end