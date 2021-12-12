% This file is for systematic resampling
%
% "3D Trajectory Prediction of Basketball Shot Using Filtering Techniques
% and Computer Vision" project. The project is a self picked topic for implementation
% in the Appied Estimation course at% the KTH Royal Institute of Technology in 2021.
%
% Authors : 
% Matthew William Lock (mwlock@kth.se)
% Miguel Garcia Naude (magn2@kth.se)

% This function performs systematic re-sampling
% Inputs:   
%           S_bar(t):       3XM         Particles
%           M                           Number of particles
% Outputs:
%           S(t):           3XM
function S = pf_systematic_resample(S_bar,M)
    
    S = zeros(3,M);
    cdf = cumsum(S_bar(3,:));

    r0 = rand()*1/M;

    for m = 1:M 
        i=find(cdf>=(r0 + (m-1)/M), 1 );
        S(:,m)=[S_bar(1:2,i);1/M];
    end

end