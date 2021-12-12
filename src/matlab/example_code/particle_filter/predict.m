% This function performs the prediction step.
% Inputs:
%           S(t-1)            4XN
%           v                 1X1
%           omega             1X1
% Outputs:   
%           S_bar(t)          4XN
function [S_bar] = predict(S, v, omega, delta_t)

    % Comment out any S_bar(3, :) = mod(S_bar(3,:)+pi,2*pi) - pi before
    % running the test script as it will cause a conflict with the test
    % function. If your function passes, uncomment again for the
    % simulation.

    global R % covariance matrix of motion model | shape 3X3
    global M % number of particles
    
    % YOUR IMPLEMENTATION
    dx = v * delta_t * cos(S(3,:));
    dy = v * delta_t * sin(S(3,:));
    dtheta = ones(1,M) * omega * delta_t;

    dzero = zeros(1,M);

    % U (4x1)
    u = [dx;
         dy;
         dtheta;
         dzero];

    % Noise (4xM) (to generate noise for each state hypothesis/particle)
    noise = [R * randn(3,M);
             zeros(1,M)];
    
    % Predicted Particle Set
    S_bar = S + u + noise;
    
end