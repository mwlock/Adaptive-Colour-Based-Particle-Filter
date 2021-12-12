% This function performs the ML data association
%           S_bar(t)                 4XM
%           z(t)                     2Xn
%           association_ground_truth 1Xn | ground truth landmark ID for
%           every measurement  1 (actually not used)
% Outputs: 
%           outlier                  1Xn
%           Psi(t)                   1XnXM
%           c                        1xnxM (actually not ever used)
function [outlier, Psi, c] = associate(S_bar, z, association_ground_truth)
    if nargin < 3
        association_ground_truth = [];
    end

    global lambda_psi % threshold on average likelihood for outlier detection
    global Q % covariance matrix of the measurement model
    global M % number of particles
    global N % number of landmarks
    global map;
    
    % YOUR IMPLEMENTATION
    n_landmarks = N;
    n_observations = size(z,2);
    n_particles = M;

    % Initiate Storage Variables
    z_hat = zeros(2, n_particles, n_landmarks);
    nu = zeros(size(z_hat)); 
    psi = zeros(n_particles, n_landmarks);
    Psi = zeros(1,n_observations,n_particles);
    outlier = zeros(1,n_observations);
    c = zeros(1,n_observations, n_particles);

    % Loop through all landmarks and use observation model to get predicted
    % measurement for each particle to specific landmark
    % This is was taken out of original loop to increase computation
    % efficiency
    for j = 1:n_landmarks
        z_hat(:,:,j) = observation_model(S_bar, j); % (2 x M x n_landmarks matrix)
    end
    
    % Manipulate Q for operations
    Q_repeated = repmat(diag(Q),[1,n_particles,n_landmarks]);

    % Loop through each observation
    for i = 1:n_observations
        nu(:,:,:) = z(:,i) - z_hat;
        nu(2,:,:) = mod(nu(2,:,:) + pi, 2*pi) - pi;
        
        psi(:,:) = 1/(2*pi*det(Q)^(1/2)) * exp(-1/2 * sum(nu.^2 ./ Q_repeated)); % (M x N)
        [Psi(1,i,:), c(1,i,:)] = max(psi,[],2); % (n_observations x n_particle (n x M))
       
    end
    outlier = mean(Psi,3) <= lambda_psi;
end
