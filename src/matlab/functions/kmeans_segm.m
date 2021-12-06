%% Kmeans Segm
%
% This file is used as an inital test script for implementing 
% "3D Trajectory Prediction of Basketball Shot Using Filtering Techniques
% and Computer Vision" project. The project is a self picked topic for implementation
% in the Appied Estimation course at% the KTH Royal Institute of Technology in 2021.
%
% Authors : 
% Matthew William Lock (mwlock@kth.se)
% Miguel Garcia Naude (magn2@kth.se)

function [ segmentation, centers ] = kmeans_segm(image, K, L, seed)

    % K means pseduo code ==================================================
    % Let X be a set of pixels and V be a set of K cluster centers in 3D (R,G,B).
    % Randomly initialize the K cluster centers
    % Compute all distances between pixels and cluster centers
    % Iterate L times
    %   Assign each pixel to the cluster center for which the distance is minimum
    %   Recompute each cluster center by taking the mean of all pixels assigned to it
    %   Recompute all distances between pixels and cluster centers
    
    % Reshape image
    width = size(image,1);
    height = size(image,2);
    Ivec = reshape(image, width*height, 3);

    % Seed random number generator
    rng(seed);
    
    % generate K random inital clusters
    % cluster_indexes = randi([1,length(Ivec)],1,3); % potentially repeating pixels

    N = length(Ivec);
    
    cluster_pixels = [];
    while size(unique(cluster_pixels,'rows'),1) ~=K
        cluster_indexes = randperm(length(Ivec),K); % the same index is not picked twice
        cluster_pixels = Ivec(cluster_indexes,:);
    end

    for i = 1:L
    
        % get distances for each pixel to each cluster
        distances = pdist2(Ivec , cluster_pixels);
    
        % Get index of cluster each pixel is closest to
        [~, indexes] = min(distances,[],2);

        % reassign clusters
        for k = 1:K
                indexes_of_pixels_belonging_to_cluster = nonzeros((indexes==k).*linspace(1,N,N)');
                pixels = Ivec(indexes_of_pixels_belonging_to_cluster,:);
                cluster_pixels(k,:) = round(mean(pixels));
        end

    end
    
    segmentation = indexes;
    segmentation = reshape(segmentation,[width,height]);
    centers = cluster_pixels;  



    
