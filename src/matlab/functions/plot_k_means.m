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

function plot_k_means(k,L,seed,scale_factor,image_sigma,Image)

figure
for i = 1:length(k) 

    K = k(i);
    I = Image;
    original=I;
    I = imresize(I, scale_factor);
    Iback = I;
    d = 2*ceil(image_sigma*2) + 1;
    h = fspecial('gaussian', [d d], image_sigma);
    I = imfilter(I, h);
    I = double(I);
    
    tic
    [ segm, centers ] = kmeans_segm(I, K, L, seed);
    toc
    Inew = mean_segments(Iback, segm);
    I = overlay_bounds(Iback, segm);
    
    subplot(length(k),3,1+(i-1)*3)
    imshow(original)
    title(sprintf('Original image'));
    
    subplot(length(k),3,2+(i-1)*3)
    imshow(I)
    title(sprintf('Bounded image, K = %d\n',K));
    
    subplot(length(k),3,3+(i-1)*3)
    imshow(Inew)
    title(sprintf('K means, K = %d\n',K));
end