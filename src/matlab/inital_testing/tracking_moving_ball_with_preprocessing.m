%% Finding ball in video
%
% This file is used as an inital test script for implementing 
% "3D Trajectory Prediction of Basketball Shot Using Filtering Techniques
% and Computer Vision" project. The project is a self picked topic for implementation
% in the Appied Estimation course at% the KTH Royal Institute of Technology in 2021.
%
% Authors : 
% Matthew William Lock (mwlock@kth.se)
% Miguel Garcia Naude (magn2@kth.se)

%% To find a ball in an image

% Read in Video
video = VideoReader("freethrow_video_cut.mp4");

duration = video.Duration;
frame_rate = video.FrameRate;

% Specify min and max radius of circles
Rmin = 3;
Rmax = 25;

while hasFrame(video)
    frame = readFrame(video);
%     blurred_frame = imgaussfilt(frame, 3);
    figure(1), clf(1);

    % Image after applying filter
    [BW, maskedRGBImage] = colorMask(frame);

%     subplot(1,2,2); imshow(maskedRGBImage); title("HSV")
    imshow(maskedRGBImage);

%     
%     [centers, radii,metric] = imfindcircles(frame,[Rmin Rmax]);
% 
%     centersStrong5 = centers(1:3,:); 
%     radiiStrong5 = radii(1:3);
%     metricStrong5 = metric(1:3);
%     
%     viscircles(centersStrong5, radiiStrong5,'EdgeColor','b');

end