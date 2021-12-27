% Demo to have the user freehand draw an irregular shape over a gray scale image.
% Then it creates new images:
% (1) where the drawn region is all white inside the region and untouched outside the region,
% (2) where the drawn region is all black inside the region and untouched outside the region,
% (3) where the drawn region is untouched inside the region and all black outside the region.
% It also (4) calculates the mean intensity value and standard deviation of the image within that shape,
% (5) calculates the perimeter, centroid, and center of mass (weighted centroid), and
% (6) crops the drawn region to a new, smaller separate image.

% Change the current folder to the folder of this m-file.
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.
fontSize = 16;

fullFileName = 'horse_v1.mp4';
v = VideoReader(fullFileName);

grayImage = read(v,15);
imshow(grayImage, []);
axis on;
title('Original Grayscale Image', 'FontSize', fontSize);

% Ask user to draw freehand mask.
message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
uiwait(msgbox(message));
hFH = imfreehand(); % Actual line of code to do the drawing.
% Create a binary image ("mask") from the ROI object.
binaryImage = hFH.createMask();
xy = hFH.getPosition;

% Now make it smaller so we can show more images.
subplot(2, 4, 1);
imshow(grayImage, []);
axis on;
drawnow;
title('Original gray scale image', 'FontSize', fontSize);

% Display the freehand mask.
subplot(2, 4, 2);
imshow(binaryImage);
axis on;
title('Binary mask of the region', 'FontSize', fontSize);
save('binaryImage.mat','binaryImage');

% Get the dimensions of the image.  
% numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(grayImage);
if numberOfColorBands > 1
	% It's not really gray scale like we expected - it's color.
	% Convert it to gray scale by taking only the green channel.
	grayImage = grayImage(:, :, 2); % Take green channel.
end

% Label the binary image and computer the centroid and center of mass.
labeledImage = bwlabel(binaryImage);
measurements = regionprops(binaryImage, grayImage, ...
    'area', 'Centroid', 'WeightedCentroid', 'Perimeter');
area = measurements.Area
centroid = measurements.Centroid
centerOfMass = measurements.WeightedCentroid
perimeter = measurements.Perimeter

% Calculate the area, in pixels, that they drew.
numberOfPixels1 = sum(binaryImage(:))
% Another way to calculate it that takes fractional pixels into account.
numberOfPixels2 = bwarea(binaryImage)

% Get coordinates of the boundary of the freehand drawn region.
structBoundaries = bwboundaries(binaryImage);
xy=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
x = xy(:, 2); % Columns.
y = xy(:, 1); % Rows.
subplot(2, 4, 1); % Plot over original image.
hold on; % Don't blow away the image.
plot(x, y, 'LineWidth', 2);
drawnow; % Force it to draw immediately.

% Burn region as white into image by setting it to 255 wherever the mask is true.
burnedImage = grayImage;
burnedImage(binaryImage) = 255;
% Display the image with the mask "burned in."
subplot(2, 4, 3);
imshow(burnedImage);
axis on;
caption = sprintf('Masked white inside region');
title(caption, 'FontSize', fontSize);

% Burn region as black into image by setting it to 255 wherever the mask is true.
burnedImage = grayImage;
burnedImage(binaryImage) = 0;
% Display the image with the mask "burned in."
subplot(2, 4, 4);
imshow(burnedImage);
axis on;
caption = sprintf('Masked black inside region');
title(caption, 'FontSize', fontSize);

% Mask the image white outside the mask, and display it.
% Will keep only the part of the image that's inside the mask, white outside mask.
whiteMaskedImage = grayImage;
whiteMaskedImage(~binaryImage) = 255;
subplot(2, 4, 5);
imshow(whiteMaskedImage);
axis on;
title('Masked white outside region', 'FontSize', fontSize);

% Mask the image outside the mask, and display it.
% Will keep only the part of the image that's inside the mask, zero outside mask.
blackMaskedImage = grayImage;
blackMaskedImage(~binaryImage) = 0;
subplot(2, 4, 6);
imshow(blackMaskedImage);
axis on;
title('Masked black outside region', 'FontSize', fontSize);
% Calculate the mean
meanGL = mean(blackMaskedImage(binaryImage));
sdGL = std(double(blackMaskedImage(binaryImage)));

% Put up crosses at the centriod and center of mass
hold on;
plot(centroid(1), centroid(2), 'r+', 'MarkerSize', 30, 'LineWidth', 2);
plot(centerOfMass(1), centerOfMass(2), 'g+', 'MarkerSize', 20, 'LineWidth', 2);

% Now crop the image.
leftColumn = min(x);
rightColumn = max(x);
topLine = min(y);
bottomLine = max(y);
width = rightColumn - leftColumn + 1;
height = bottomLine - topLine + 1;
croppedImage = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]);
% Display cropped image.
subplot(2, 4, 7:8);
imshow(croppedImage);
axis on;
title('Cropped image', 'FontSize', fontSize);

% Put up crosses at the centriod and center of mass
hold on;
plot(centroid(1)-leftColumn, centroid(2)-topLine, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
plot(centerOfMass(1)-leftColumn, centerOfMass(2)-topLine, 'g+', 'MarkerSize', 20, 'LineWidth', 2);

% Report results.
message = sprintf('Mean value within drawn area = %.3f\nStandard deviation within drawn area = %.3f\nNumber of pixels = %d\nArea in pixels = %.2f\nperimeter = %.2f\nCentroid at (x,y) = (%.1f, %.1f)\nCenter of Mass at (x,y) = (%.1f, %.1f)\nRed crosshairs at centroid.\nGreen crosshairs at center of mass.', ...
meanGL, sdGL, numberOfPixels1, numberOfPixels2, perimeter, ...
centroid(1), centroid(2), centerOfMass(1), centerOfMass(2));
msgbox(message);

