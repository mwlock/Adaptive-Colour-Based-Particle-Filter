function [BW,maskedRGBImage] = colorMask(image)
    % Convert RGB Image to HSV Image
    image_hsv = rgb2hsv(image);

    % Thresholds for HUE (Percentage of where it lies on color wheel)

    hue_min = 15/360;
    hue_max = 35/360;

    % Thresholds for Saturation
    sat_min = 0.5;
    sat_max = 1;

    % Thresholds for Value
    val_min = 0.5;
    val_max = 1;

    % Create mask based on chosen histogram thresholds
    BW = ( (image_hsv(:,:,1) >= hue_min) | (image_hsv(:,:,1) <= hue_max) ) & (image_hsv(:,:,2) >= sat_min ) & (image_hsv(:,:,2) <= sat_max) & (image(:,:,3) >= val_min ) & (image_hsv(:,:,3) <= val_max);
    
    % Initialize output masked image based on input image.
    maskedRGBImage = image;

    % Set background pixels where BW is false to zero.
    maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

end