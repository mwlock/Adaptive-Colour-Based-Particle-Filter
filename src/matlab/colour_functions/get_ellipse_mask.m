% This file is used to generate image masks from samples
% Mask returned is a ellipsoid of specified domensions, with the center
% determined by the sample position
%
% Authors : 
% Matthew William Lock (mwlock@kth.se)
% Miguel Garcia Naude (magn2@kth.se)

function [mask] = get_ellipse_mask(x,y,Hx,Hy,xgrid,ygrid)

    mask = ( (xgrid-x).^2/Hx^2 + (ygrid-y).^2/Hy^2 )<= 1;

end