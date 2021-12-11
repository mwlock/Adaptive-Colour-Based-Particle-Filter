% This function is the implementation of the measurement model.
% The bearing should be in the interval [-pi,pi)
% Inputs:
%           S(t)                           4XM
%           j                              1X1
% Outputs:  
%           z_j                              2XM
function z_j = observation_model(S, j)

    global map % map including the coordinates of all landmarks | shape 2Xn for n landmarks
    global M % number of particles

    % YOUR IMPLEMENTATION
    h_theta = atan2(map(2,j) - S(2,:), map(1,j) - S(1,:)) - S(3,:);
    h = [sqrt((map(1,j) - S(1,:)).^2 + (map(2,j) - S(2,:)).^2);
         mod(h_theta + pi, 2 * pi) - pi];

    z_j = h;

end

%     global map % map including the coordinates of all landmarks | shape 2Xn for n landmarks
%     global M % number of particles
% 
%     % YOUR IMPLEMENTATION
% 
%     x = S;
% 
%     % calculate angle
%     theta = mod(x(3,:) + pi, 2*pi) - pi;
% 
%     % YOUR IMPLEMENTATION %
%     a = sqrt( (map(1,j)-x(1,:)).^2 +  (map(2,j)-x(2,:)).^2 );
%     b = atan2( map(2,j)-x(2,:), map(1,j)-x(1,:)) -theta;
%     b = mod(b + pi, 2*pi) - pi;
% 
%     h = [
%         a;
%         b
%     ];
% 
%     z_j = h;
% 
% 
% end
