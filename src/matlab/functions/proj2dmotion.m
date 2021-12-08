%% Simulate Basket ball motion]
% heavily adapter from https://se.mathworks.com/matlabcentral/answers/497819-basket-shot-3d-animation

function [X,Y] = proj2dmotion (x0,y0,v0,theta,seconds,delta_t)

% Constants and environment
e = 0.8;                                        % co-effecient of restition    
g = 9.81;                                       % gravity
t = 0:delta_t:seconds;                          % simulation time
angle = theta*(pi./180);
r = 0.12;                                        % radius of the ball                              

if x0-r<0
    x0 = r;
end

if y0 - r < 0
    y0 = r;
end

vix = v0*cos(angle);
viy_0 = v0*sin(angle);

if theta >90
    error('Please select angle value of 90 degrees or less')
end

% Declare arrays
t_length = length(t);
X = zeros(1,t_length); X(1) = x0;
Y = zeros(1,t_length); Y(1) = y0;
VY = zeros(1,t_length); VY(1) = viy_0;

% Simulate motion
for k=2:t_length
    
    % update y velocity
    VY(k) = VY(k-1)-g*delta_t;

    X(k) = X(k-1)+delta_t*vix;
    Y(k) = Y(k-1)+delta_t*VY(k-1); 

    % check for bounce
    if (Y(k) - r) < 0 && VY(k-1) < 0
        VY(k) = abs(VY(k))*e;
    end

end

% Determine maximum displacement
x_max = max(X)+1;
y_max = max(Y)+1;

if x_max > y_max*2
    y_max = x_max/2;
end

y_max = 2.5;

for k=1:t_length
    
    % Plot center position of ball
    h = plot(X(k),Y(k),'.');
    xlabel('X (meters)');
    ylabel('Y - Height (meters)');
    title('Basket shot animation');
    set(h,'MarkerSize',10);
    set(h,'Color','g');
    
    % Axes limits
    xlim([-1 x_max]);
    ylim([-1 y_max]);
    daspect([1 1 1]);
    
    % Plot ball
    grid on
    hold on;
    pause(1/60);
     
    [h_new,fill_new] = circle(X(k), Y(k), r, [0.8500 0.3250 0.0980]);
    if k > 1
        % delete previous circle
        delete(h_old);
        delete(fill_old);
    end
    h_old = h_new;  
    fill_old = fill_new;

end