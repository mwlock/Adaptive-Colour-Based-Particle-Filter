%% Simulate Basket ball motion]
% heavily adapter from https://se.mathworks.com/matlabcentral/answers/497819-basket-shot-3d-animation

function projmotion(x0,y0,z0,v0,theta,simulation_seconds,delta_t)

g = 9.81;

angle = theta*(pi./180);
hangtime = 2*v0*sin(angle)/g;
t = hangtime;

vix = v0*cos(angle);
viy = v0*sin(angle);
viz = v0*cos(angle)*sin(angle)*t


% x = x0+vix*t;
% z = z0+viz*t;
% y = y0+viy*t-(g*t.^2)/2;
% 
% 
% maxheight = y0 + (viy)^2./(2*g)
% xheight = x0 + vix*(t/2);
% zheight = z0 + viz*t + 2

if theta >90
    error('Please select angle value of 90 degrees or less')
end;

% figure('Color', [1 1 1]);
% for k=0:t/100:t
%     x = x0 +(vix*k);
%     y = y0 +(viy*k)-0.5*g*(k^2);
%     Z = z0 +(viz*k);
% end

T = 0:t/100:t;
t_length= length(T);
X = zeros(1,t_length);
Y = zeros(1,t_length);
Z = zeros(1,t_length);
  
index = 1;
for k=0:t/100:t
    X(index) = x0 +(vix*k);
    Z(index) = z0 +(viz*k);
    Y(index) = y0 +(viy*k)-0.5*g*(k^2); 
    index = index +1;

end

x_max = max(X)+1;
y_max = max(Y)+1;
z_max = max(Z)+1;
axes_max = max(x_max,z_max);
height_max = max(y_max);


index = 1;
for k=0:t/100:t

    h = plot3(X(index),Z(index),Y(index),'.');
    xlabel('Y (meters)');
    ylabel('X (meters)');
    zlabel('Height (meters)');
    title('Basket shot animation');
    set(h,'MarkerSize',10);
    set(h,'Color','g');

    xlim([-1 axes_max]);
    ylim([-1 axes_max]);
    zlim([-1 height_max]);
    daspect([1 1 1])
    
    
    grid on
    hold on;
    pause(0.02);
   
    
    [Xs,Ys,Zs] = sphere(100);
    r = 0.2;
    Xs = Xs * r;
    Zs = Zs * r;
    Ys = Ys * r;
    
    h_new = surf(Xs+X(index),Zs+Z(index),Ys+Y(index));
    set(h_new,'FaceColor',[1 0.65 0],'FaceAlpha',0.5);
    set(h_new,'LineStyle','none')
    if index > 1
        delete(h_old);
    end
    h_old = h_new;
    index = index +1;
    

end