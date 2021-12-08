%% source : https://se.mathworks.com/matlabcentral/answers/437523-how-to-create-a-filled-circle

function [circles, fill_colour] = circle(x,y,r,c)
th = 0:pi/50:2*pi;
x_circle = r * cos(th) + x;
y_circle = r * sin(th) + y;
circles = plot(x_circle, y_circle,'w');
fill_colour = fill(x_circle, y_circle, c);
end