% x0 = 0.18;
% x1 = 0.25;
% x = linspace(x0-0.1,x1+0.1)
% y = sigmoidn(x, x0, x1, 0, 1)
% plot(x,y)


function y = sigmoidn(x, x0, x1, y0, y1)

scale = 10 / (x1-x0);
shift = (x1+x0) / 2;
y = 1 ./ (1 + exp(-scale .* (x-shift)));
y(x<x0) = 0;
y(x>x1) = 1;
y = y0 + y*(y1-y0);






