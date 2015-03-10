X = -1:0.05:3;
p0 = 50; p1 = 100; p2 = 250; p3 = 115;
%p0 = 5; p1 = 10; p2 = 15; p3 = 20;

a = p3 - p2 + p1 - p0;
b = p0 - p1 - a;
c = p2 - p0;
d = p1;

Y = a * X.^3 + b * X.^2 + c * X + d;
Ys = 3*a * X.^2 + 2*b * X + c;

plot([-1 0 1 2], [p0 p1 p2 p3], 'or',  X, Y,  X, Ys, 'g');
%ylim([p0, p3]);
hold on;

% plot slope
plot(X, Ys(X == 0)*X + p1, 'r-');

a = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
b = p0 - 2.5 * p1 + 2 * p2 - 0.5 * p3;
c = -0.5 * p0 + 0.5 * p2;
d = p1;

Y = a * X.^3 + b * X.^2 + c * X + d;
Ys = 3*a * X.^2 + 2*b * X + c;

plot(X, Y, '.-',  X, Ys, '.-g');
