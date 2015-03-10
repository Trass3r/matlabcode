function [ output_args ] = Kalman( input_args )
%KALMAN Summary of this function goes here
%   Detailed explanation goes here

close all;

x = [0 0 0]'; % x y phi
x_gt = x;

C = [0 0; 0 0];

figure;
%daspect([1,1,1]);
axis equal;
trajectory = x;
dt = 1; % secs
v = 5; % m/s
for i = 1:500
	hold off;
%	plotPos(x_gt(1), x_gt(2), x_gt(3));
	hold on;
	
%	phi = (rand-0.5) * 0.2;
%	x_gt = x_gt + [cos(phi) * v*dt; sin(phi) * v*dt; phi];
	dphi = (rand-0.5) * 0.3;
	phi = x_gt(3);
	x_gt = x_gt + [cos(phi) * v*dt; sin(phi) * v*dt; dphi];

	trajectory(:, end) = x_gt;
	line(trajectory(1, :), trajectory(2, :),'LineStyle','-');
	drawnow;
%	x_p = 
%	measurement = 
end

function h = circle2(x,y,r)
d = r*2;
px = x-r;
py = y-r;
h = rectangle('Position',[px py d d],'Curvature',[1,1]);
%daspect([1,1,1])

function h = plotPos(x,y,phi)
radius = 0.20; % meters
h = circle2(x, y, radius);
line([x x+cos(phi)*radius], [y y+sin(phi)*radius]);