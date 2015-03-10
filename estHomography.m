%function [ output_args ] = estHomography( input_args )
%ESTHOMOGRAPHY Summary of this function goes here
%   Detailed explanation goes here

A = imread('Lena.tif');
%A = single(A) ./ single(max(A(:)));

subplot 221;
imagesc(A); colorbar; title original;

% links oben 100,230
% rechts oben ~400, 230
% links unten 0, 512?
% rechts unten 512, 512

correctMatrix = [0.58464 -0.19531  100
                 0        0.13542  230
                 0       -0.00081    1];
			 
% x1  - 3xN set of homogeneous points
% x2  - 3xN set of homogeneous points such that x1<->x2

x1 = [  1   1  1    % links oben
      512   1  1    % rechts oben
        1 512  1    % links unten
      512 512  1]';  % rechts unten

x2 = [100 230  1    % links oben
      400 230  1    % rechts oben
        1 512  1    % links unten
      512 512  1]';  % rechts unten

% fast normalize instead of sqrt(2) from centroid
minX1 = min(x1, 1); minX2 = min(x2, 1);
maxX1 = max(x1, 1); maxX2 = max(x2, 1);

if abs(

Npts = length(x1);
A = zeros(3*Npts,9);

O = [0 0 0];
for n = 1:Npts
	X = x1(:,n)';
	u = x2(1,n); v = x2(2,n); w = x2(3,n);
	A(3*n-2,:) = [  O  -w*X  v*X];
	A(3*n-1,:) = [ w*X   O  -u*X];
	A(3*n  ,:) = [-v*X  u*X   O ];
end

[U,D,V] = svd(A, 0); % 'Economy' decomposition for speed

% Extract homography
H = reshape(V(:,9),3,3)';

% Denormalise
H = T2\H*T1;
