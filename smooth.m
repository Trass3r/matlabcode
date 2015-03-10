function Y = smooth(X)
%SMOOTH Summary of this function goes here
%   Detailed explanation goes here
I = -2:2;
sigma = 2;
gaussian = exp(-I.^2 / (2*sigma^2)) / (sigma * sqrt(2*pi));
Y = conv(X, gaussian);
Y = Y(3:end-2);%, 'same');