function [ output_args ] = scalespace( input_args )
%SCALESPACE Summary of this function goes here
%   Detailed explanation goes here

A = imread('Lena.tif');
% rgb2gray
A = 0.2989 * A(:,:,1) + 0.5870 * A(:,:,2) + 0.1140 * A(:,:,3);
subplot(221);
imagesc(A);

[X Y] = meshgrid(-2:2, -2:2);
K = gaussian(X, Y, 1);
A2 = conv2(A, K);
subplot(222);
imagesc(A2);
colormap(gray);

function g = gaussian(x,y,t)
g = 1/(2*pi*t.^2)*exp(-(x.^2 + y.^2)/(2*t.^2));