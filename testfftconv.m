%function [ output_args ] = testfftconv( input_args )
%TESTFFTCONV Summary of this function goes here
%   Detailed explanation goes here
[X, Y] = meshgrid(-5:0.05:5);
K = exp(-(X.^2 + Y.^2)/2) / (sqrt(2*pi));
Z = fftshift(fft2(K) / (size(K,1) * size(K,2)));
%surf(real(Z)); hold on; surf(imag(Z) + 5); hold off;
subplot 221;
imagesc(real(Z)); colorbar; axis image; title real;
subplot 222;
imagesc(imag(Z)); colorbar; axis image; title imag;
subplot 223;
imagesc(abs(Z)); colorbar; axis image; title abs;
subplot 224;
imagesc(angle(Z)); colorbar; axis image; title angle;

pause(0.001); % jframe nastiness
jframe = get(handle(gcf), 'JavaFrame');
jframe.setMaximized(true); % maximize figure
