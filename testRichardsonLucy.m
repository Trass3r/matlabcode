function [ output_args ] = testRichardsonLucy( input_args )
%TESTRICHARDSON-LUCY Summary of this function goes here
%   Detailed explanation goes here

A = imread('Lena.tif');
A = single(A) ./ single(max(A(:)));

sigma = 2;
%range = -3:0.1:3;
range = -3:3;
[X,Y] = meshgrid(range);

psf = 1/(2*pi*sigma^2) * exp(-(X.^2 + Y.^2) / (2 * sigma^2));

surf(X, Y, psf);

subplot 231;
imagesc(A(:,:,1)); colormap(gray); colorbar; title original;
subplot 232;
B = cat(3, conv2(A(:,:,1), psf, 'same'), conv2(A(:,:,2), psf, 'same'), conv2(A(:,:,3), psf, 'same'));

B = B(:,:,1);
imagesc(B); colormap(gray); colorbar; title blurred;

estm = 0.5*ones(size(B));
psf_hat = psf(end:-1:1, end:-1:1);
for i = 1:50
	estmconv = conv2(estm, psf, 'same');
	relblur  = B ./ estmconv;
	mean(abs(1 - relblur(:)))
	if (mean(abs(1 - relblur(:))) < 0.001)
		disp('Iteration stopped at ', i);
	end
	errorestm = conv2(relblur, psf_hat, 'same');

	subplot 234; imagesc(estmconv); colorbar; title estmconv;
	subplot 235; imagesc(relblur, [0.9 1.25]); colorbar; title relblur;
	subplot 236; imagesc(errorestm, [0.75 1]); colorbar; title errorestm;
	
	subplot 233; imagesc(estm, [0 1.02]); colorbar; title estm; colormap(gray);
	drawnow;

	estm = estm .* errorestm;
end

end