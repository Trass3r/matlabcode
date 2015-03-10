function calibpatterntest

function y = f(x)
	y = 0.5 * sin(x);
end

T = 2*pi;
% TODO: what if T is not pi??

function c = colorAt(x,y)
	c = 1; % white by default
	if y > f(x) + 0.5*T
		% upper half

		if x > (-f(y) + 0.5*T) && x < (f(y) + T)
			c = 0;
		end
		
		if y > -f(x) + T
			c = 0;
		end
 	else
		if x < -f(y) + 0.5*T
			if x > f(y)
				c = 0;
			end
 		else
 			if y < -f(x)
 				c = 0;
 			end
 		end
 	end
end

N = 2000; % samples per T
Lx = 4; % how many multiples of T
Ly = 4*sqrt(2);

[X,Y] = meshgrid(1:N);
%A = zeros(N, N);
A = arrayfun(@colorAt, (X-1)/N * T, (Y-1)/N * T);

% for j = 1:N
% 	for i = 1:N
% 		x = (i-1)/N * T;
% 		y = (j-1)/N * T;
% 		A(j, i) = colorAt(x, y);
% 	end
% end
A = repmat(A, Lx, Ly);

%imagesc(0:L/N:L, 0:L/N:L, A);
imagesc(A, [0 1]);
axis xy equal
colormap(gray);
xlabel 'x'
ylabel 'y'
colorbar;

imwrite(A, 'test.png')

end