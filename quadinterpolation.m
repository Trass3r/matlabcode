function quadinterpolation()

px = [1 4 7 2];
py = [1 1.5 5 4];
pz = [1 0 4 2];
A = [1 0 0 0;1 1 0 0;1 1 1 1;1 0 1 0];
AI = inv(A);
a = AI*px';
b = AI*py';
plotinter(px, py, pz, a, b);

function z = lerp(z1, z2, w)
	z = (1-w)*z1 + w*z2;
end

% plot surface of the quad
function []=plotinter(px, py, pz, a, b)
	x = [];
	y = [];
	z = [];
	% convert square coordinates to physical
    for l=0:0.01:1
		for m=0:0.01:1
			x(end+1) = a(1) + a(2)*l + a(3)*m + a(4)*l*m;
			y(end+1) = b(1) + b(2)*l + b(3)*m + b(4)*l*m;
			
			z00 = pz(1); z10 = pz(2);
			z01 = pz(4); z11 = pz(3);
			
			zi = lerp(lerp(z00, z10, l), lerp(z01, z11, l), m);
			z(end+1) = zi;
		end
    end
 
	plot3([px px(1)],[py py(1)], [pz pz(1)]);
	hold on;
%	plot3(x,y,z,'x');
	[XI, YI] = meshgrid(0:0.05:10);
	surf(XI, YI, griddata(x,y,z,XI,YI));
	xlabel 'x';
	ylabel 'y';
	zlabel 'z';
end

end