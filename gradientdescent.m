function gradientdescent()

spacing = 0.05;

[X,Y] = meshgrid(0:spacing:2*pi);
Z = sin(X.*Y) + cos(Y) + sin(X);
%surf(X,Y,Z, 'EdgeAlpha', 0.99);
%axis([-3 3 -3 3 -10 5])

[gx,gy] = gradient(Z, spacing, spacing);

subplot(1,2,1);
%imagesc(Z); colorbar;
contour(X,Y,Z);
hold on
quiver(X, Y, gx, gy);

subplot(1,2,2);
imagesc(abs(gx) + abs(gy)); colorbar;
hold on;
quiver(gx, gy);