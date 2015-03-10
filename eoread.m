fileID = fopen('t.asc');
fclose(fileID);

% gpstime pos[xyz] attitude[omega phi kappa] quat
range = 1:size(A,1);
%range = 17000:10:22000;

gps   = A(range, 1);
pos   = A(range, 2:4);
omega = A(range, 5);
phi   = A(range, 6);
kappa = A(range, 7);
quat  = A(range, 8:11);

num = length(omega);

% for i = 1:length(omega)
% end
% R = [
% 	cos(omega).*cos(kappa), sin(phi).*sin(omega).*cos(kappa) - cos(phi).*sin(kappa), cos(phi).*sin(omega).*cos(kappa) + sin(phi).*sin(kappa)
% 	cos(omega).*sin(kappa), sin(phi).*sin(omega).*sin(kappa) + cos(phi).*cos(kappa), cos(phi).*sin(omega).*sin(kappa) - sin(phi).*cos(kappa)
% 	-sin(omega),           sin(phi).*cos(omega),                                  cos(phi).*cos(omega)
% 	];

eoprange = 1:100:num;
vecs = repmat([1 1 0], length(eoprange), 1);
C = cross(quat(eoprange,2:4), vecs, 2);
C = C+C;
vecs = vecs + bsxfun(@times, quat(eoprange, 1), C) + cross(quat(eoprange,2:4), C, 2);

%quats = quaternion(quat);
%quats2 = quats(eoprange);
%vecs = repmat([1 0 0], length(eoprange), 1);
%vecs = vecs * quats2;

plot3(A(:,2), A(:,3), A(:,4));
hold on;
quiver3(pos(eoprange,1), pos(eoprange,2), pos(eoprange,3), vecs(:,1), vecs(:,2), vecs(:,3), 0.1);
xlabel 'x'; ylabel 'y'; zlabel 'z';
grid on;
text(pos(1,1), pos(1,2), pos(1,3), 'start');
%axis equal;