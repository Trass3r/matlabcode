x=(-5:1/32:5)';
y=gamma(x);
y=sin(x)./x;
[N, D]=ratpolyfit(x,y,3,3);
figure(1); plot(roots(N),'ob'); hold on; plot(roots(D),'xr'); grid on   

yy=polyval(N,x)./polyval(D,x);

figure(2);plot(x,y,'b', x,yy,'dr'); grid on;
axis([min(x) max(x) -0.5 1]);
%axis([min(x) max(x) -25 25]);