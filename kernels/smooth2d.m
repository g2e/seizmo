function [sp,sa]=smooth2d(x,y,p,a,lambda)
% smooth a 2d sensitivity kernel w/ truncation

% create 2D gaussian
% assumes x,y are equally evenly spaced
dx=x(2)-x(1);
n=ceil(3*lambda/dx);
g=exp(-(abs(-n:n)*dx/lambda).^2);
gtot=sum(sum(g'*g));

% get truncation surface
t=conv2(g,g,ones(size(p)),'same');

% now get smoothed kernel
sp=conv2(g,g,p,'same')./gtot;
sa=conv2(g,g,a,'same')./gtot;

figure;
surface(x,y,sp);
colormap(jet(1024));
shading interp;
hold on
title('Smoothed Phase Kernel')
figure;
surface(x,y,sa);
colormap(jet(1024));
shading interp;
hold on
title('Smoothed Amplitude Kernel')

end