% An example, plotting spherical harmonics of order (5,2) with various 
% resolutions of the sphere grid but with relatively high resolution of the
% texture map.
% 
% 
c  = inline('real(spharm(theta,phi,5,2))','theta','phi');
%c  = inline('cos(phi).^5 + cos(theta).^7','theta','phi');

r1 = inline('1','theta','phi');
r2 = inline('5+real(spharm(theta,phi,5,2))','theta','phi');
r3 = inline('abs(real(spharm(theta,phi,5,2)))','theta','phi');

figure
subplot(3,3,1);
h = sphereplot([0 0 0],20,100,r1,c);
axis equal; axis off; zoom(1.8);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.1);

subplot(3,3,2);
h = sphereplot([0 0 0],20,100,r2,c);
axis equal; axis off;  zoom(1.8);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.1);

subplot(3,3,3);
h = sphereplot([0 0 0],20,100,r3,c);
axis equal; axis off;  zoom(1.8);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.1);

subplot(3,3,4);
h = sphereplot([0 0 0],10,100,r1,c);
axis equal; axis off; zoom(1.8);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.1);

subplot(3,3,5);
h = sphereplot([0 0 0],10,100,r2,c);
axis equal; axis off;  zoom(1.8);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.1);

subplot(3,3,6);
h = sphereplot([0 0 0],10,100,r3,c);
axis equal; axis off;  zoom(1.8);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.1);

subplot(3,3,7);
h = sphereplot([0 0 0],5,100,r1,c);
axis equal; axis off; zoom(1.8);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.1);

subplot(3,3,8);
h = sphereplot([0 0 0],5,100,r2,c);
axis equal; axis off;  zoom(1.8);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.1);

subplot(3,3,9);
h = sphereplot([0 0 0],5,100,r3,c);
axis equal; axis off;  zoom(1.8);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.1);



