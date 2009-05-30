function [p] = sphereplot(xc,yc,zc,vertnum,imnum,radscale,radfun,colfun,coordsys)

% Plot a function defined on x,y,z that belongs to the unit sphere.
% xc,yc,zc: center of the sphere when plotted
% vertnum: geometric resolution of sphere
% imnum: texturemapping resolution of sphere
% radscale: scales radius of sphere
% radfun(X,Y,Z): function that gives a radius
% colfun(X,Y,Z): function that gives a texturemap

% Make bitmap

  
if nargin<9
  coordsys = 'cartesian';
end


if coordsys == 'spherical'
  [X,Y,Z] = sphere(imnum);
  [theta,phi,r] = cart2sph(X,Y,Z);
  im = colfun(theta,phi,r);
  [X,Y,Z] = sphere(vertnum);
  [theta,phi,r] = cart2sph(X,Y,Z);
  R = radfun(theta,phi,r) * radscale;
else % cartesian ...
  [X,Y,Z] = sphere(imnum);
  im = colfun(X,Y,Z);
  [X,Y,Z] = sphere(vertnum);
  R = radfun(X,Y,Z) * radscale;
end

% Make vertices

p = surf(xc+X.*R,yc+Y.*R,zc+Z.*R);

h = get(p,'VertexNormals');
h(:,1,:) = (h(:,1,:) + h(:,end,:))/2;
h(:,end,:) = h(:,1,:);

set(p,'FaceColor','texturemap',...
      'EdgeColor','none',...
      'CData',im,...
      'CDataMapping','scaled',...
      'VertexNormals',h);
