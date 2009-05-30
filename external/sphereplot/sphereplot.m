function patches = sphereplot(pos,vertnum,imnum,radfun,colfun)
% SPHEREPLOT Plot a function in spherical coordinates using texturemapping
%   The sphereplot function plots a function of spherical coordinates. The
%   output is a sphere, with augmented radius according to the radfun
%   parameter, that is colored using the colfun parameter. The coloring is
%   performed using texturemapping for efficiency.
%
%   The sphere is plotted using 6 surf-objects, derived from a cube, to 
%   avoid degenerate coordinates at the "north and south poles" of the 
%   sphere.
%
%   Setting radfun = inline('1','theta','phi') yields a perfect sphere.
%
% Parameters:
%   pos: a 3-D vector specifying the center of the sphere when plotted
%   vertnum: vertex resolution
%   imnum: texture mapping resolution 
%   radfun(theta,phi): a function that gives a radius
%   colfun(theta,phi): a function that gives a texturemap color
%
% Example:    
%   % Notice the low number of vertices and the high texture accuracy.
%   rfun = inline('cos(theta*3).*cos(phi*3)*0.2 + 1','theta','phi');
%   cfun = ['cos(51*sin(phi)).*cos(51*cos(phi).*cos(theta)).*', ...
%           'cos(51*cos(phi).*sin(theta))*0.3 + 1'];
%   cfun = inline(cfun,'theta','phi');
%   h = sphereplot([0 0 0],16,500,rfun,cfun);
%   axis equal
%   set(h,'EdgeColor','black');
%   set(h,'EdgeAlpha',0.3);
%
%   
%
%
% See also SPH2CART.
%
% Author: Anders Brun
%         anders@cb.uu.se


C = [];
C(1,:,:) = [1 0 0;
            0 1 0;
            0 0 1];
C(2,:,:) = [1 0 0;
            0 1 0;
            0 0 -1];
C(3,:,:) = [1 0 0;
            0 0 1;
            0 1 0];
C(4,:,:) = [1 0 0;
            0 0 -1;
            0 1 0];
C(5,:,:) = [0 0 1;
            1 0 0;
            0 1 0];
C(6,:,:) = [0 0 -1;
            1 0 0;
            0 1 0];

patches = [];
for k = 1:6
    % generate impatch
    [uu,vv] = meshgrid((0:imnum)/imnum*2-1,(0:imnum)/imnum*2-1);
    X = uu*C(k,1,1)+vv*C(k,1,2)+C(k,1,3);
    Y = uu*C(k,2,1)+vv*C(k,2,2)+C(k,2,3);
    Z = uu*C(k,3,1)+vv*C(k,3,2)+C(k,3,3);
    D = sqrt(X.^2 + Y.^2 + Z.^2);
    X = X./D;
    Y = Y./D;
    Z = Z./D;
    [theta,phi,r] = cart2sph(X,Y,Z);
    
    im = colfun(theta,phi);

    [uu,vv] = meshgrid((0:vertnum)/vertnum*2-1,(0:vertnum)/vertnum*2-1);
    X = uu*C(k,1,1)+vv*C(k,1,2)+C(k,1,3);
    Y = uu*C(k,2,1)+vv*C(k,2,2)+C(k,2,3);
    Z = uu*C(k,3,1)+vv*C(k,3,2)+C(k,3,3);
    D = sqrt(X.^2 + Y.^2 + Z.^2);
    X = X./D;
    Y = Y./D;
    Z = Z./D;
    [theta,phi,r] = cart2sph(X,Y,Z);
    R = radfun(theta,phi);
    p = surf(pos(1)+X.*R,pos(2)+Y.*R,pos(3)+Z.*R);
    set(p,'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CData',im,...
        'CDataMapping','scaled');
    patches(end+1) = p;
    hold on;
end

% Make vertices

