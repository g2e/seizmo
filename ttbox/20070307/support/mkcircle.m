function [xp,yp,handle]=mkcircle(x,y,r,phi0,phi360,n,style);
% mkcircle..............plots a circle into the current plot
%
% plots a circle of radius r around (x,y). The circle will be approximated
% by a polygon.
%
% call: [xp,yp,handle]=mkcircle(x,y,r,phi0,phi360,n);
%       [xp,yp,handle]=mkcircle(x,y,r,phi0,phi360,n,style);
%
%           x: x of center
%           y: y of center
%           r: radius
%           phi0: smallest angle - lower angular bound
%           phi360: largest angle - higher angular bound
%                   this has to be the larger angle!
%           The angular parms have to be in radians.
%           n: number of corners for the approximating polygon
%           style: a line-style-string as in PLOT. default:'k'
%
% result: xp: x coordinates of plotted points
%         yp: y coordinates of plotted points
%         handle: handle to the created line object (perimeter only)
%
% to obtain a full circle, set phi0=0, phi360=2*pi
%
% Martin Knapmeyer, 03.06.1995, 16.02.2006

% returning of handle added 16.02.2006


if nargin<7
   style='k';
end; % if nargin<7

% x plot coordinates
xp=r*cos((1:n).*(abs(phi360-phi0)/n)+phi0)+x;

% y plot coordinates
yp=r*sin((1:n).*(abs(phi360-phi0)/n)+phi0)+y;

% plot
%hold on;
plot(x,y,style);
handle=plot(xp,yp,style);
%hold off;



