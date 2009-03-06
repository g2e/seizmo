function dmy=mkfillcircle(x,y,r,phi0,phi360,n,style,col);
% mkfillcircle..............plots a filled circle into the current plot
%
% plots a circle of radius r around (x,y). The circle will be approximated
% by a polygon and filled with the given color.
%
% call: dmy=mkcircle(x,y,r,phi0,phi360,n,style,col);
%
%           x: x of center
%           y: y of center
%           r: radius
%           phi0: smallest angle - lower angular bound
%           phi360: largest angle - higher angular bound
%                   this has to be the larger angle!
%           The angular parms have to be in radians.
%           n: number of corners for the approximating polygon
%           style: a line-style-string as in PLOT.
%           col: fill color (in colspec format)
%                if col is NaN, a white circle is produced.
%
% result: always zero
%
% to obtain a full circle, set phi0=0, phi360=2*pi
%
% Martin Knapmeyer, 05.04.2002




% x plot coordinates
xp=r*cos((0:n).*(abs(phi360-phi0)/n)+phi0)+x;

% y plot coordinates
yp=r*sin((0:n).*(abs(phi360-phi0)/n)+phi0)+y;

% add circle center to point list to get a piece of cake
xp=[x xp];
yp=[y yp];

% plot
   hold on
   plot(xp,yp,style);
   if (sum(isnan(col))~=0)|(~isempty(find(col==-1)))|(~isempty(find(abs(col)==inf))) %(isnan(col)|(col==-1)|abs(col)==-inf)
      col=[1 1 1];
   end; % if !isnan(col)
   if style=='w'
      fill(xp,yp,col,'EdgeColor','none');
   else
      fill(xp,yp,col);
   end; % if style

   hold off;

dmy=0;


