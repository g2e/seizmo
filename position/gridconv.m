function [theta]=gridconv(clat,clon,lat,lon)
%GRIDCONV    Calculates grid convergence for orthographic projection
%
%    Usage:    theta=gridconv(clat,clon,lat,lon)
%
%    Description:
%     THETA=GRIDCONV(CLAT,CLON,LAT,LON) finds the grid convergence (that
%     is, the angle from true north to grid north where cw is positive) for
%     points specified by LAT & LON given the grid is an orthographic
%     projection (aka local tangent plane) with a tangent location given by
%     CLAT & CLON.  All inputs must be in degrees!  The output angle THETA
%     is in degrees.
%
%    Notes:
%     - Poorly written.  Please make it better.
%
%    Examples:
%     % Map of grid convergence:
%     [lat,lon]=meshgrid(-90:90,-180:180);
%     theta=gridconv(45,45,lat,lon);
%     figure; imagesc(-180:180,-90:90,theta')
%     set(gca,'ydir','normal')
%
%    See also: GEOGRAPHIC2ENU, ENU2GEOGRAPHIC

%     Version History:
%        July 22, 2010 - initial version
%        Feb.  9, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2012 at 16:50 GMT

% todo:

% check nargin
error(nargchk(4,4,nargin));

% size up inputs
sz{1}=size(clat); sz{2}=size(clon);
sz{3}=size(lat); sz{4}=size(lon);
n(1)=prod(sz{1}); n(2)=prod(sz{2});
n(3)=prod(sz{3}); n(4)=prod(sz{4});

% check inputs
if(~isreal(clat) || ~isreal(clon) || ~isreal(lat) || ~isreal(lon))
    error('seizmo:gridconv:badInput',...
        'All position inputs must be real arrays!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:gridconv:badSize',...
        'All inputs must be equal sized or scalar!');
end

% get grid convergence
% - This is terrible.  I use the difference between each point and a
%   point 1km north of them to get the local north direction.
% - ok for orthographic projection
[e1,n1]=geographic2enu(lat,lon,0,clat,clon,0);
[lat2,lon2]=vincentyfwd(lat,lon,1,0); % 1km north
[e2,n2]=geographic2enu(lat2,lon2,0,clat,clon,0);
theta=180/pi*atan2(e2-e1,n2-n1);

%theta=atand(-sind(lat).*tand(lon-clon)); % utm?

end
