function [x,y,z]=geocentric2xyz(lat,lon,radius,r)
%GEOCENTRIC2XYZ    Converts coordinates from geocentric to cartesian
%
%    Usage:    [x,y,z]=geocentric2xyz(lat,lon,radius)
%              [x,y,z]=geocentric2xyz(lat,lon,depth,r)
%
%    Description: [X,Y,Z]=GEOCENTRIC2XYZ(LAT,LON,RADIUS) converts
%     coordinates in geocentric latitude, longitude, radius to 
%     Earth-centered, Earth-Fixed (ECEF).  LAT and LON are in degrees.  
%     X, Y, Z will match the units of RADIUS.
%
%     GEOCENTRIC2XYZ(LAT,LON,DEPTH,R) allows specifying the radius R of
%     the sphere and DEPTH rather than RADIUS.  X, Y, Z will match the
%     units of DEPTH and R (must have the same units).
%
%    Notes:
%     - the ECEF coordinate system has the X axis passing through the
%       equator at the prime meridian, the Z axis through the north pole
%       and the Y axis through the equator at 90 degrees longitude.
%
%    Tested on: Matlab r2007b
%
%    Examples:
%     Find out how far a position is from the equatorial plane (z):
%      [x,y,z]=geocentric2xyz(lat,lon,depth,r)
%
%    See also: xyz2geocentric, xyz2geodetic, geodetic2xyz

%     Version History:
%        Oct. 14, 2008 - initial version
%        Nov. 10, 2008 - scalar expansion, doc update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 21:45 GMT

% todo:

% require 3 or 4 inputs
msg=nargchk(3,4,nargin);
if(~isempty(msg)); error(msg); end

% size up inputs
sx=size(lat); sy=size(lon); sz=size(radius);
nx=prod(sx); ny=prod(sy); nz=prod(sz);

% basic check inputs
if(~isnumeric(lat) || ~isnumeric(lon) || ~isnumeric(radius))
    error('seizmo:geocentric2xyz:nonNumeric','All inputs must be numeric!');
elseif(any([nx ny nz]==0))
    error('seizmo:geocentric2xyz:unpairedCoord',...
        'Coordinate inputs must be nonempty arrays!');
elseif((~isequal(sx,sy) && all([nx ny]~=1)) ||...
       (~isequal(sx,sz) && all([nx nz]~=1)) ||...
       (~isequal(sz,sy) && all([nz ny]~=1)))
    error('seizmo:geocentric2xyz:unpairedCoord',...
        'Coordinate inputs must be scalar or equal sized arrays!');
end

% expand scalars
if(all([nx ny nz]==1))
    % do nothing
elseif(all([nx ny]==1))
    lat=repmat(lat,sz); lon=repmat(lon,sz);
elseif(all([nx nz]==1))
    lat=repmat(lat,sy); radius=repmat(radius,sy);
elseif(all([ny nz]==1))
    lon=repmat(lon,sx); radius=repmat(radius,sx);
elseif(nx==1)
    lat=repmat(lat,sz);
elseif(ny==1)
    lon=repmat(lon,sz);
elseif(nz==1)
    radius=repmat(radius,sy);
end

% check input (converts depth to radius)
if(nargin==4)
    if(~isnumeric(r) || ~isscalar(r))
        error('seizmo:geocentric2xyz:badR',...
        'Sphere radius input must be numerical scalar!');
    end
    radius=r-radius; 
end

% convert to xyz
z=radius.*sind(lat);
x=radius.*cosd(lon).*cosd(lat);
y=radius.*sind(lon).*cosd(lat);

end
