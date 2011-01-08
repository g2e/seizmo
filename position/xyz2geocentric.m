function [lat,lon,radius]=xyz2geocentric(x,y,z,r)
%XYZ2GEOCENTRIC    Converts coordinates from cartesian to geocentric
%
%    Usage:    [lat,lon,radius]=xyz2geocentric(x,y,z)
%              [lat,lon,depth]=xyz2geocentric(x,y,z,r)
%
%    Description: [LAT,LON,RADIUS]=XYZ2GEOCENTRIC(X,Y,Z) converts arrays of
%     coordinates in Earth-centered, Earth-Fixed (ECEF) to geocentric
%     latitude, longitude, radius.  LAT and LON are in degrees.  X, Y and Z
%     must have the same units (so RADIUS will be in those units) and must
%     be same size arrays or scalars.
%
%     [LAT,LON,DEPTH]=XYZ2GEOCENTRIC(X,Y,Z,R) allows specifying the radius
%     R of the sphere so depth is returned rather than radius.  In this 
%     case, the units of X, Y, Z must match those of R (so DEPTH will be in
%     those units).  R must be a scalar.
%
%    Notes:
%     - the ECEF coordinate system has the X axis passing through the
%       equator at the prime meridian, the Z axis through the north pole
%       and the Y axis through the equator at 90 degrees longitude.
%
%    Examples:
%     Find out how far an xyz position is from the earth's center:
%      [lat,lon,radius]=xyz2geocentric(x,y,z)
%
%    See also: GEOCENTRIC2XYZ, XYZ2GEOGRAPHIC, GEOGRAPHIC2XYZ,
%              GEOGRAPHIC2GEOCENTRIC, GEOCENTRIC2GEOGRAPHIC

%     Version History:
%        Oct. 14, 2008 - initial version
%        Oct. 26, 2008 - scalar expansion, doc and comment update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Sep.  5, 2009 - minor doc update
%        Nov. 13, 2009 - name change: geodetic to geographic
%        Jan.  4, 2011 - minor code improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan.  4, 2011 at 20:15 GMT

% todo:

% require 3 or 4 inputs
error(nargchk(3,4,nargin));

% size up inputs
sx=size(x); sy=size(y); sz=size(z);
nx=prod(sx); ny=prod(sy); nz=prod(sz);

% basic check inputs
if(~isreal(x) || ~isreal(y) || ~isreal(z))
    error('seizmo:xyz2geocentric:nonReals','All inputs must be reals!');
elseif(any([nx ny nz]==0))
    error('seizmo:xyz2geocentric:unpairedCoord',...
        'Coordinate inputs must be nonempty arrays!');
elseif((~isequal(sx,sy) && all([nx ny]~=1)) ||...
       (~isequal(sx,sz) && all([nx nz]~=1)) ||...
       (~isequal(sz,sy) && all([nz ny]~=1)))
    error('seizmo:xyz2geocentric:unpairedCoord',...
        'Coordinate inputs must be scalar or equal sized arrays!');
end

% expand scalars
if(all([nx ny nz]==1))
    % do nothing
elseif(all([nx ny]==1))
    x=repmat(x,sz); y=repmat(y,sz);
elseif(all([nx nz]==1))
    x=repmat(x,sy); z=repmat(z,sy);
elseif(all([ny nz]==1))
    y=repmat(y,sx); z=repmat(z,sx);
elseif(nx==1)
    x=repmat(x,sz);
elseif(ny==1)
    y=repmat(y,sz);
elseif(nz==1)
    z=repmat(z,sy);
end

% convert to geocentric
radius=sqrt(x.^2+y.^2+z.^2);
lon=atan2(y,x).*(180/pi);
lat=asind(z./radius);

% check input (converts radius to depth)
if(nargin==4)
    if(~isreal(r) || ~isscalar(r))
        error('seizmo:xyz2geocentric:badR',...
        'R input must be a real-valued scalar!');
    end
    radius=r-radius; 
end

end
