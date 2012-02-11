function [radius]=arrayradius(lat,lon,azrng,ellipsoid)
%ARRAYRADIUS    Returns the radius of an array
%
%    Usage:    radius=arrayradius(lat,lon)
%              radius=arrayradius(lat,lon,azrng)
%              radius=arrayradius(lat,lon,azrng,[a f])
%
%    Description:
%     RADIUS=ARRAYRADIUS(LAT,LON) returns the maximum radius of an
%     array of positions given by latitudes LAT and longitudes LON.  LAT &
%     LON must be equal sized or scalar, are assumed to be in degrees, and
%     are based in the WGS-84 reference ellipsoid.  RADIUS is in
%     kilometers from the array's center position (see ARRAYCENTER).
%
%     RADIUS=ARRAYRADIUS(LAT,LON,AZRNG) limits the azimuth to which the
%     radius corresponds to AZRNG.  AZRNG should be [AZMIN AZMAX] in
%     degrees.  Multiple ranges may be specified.  Azimuths are from the
%     array's center position.
%
%     RADIUS=ARRAYRADIUS(LAT,LON,AZRNG,[A F]) allows specifying the
%     ellipsoid parameters A (equatorial radius in kilometers) and F
%     (flattening).  The default corresponds to WGS-84.  This is compatible
%     with output from Matlab's Mapping Toolbox function ALMANAC.
%
%    Notes:
%
%    Examples:
%     % Get the radius variation in 30deg increments:
%     arrayradius(lat,lon,[(0:30:150)' (30:30:180)'])
%
%    See also: ARRAYCENTER, ARRAYAPERTURE

%     Version History:
%        Feb.  9, 2012 - initial version (based on ARRAYAPERTURE)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2012 at 19:25 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% default azrng
if(nargin<3); azrng=[0 360]; end

% default - WGS-84 Reference Ellipsoid
if(nargin<4 || isempty(ellipsoid))
    % a=radius at equator (major axis)
    % f=flattening
    a=6378.137;
    f=1/298.257223563;
else
    % manually specify ellipsoid (will accept almanac output)
    if(isreal(ellipsoid) && numel(ellipsoid)==2 && ellipsoid(2)<1)
        a=ellipsoid(1);
        f=ellipsoid(2);
    else
        error('seizmo:arrayradius:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% size up inputs
sz{1}=size(lat); sz{2}=size(lon);
n(1)=prod(sz{1}); n(2)=prod(sz{2});

% check inputs
if(~isreal(lat) || ~isreal(lon) || ~isreal(azrng))
    error('seizmo:arrayradius:badInput',...
        'Inputs must be real valued arrays!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:arrayradius:badSize',...
        'All inputs must be equal sized or scalar!');
elseif(size(azrng,2)~=2 || ndims(azrng)~=2 || any(azrng(:,1)>=azrng(:,2)))
    error('seizmo:arrayradius:badInput',...
        'AZRNG must be a Nx2 array of [AZMIN AZMAX]!');
end

% get center position of array
[clat,clon]=arraycenter(lat,lon,[a f]);

% get station pair distance & azimuth
[d,az]=vincentyinv(clat,clon,lat,lon,[a f]);

% column vectors
d=d(:); az=az(:);

% get azrng setup
azrng=mod(azrng,360);
flip=(azrng(:,2)-azrng(:,1))<=0;
azrng(flip,2)=azrng(flip,2)+360;

% loop of ranges
nrng=size(azrng,1);
radius=nan(nrng,1);
for i=1:nrng
    in=((az-azrng(i,1))>=0 & (az-azrng(i,2))<=0) ...
        | ((az+360-azrng(i,1))>=0 & (az+360-azrng(i,2))<=0);
    if(~any(in))
        warning('seizmo:arrayradius:badAzRange',...
            'No data in azimuthal range: %d to %d!',...
            azrng(i,1),azrng(i,2));
        radius(i)=nan;
    else
        radius(i)=max(d(in));
    end
end

end
