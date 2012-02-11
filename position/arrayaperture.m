function [aperture]=arrayaperture(lat,lon,azrng,ellipsoid)
%ARRAYAPERTURE    Returns the aperture of an array
%
%    Usage:    aperture=arrayaperture(lat,lon)
%              aperture=arrayaperture(lat,lon,azrng)
%              aperture=arrayaperture(lat,lon,azrng,[a f])
%
%    Description:
%     APERTURE=ARRAYAPERTURE(LAT,LON) returns the maximum aperture of an
%     array of positions given by latitudes LAT and longitudes LON.  LAT &
%     LON must be equal sized or scalar, are assumed to be in degrees, and
%     are based in the WGS-84 reference ellipsoid.  APERTURE is in
%     kilometers.  Unlike ARRAYRADIUS this does not use the array's center
%     for the measurements.
%
%     APERTURE=ARRAYAPERTURE(LAT,LON,AZRNG) limits the azimuth to which
%     the aperture corresponds to AZRNG.  AZRNG should be [AZIMIN AZIMAX]
%     in degrees.  Multiple ranges may be specified.
%
%     APERTURE=ARRAYAPERTURE(LAT,LON,AZRNG,[A F]) allows specifying the
%     ellipsoid parameters A (equatorial radius in kilometers) and F
%     (flattening).  The default corresponds to WGS-84.  This is compatible
%     with output from Matlab's Mapping Toolbox function ALMANAC.
%
%    Notes:
%     - This algorithm is not good for azimuthal variations of pole
%       encompassing arrays.  ARRAYRADIUS may work better in those cases.
%
%    Examples:
%     % Get the aperture variation in 30deg increments:
%     arrayaperture(lat,lon,[(0:30:150)' (30:30:180)'])
%
%    See also: ARRAYCENTER, ARRAYRADIUS

%     Version History:
%        July 23, 2010 - initial version
%        Feb.  9, 2012 - doc update, warn/nan for azrng w/o data
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2012 at 19:25 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% default azrng
if(nargin<3); azrng=[0 360]; end

% default ellipsoid - WGS-84 Reference Ellipsoid
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
        error('seizmo:arrayaperture:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% size up inputs
sz{1}=size(lat); sz{2}=size(lon);
n(1)=prod(sz{1}); n(2)=prod(sz{2});

% check inputs
if(~isreal(lat) || ~isreal(lon) || ~isreal(azrng))
    error('seizmo:arrayaperture:badInput',...
        'Inputs must be real valued arrays!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:arrayaperture:badSize',...
        'LAT & LON must be equal sized or scalar!');
elseif(size(azrng,2)~=2 || ndims(azrng)~=2)
    error('seizmo:arrayaperture:badInput',...
        'AZRNG must be a Nx2 array of [AZMIN AZMAX]!');
end

% vectorize & expand
lat=lat(:); lon=lon(:);
if(n(1)==1); lat=lat(ones(n(2),1)); end
if(n(2)==1); lon=lon(ones(n(1),1)); end
nsta=numel(lat);

% get station pair distance & azimuth
[d,az]=vincentyinv(lat(:,ones(nsta,1)),lon(:,ones(nsta,1)),...
    lat(:,ones(nsta,1))',lon(:,ones(nsta,1))',[a f]);

% column vectors
d=d(:); az=az(:);

% get azrng setup
azrng=mod(azrng,360);
flip=(azrng(:,2)-azrng(:,1))<=0;
azrng(flip,2)=azrng(flip,2)+360;

% loop of ranges
nrng=size(azrng,1);
aperture=nan(nrng,1);
for i=1:nrng
    in=((az-azrng(i,1))>=0 & (az-azrng(i,2))<=0) ...
        | ((az+360-azrng(i,1))>=0 & (az+360-azrng(i,2))<=0);
    if(~any(in))
        warning('seizmo:arrayaperture:badAzRange',...
            'No data in azimuthal range: %d to %d!',...
            azrng(i,1),azrng(i,2));
        aperture(i)=nan;
    else
        aperture(i)=max(d(in));
    end
end

end
