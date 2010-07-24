function [aperture]=arrayaperture(lat,lon,azirng,ellipsoid)
%ARRAYAPERTURE    Returns the aperture of an array
%
%    Usage:    aperture=arrayaperture(lat,lon)
%              aperture=arrayaperture(lat,lon,azirng)
%              aperture=arrayaperture(lat,lon,azirng,[a f])
%
%    Description: APERTURE=ARRAYAPERTURE(LAT,LON) returns the maximum
%     aperture of an array of positions given by latitudes LAT and
%     longitudes LON.    LAT & LON must be equal sized or scalar, are
%     assumed to be in degrees, and are based in the WGS-84 reference
%     ellipsoid.  APERTURE is in kilometers.
%
%     APERTURE=ARRAYAPERTURE(LAT,LON,AZIRNG) limits the azimuth to which
%     the aperture corresponds to AZIRNG.  AZIRNG should be [AZIMIN AZIMAX]
%     in degrees.  Multiple ranges may be specified.
%
%     APERTURE=ARRAYAPERTURE(LAT,LON,AZIRNG,[A F]) allows specifying the
%     ellipsoid parameters A (equatorial radius in kilometers) and F
%     (flattening).  The default corresponds to WGS-84.  This is compatible
%     with output from Matlab's Mapping Toolbox function ALMANAC.
%
%    Notes:
%
%    Examples:
%     Get the aperture variation in 30deg increments:
%      arrayaperture(lat,lon,[(0:30:150)' (30:30:180)'])
%
%    See also: ARRAYCENTER

%     Version History:
%        July 23, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 23, 2010 at 19:25 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% default azirng
if(nargin<3); azirng=[0 360]; end

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
        error('seizmo:arrayaperture:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% check inputs
if(~isreal(lat) || ~isreal(lon) || ~isreal(azirng))
    error('seizmo:arrayaperture:badInput',...
        'Inputs must be real valued arrays!');
elseif(~isequalsizeorscalar(lat,lon))
    error('seizmo:arrayaperture:badInput',...
        'LAT & LON must be equal sized or scalar!');
elseif(size(azirng,2)~=2 || ndims(azirng)~=2)
    error('seizmo:arrayaperture:badInput',...
        'AZIRNG must be a Nx2 array of [AZIMIN AZIMAX]!');
end

% expand and vectorize
[lat,lon]=expandscalars(lat,lon);
lat=lat(:); lon=lon(:);
nsta=numel(lat);

% get station pair distance & azimuth
[dist,azi]=vincentyinv(lat(:,ones(nsta,1)),lon(:,ones(nsta,1)),...
    lat(:,ones(nsta,1))',lon(:,ones(nsta,1))',[a f]);

% extract upper triangle
%dist=dist(triu(true(nsta),1));
%azi=azi(triu(true(nsta),1));

% column vectors
dist=dist(:);
azi=azi(:);

% loop of ranges
nrng=size(azirng,1);
aperture=nan(nrng,1);
for i=1:nrng
    tmp=max(dist((azi>=azirng(i,1) & azi<=azirng(i,2)) ...
        | (azi>=azirng(i,1)-360 & azi<=azirng(i,2)-360) ...
        | (azi>=azirng(i,1)+360 & azi<=azirng(i,2)+360)));
    if(isempty(tmp))
        error('seizmo:arrayaperture:badAziRange',...
            'No data in azimuthal range: %d to %d!',...
            azirng(i,1),azirng(i,2));
    end
    aperture(i)=tmp;
end

end
