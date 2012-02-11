function [clat,clon]=arraycenter(lat,lon,ellipsoid)
%ARRAYCENTER    Returns the geographic center of an array
%
%    Usage:    [clat,clon]=arraycenter(lat,lon)
%              [clat,clon]=arraycenter(lat,lon,[a f])
%
%    Description:
%     [CLAT,CLON]=ARRAYCENTER(LAT,LON) returns the geographic center of an
%     array of positions given by latitudes LAT and longitudes LON.  LAT &
%     LON must be equal sized or scalar, are assumed to be in degrees, and
%     are based in the WGS-84 reference ellipsoid.  CLAT & CLON are in
%     degrees.  ALL POSITIONS ARE ASSUMED TO BE SEALEVEL!
%
%     [CLAT,CLON]=ARRAYCENTER(LAT,LON,[A F]) allows specifying the
%     ellipsoid parameters A (equatorial radius in kilometers) and F
%     (flattening).  The default corresponds to WGS-84.  This is compatible
%     with output from Matlab's Mapping Toolbox function ALMANAC.
%
%    Notes:
%     - Finds the average cartesian position and then gets the
%       corresponding latitude and longitude.
%
%    Examples:
%     % Find the geographic center of an array encompassing the North Pole:
%     lat=80+10*(rand(100)-0.5);
%     lon=360*(rand(100)-0.5);
%     [clat,clon]=arraycenter(lat,lon)
%
%    See also: GEOGRAPHIC2XYZ, XYZ2GEOGRAPHIC

%     Version History:
%        May   1, 2010 - initial version
%        July 23, 2010 - nargchk fix
%        Feb.  9, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2012 at 19:00 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% ellipsoid given or not
if(nargin==2 || isempty(ellipsoid))
    [x,y,z]=geographic2xyz(lat,lon,0);
    n=numel(x);
    [clat,clon]=xyz2geographic(sum(x(:))/n,sum(y(:))/n,sum(z(:))/n);
else
    [x,y,z]=geographic2xyz(lat,lon,0,ellipsoid);
    n=numel(x);
    [clat,clon]=xyz2geographic(...
        sum(x(:))/n,sum(y(:))/n,sum(z(:))/n,ellipsoid);
end

end
