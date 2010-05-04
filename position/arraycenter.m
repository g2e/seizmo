function [clat,clon]=arraycenter(lat,lon,ellipsoid)
%ARRAYCENTER    Returns the geographic center of an array of positions
%
%    Usage:    [clat,clon]=arraycenter(lat,lon)
%              [clat,clon]=arraycenter(lat,lon,[a f])
%
%    Description: [CLAT,CLON]=ARRAYCENTER(LAT,LON) returns the geographic
%     center of an array of positions given by latitudes LAT and longitudes
%     LON.  LAT & LON must be equal sized or scalar, are assumed to be in
%     degrees, and are based in the WGS-84 reference ellipsoid.  CLAT &
%     CLON are in degrees.  ALL POSITIONS ARE ASSUMED TO BE AT THE SURFACE.
%
%     [CLAT,CLON]=ARRAYCENTER(LAT,LON,[A F]) allows specifying the
%     ellipsoid parameters A (equatorial radius in kilometers) and F
%     (flattening).  The default corresponds to WGS-84.  This is compatible
%     with output from Matlab's Mapping Toolbox function ALMANAC.
%
%    Notes:
%     - finds the average cartesian position and then gets the
%       corresponding latitude and longitude
%
%    Examples:
%     Find the geographic center of an array encompassing the North Pole:
%      lat=80+10*(rand(100)-0.5);
%      lon=360*(rand(100)-0.5);
%      [clat,clon]=arraycenter(lat,lon)
%
%    See also: GEOGRAPHIC2XYZ, XYZ2GEOGRAPHIC

%     Version History:
%        May   1, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   1, 2010 at 10:00 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
if(~isempty(msg)); error(msg); end

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
