function [lat,lon]=fixlatlon(lat,lon)
%FIXLATLON    Returns latitudes & longitudes in reasonable ranges
%
%    Usage:    [lat,lon]=fixlatlon(lat,lon)
%
%    Description:
%     [LAT,LON]=FIXLATLON(LAT,LON) returns latitudes within the range
%     +/-90 and longitudes within +/-180 taking care to preserve the actual
%     corresponding location.  The inputs LAT and LON must be
%     real arrays of the same size, or real scalars.
%
%    Notes:
%
%    Examples:
%     % Some dumb programs may go "over the pole" in terms of latitude.
%     % FIXLATLON can fix this while handling the accompanying shift in
%     % longitude:
%     [lat,lon]=fixlatlon(100,0)
%     % returns:
%     %  lat =   80
%     %  lon = -180
%
%    See also: LATMOD, LONMOD, MOD, REM

%     Version History:
%        Feb. 16, 2010 - initial version
%        Feb. 11, 2011 - mass nargchk fix, fix see also section
%        Feb.  9, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check inputs
if(~isreal(lat) || ~isreal(lon))
    error('seizmo:fixlatlon:badInput',...
        'LAT & LON must be real arrays!');
elseif(~isscalar(lat) && ~isscalar(lon) && ~isequal(size(lat),size(lon)))
    error('seizmo:fixlatlon:badInput',...
        'LAT & LON must be scalar or equal sized arrays!');
end

% wrap position
[lat,p]=latmod(lat,90);
lon=lonmod(lon+mod(p,2)*180,360);

end
