function [gcarc,az,baz]=sphericalinv(evla,evlo,stla,stlo)
%SPHERICALINV    Return distance and azimuths between 2 locations on sphere
%
%    Usage:    [gcarc,az,baz]=sphericalinv(lat1,lon1,lat2,lon2)
%
%    Description:
%     [GCARC,AZ,BAZ]=SPHERICALINV(LAT1,LON1,LAT2,LON2) returns the
%     great-circle-arc degree distances GCARC, forward azimuths AZ and
%     backward azimuths BAZ between initial point(s) with geocentric
%     latitudes LAT1 and longitudes LON1 and final point(s) with geocentric
%     latitudes LAT2 and longitudes LON2 on a sphere.  All inputs must be
%     in degrees.  Outputs are also all in degrees.
%
%    Notes:
%     - Will always return the shorter great-circle-arc (GCARC<=180)
%     - GCARC accuracy degrades when < 3000km (see HAVERSINE example!)
%     - Azimuths are returned in the range 0<=az<=360
%
%    Examples:
%     % St. Louis, MO USA to Yaounde, Cameroon:
%     [dist,az,baz]=sphericalinv(38.649,-90.305,3.861,11.521)
%
%     % St. Louis, MO USA to Isla Isabella, Galapagos:
%     [dist,az,baz]=sphericalinv(38.649,-90.305,-0.823,-91.097)
%
%    See also: HAVERSINE, SPHERICALFWD, VINCENTYINV, VINCENTYFWD

%     Version History:
%        Oct. 14, 2008 - initial version
%        Oct. 26, 2008 - improved scalar expansion, doc and comment update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Apr. 10, 2010 - fix for colocated positions giving complex gcarc
%        Jan. 22, 2011 - use degrees functions, nargchk fix
%        Feb. 10, 2011 - force equal positions to give gcarc=0,az=0,baz=0
%        Feb.  9, 2012 - doc update, drop replication for speed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2012 at 13:15 GMT

% todo:

% require 4 inputs
error(nargchk(4,4,nargin));

% size up inputs
sz{1}=size(evla); sz{2}=size(evlo);
sz{3}=size(stla); sz{4}=size(stlo);
n(1)=prod(sz{1}); n(2)=prod(sz{2});
n(3)=prod(sz{3}); n(4)=prod(sz{4});

% basic check inputs
if(~isnumeric(evla) || ~isnumeric(evlo) ||...
        ~isnumeric(stla) || ~isnumeric(stlo))
    error('seizmo:sphericalinv:nonNumeric','All inputs must be numeric!');
elseif(sum(n~=1)>1 && ~isequal(sz{n~=1}))
    error('seizmo:sphericalinv:badSize',...
        'All inputs must be equal sized or scalar!');
end

% shortcut for just distance
if(nargout<2)
    % get law-of-cosines distance
    % - use real to avoid occasional complex values when points coincide
    gcarc=real(acosd(sind(evla).*sind(stla)...
        +cosd(evla).*cosd(stla).*cosd(stlo-evlo)));
    
    % force equal points to be 0,0,0
    eqpo=(evla==stla & evlo==stlo);
    if(any(eqpo(:))); gcarc(eqpo)=0; end
    return;
end

% optimize computation over memory
sinlat1=sind(evla);
sinlat2=sind(stla);
coslat1=cosd(evla);
coslat2=cosd(stla);
coslo=cosd(stlo-evlo);
sinlo=sind(stlo-evlo);

% get law-of-cosines distance
% - use real to avoid occasional complex values when points coincide
gcarc=real(acosd(sinlat1.*sinlat2+coslat1.*coslat2.*coslo));

% for conversion
r2d=180/pi;

% azimuths
az=mod(r2d.*atan2(sinlo.*coslat2,...
    coslat1.*sinlat2-sinlat1.*coslat2.*coslo),360);
baz=mod(r2d.*atan2(-sinlo.*coslat1,...
    coslat2.*sinlat1-sinlat2.*coslat1.*coslo),360);

% force equal points to be 0,0,0
eqpo=(evla==stla & evlo==stlo);
if(any(eqpo(:))); gcarc(eqpo)=0; az(eqpo)=0; baz(eqpo)=0; end

end
