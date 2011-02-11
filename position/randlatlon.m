function [varargout]=randlatlon(n)
%RANDLATLON    Returns lat/lon points uniformly distributed on a sphere
%
%    Usage:    latlon=randlatlon(n)
%              [lat,lon]=randlatlon(n)
%
%    Description:
%     LATLON=RANDLATLON(N) returns the random position of N points
%     uniformly distributed on a sphere.  The output LATLON is formatted as
%     a Nx2 matrix of [LAT LON] where the latitudes and longitudes are in
%     degrees.  Longitudes are in the range -180 to 180 degrees.
%
%     [LAT,LON]=RANDLATLON(N) returns the latitude and longitude as
%     separate Nx1 column vectors.
%
%    Notes:
%     - uses RANDSPHERE from Matlab File Exchange (included as subfunction)
%
%    Examples:
%     % Get a single point in latitude & longitude:
%     randlatlon
%
%     % Now return 1000 points & plot them on a map:
%     mmap('events',randlatlon(1000),'proj','hammer')
%
%    See also: RANDSPHERE, XYZ2GEOCENTRIC, RAND

%     Version History:
%        Jan.  4, 2011 - initial version
%        Feb. 10, 2011 - fix example for mmap
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2011 at 02:15 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% default nargin
if(nargin<1); n=1; end

% check n
if(~isscalar(n) || ~isreal(n) || n~=fix(n))
    error('seizmo:randlatlon:badInput',...
        'N must be a scalar integer value!');
end

% random lat/lon
xyz=randsphere(n,3,1);
[lat,lon]=xyz2geocentric(xyz(:,1),xyz(:,2),xyz(:,3));

% output style
if(nargout<2)
    varargout{1}=[lat lon];
else
    varargout={lat lon};
end

end

function x=randsphere(m,n,r)
% Roger Stafford - 12/23/05
x=randn(m,n);
s2=sum(x.^2,2);
x=x.*repmat(r*(gammainc(s2/2,n/2).^(1/n))./sqrt(s2),1,n);
end
