function [c]=lonmod(a,b)
%LONMOD    Returns a longitude modulus (ie unwraps longitudes)
%
%    Usage:    lon=lonmod(lon)
%              c=lonmod(a,b)
%
%    Description:
%     LON=LONMOD(LON) returns the longitudes LON within the range +/-180.
%     For example, LON=181 is output as -179.  See the other usage forms
%     for more algorithm details.
%
%     LONMOD(A,B) is A-N.*B where N=round(A./B) if B~=0.  Thus LONMOD(A,B)
%     is always within the range +/-(B/2).  The inputs A and B must be real
%     arrays of the same size, or real scalars.
%
%    Notes:
%     - By convention:
%        lonmod(A,0) returns A
%        lonmod(A,A) returns 0
%
%    Examples:
%     % Use FIXLATLON to unwrap both the latitude and longitude while
%     % preserving the actual position:
%     [LAT,LON]=fixlatlon(LAT,LON);
%
%    See also: LATMOD, FIXLATLON, MOD, REM

%     Version History:
%        Mar. 24, 2009 - initial version
%        Apr. 23, 2009 - move usage up
%        Sep.  7, 2009 - minor doc update
%        Sep. 30, 2009 - changed name from CMOD to LONMOD
%        Feb. 16, 2010 - doc update, preserve sign of -180/180
%        Feb.  9, 2012 - B=360 by default, doc update
%        Feb. 29, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 29, 2012 at 13:45 GMT

% todo:

% default range for longitudes
if(nargin<2 || isempty(b)); b=360; end

% this always flips the sign of -180/180
%c=a-round(a./b).*b;

% this preserves the sign of -180/180
c=a-sign(a).*ceil((abs(a)-b./2)./b).*b;

% handle zero case
d=b==0;
if(isscalar(b))
    if(d); c=a; end
else
    c(d)=a(d);
end

end
