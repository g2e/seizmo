function [kx,ky]=slowbaz2kxy(s,baz,f)
%SLOWBAZ2KXY    Converts slowness and back-azimuth to wavenumbers in x & y
%
%    Usage:    [kx,ky]=slowbaz2kxy(s,baz,f)
%
%    Description:
%     [Kx,Ky]=SLOWBAZ2KXY(S,BAZ,F) converts slowness S, backazimuth BAZ and
%     frequency F into wavenumbers in the East (Kx) and North (Ky)
%     directions.  Positive Kx correspond to waves from the East heading
%     West and positive Ky to waves from the North going South.  The units
%     of Kx & Ky depend on the units of S & F.  The distance units of Kx &
%     Ky will match that of S, so if S is in seconds per degree then Kx &
%     Ky will be in wavelengths per degree.  If F is an angular frequency
%     then Kx & Ky will be in radians per distance unit.  If F is not
%     given, it is assumed to be 1Hz.
%
%    Notes:
%     - S, BAZ & F must match in size or be scalar
%
%    Examples:
%     % Get wavenumbers for a wave coming from the SouthEast at 3km/s:
%     [kx,ky]=slowbaz2kxy(1/3,135)
%
%     % Wavenumbers in rad/km for a 10s wave moving 4.43s/deg from 312deg:
%     [kx,ky]=slowbaz2kxy(4.43*180/pi/6371,312,2*pi/10)
%
%    See also: KXY2SLOWBAZ, FKARF

%     Version History:
%        May   1, 2010 - initial version
%        June 16, 2010 - fixed nargchk
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 14:15 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% check inputs
if(nargin==2 || isempty(f)); f=1; end
if(~isreal(s) || ~isreal(baz) || ~isreal(f))
    error('seizmo:slowbaz2kxy:badInput',...
        'S, BAZ, & F must all be real valued!');
end
if(~isequalsizeorscalar(s,baz,f))
    error('seizmo:slowbaz2kxy:badInput',...
        'S, BAZ, & F must all be equal sized or scalar!');
end

% get wavenumbers
baz=baz*pi/180;
kx=f.*s.*sin(baz);
ky=f.*s.*cos(baz);

end
