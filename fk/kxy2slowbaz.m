function [s,baz]=kxy2slowbaz(kx,ky,f)
%KXY2SLOWBAZ    Converts wavenumbers in x & y to slowness and back-azimuth
%
%    Usage:    [s,baz]=kxy2slowbaz(kx,ky,f)
%
%    Description:
%     [S,BAZ]=KXY2SLOWBAZ(KX,KY,F) converts wavenumbers in the East/West
%     direction (Kx) and North/South direction (Ky) to slowness S and
%     backazimuth BAZ for frequency F.  Positive Kx correspond to waves
%     from the East heading West and positive Ky to waves from the North
%     going South.  The units of S depend on the units of Kx, Ky & F.  The
%     distance units of Kx & Ky will match that of S, so if Kx & Ky are
%     in wavelengths per degree, then S is in seconds per degree.  If F is
%     an angular frequency then Kx & Ky must be angular wavenumbers so S is
%     in seconds per distance unit.  If F is not given, it is assumed to be
%     1Hz.
%
%    Notes:
%     - Kx, Ky & F must match in size or be scalar
%
%    Examples:
%     % Get slowness and backazimuth for a wave with a period of 10s and a
%     % wavelength of 100km North/South and 20km East/West:
%     [s,baz]=kxy2slowbaz(1/20,1/100,1/10)
%
%    See also: SLOWBAZ2KXY, FKARF

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
if(~isreal(kx) || ~isreal(ky) || ~isreal(f))
    error('seizmo:kxy2slowbaz:badInput',...
        'Kx, Ky, & F must all be real valued!');
end
if(~isequalsizeorscalar(kx,ky,f))
    error('seizmo:kxy2slowbaz:badInput',...
        'Kx, Ky, & F must all be equal sized or scalar!');
end

% get slowness and backazimuth
s=sqrt(kx.^2+ky.^2)./f;
baz=atan2(kx,ky)*180/pi;

end
