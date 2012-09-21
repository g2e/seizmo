function [s]=fsspol2cart(s,slowpnts,method)
%FSSPOL2CART    Converts a polar grid fss spectra to a cartesian grid
%
%    Usage:    s=fsspol2cart(s)
%              s=fsspol2cart(s,slowpnts)
%              s=fsspol2cart(s,slowpnts,method)
%
%    Description:
%     S=FSSPOL2CART(S) converts an frequency-slowness spectra regularly
%     sampled in polar space (ie. backazimuth & |S|) to one regularly
%     sampled in cartesian space (ie. Seast & Snorth).  Please note that
%     this uses 2D-based interpolation!
%
%     S=FSSPOL2CART(S,SLOWPNTS) specifies the number of horizontal slowness
%     points of the output cartesian grid (SLOWPNTSxSLOWPNTS).  The number
%     of slowness points defaults to give about the same resolution as that
%     of the input polar spectra.
%
%     S=FSSPOL2CART(S,SLOWPNTS,METHOD) specifies the 2D interpolation
%     method.  See INTERP2D for possible options.
%
%    Notes:
%     - Expands to contain the whole input grid (no data loss).
%     - Cartesian grids are passed through without action.
%
%    Examples:
%     % Plot a polar spectra, convert to cartesian,
%     % and plot again to compare:
%     plotfss(s)
%     plotfss(fsspol2cart(s))
%
%    See also: FSSCART2POL, FSS, FSSXC, FSSHORZ, FSSHORZXC, ARF

%     Version History:
%        Sep. 14, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 14, 2012 at 16:05 GMT

% todo:

% check nargin
error(nargchk(1,4,nargin));

% check fss struct
error(chkfss(s));
nfss=numel(s);

% defaults
if(nargin<2 || isempty(slowpnts))
    slowpnts=nan(nfss,1);
    for i=1:nfss
        slowpnts(i)=fix(2*max(s(i).y)/abs(min(diff(s(i).y)))+1);
    end
end
if(nargin<3 || isempty(method)); method='linear'; end

% expand step sizes
if(isscalar(slowpnts)); slowpnts=slowpnts(ones(nfss,1),1); end

% check optional inputs
if(~isreal(slowpnts) || ~isequal(numel(s),numel(slowpnts)) ...
        || any(slowpnts<=0))
    error('seizmo:fsspol2cart:badInput',...
        'SLOWPNTS must be a scalar or equal in size to S & positive!');
end
if(~ischar(method) || size(method,1)~=1 || ndims(method)>2)
    error('seizmo:fsspol2cart:badInput',...
        'METHOD must be a string!');
end

% loop over each fss element
for i=1:nfss
    % skip if cartesian
    if(~s(i).polar); continue; end
    
    % previous x/y grid
    [baz,slow]=meshgrid(s(i).x,s(i).y);
    
    % azi/slow grid
    smax=max(s(i).y);
    s(i).x=-smax:2*smax/(slowpnts(i)-1):smax;
    s(i).y=fliplr(s(i).x).';
    [east,north]=meshgrid(s(i).x,s(i).y);
    
    % convert to x/y (baz/slow)
    [slowi,bazi]=kxy2slowbaz(east,north);
    
    % convert bazi, -180 to 180 => 0 to 360
    bazi(bazi<0)=bazi(bazi<0)+360;
    
    % interpolate
    nfreq=size(s(i).spectra,3);
    zi=nan(slowpnts(i),slowpnts(i),nfreq);
    for j=1:nfreq
        zi(:,:,j)=interp2(baz,slow,s(i).spectra(:,:,j),...
            bazi,slowi,method);
    end
    
    % assign
    s(i).spectra=zi;
    s(i).polar=false;
end

end
