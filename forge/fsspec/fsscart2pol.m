function [s]=fsscart2pol(s,azipnts,slowpnts,method)
%FSSCART2POL    Converts a cartesian grid fss spectra to a polar grid
%
%    Usage:    s=fsscart2pol(s)
%              s=fsscart2pol(s,azipnts,slowpnts)
%              s=fsscart2pol(s,azipnts,slowpnts,method)
%
%    Description:
%     S=FSSCART2POL(S) converts an frequency-slowness spectra regularly
%     sampled in cartesian space (ie. Seast & Snorth) to one regularly
%     sampled in polar space (ie. backazimuth & |S|).  Please note that
%     this uses 2D-based interpolation!
%
%     S=FSSCART2POL(S,AZIPNTS,SLOWPNTS) specifies the number of azimuthal
%     and horizontal slowness samples of the polar grid.  The number of
%     azimuthal points AZIPNTS is by default 361 while the number of
%     slowness points defaults to give about the same resolution as that of
%     the cartesian spectra.
%
%     S=FSSCART2POL(S,AZIPNTS,SLOWPNTS,METHOD) specifies the 2D
%     interpolation method.  See INTERP2D for possible options.
%
%    Notes:
%     - Expands to contain the whole grid (no data loss).
%     - Polar grids are passed through without action.
%
%    Examples:
%     % Plot a cartesian spectra, convert to polar,
%     % and plot again to compare:
%     plotfss(s)
%     plotfss(fsscart2pol(s))
%
%    See also: FSSPOL2CART, FSS, FSSXC, FSSHORZ, FSSHORZXC, ARF

%     Version History:
%        July 14, 2010 - initial version
%        Apr.  4, 2012 - minor doc update
%        Sep. 14, 2012 - adapted from fkcart2pol, use nsteps, no data loss
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
if(nargin<2 || isempty(azipnts)); azipnts=361; end
if(nargin<3 || isempty(slowpnts))
    slowpnts=nan(nfss,1);
    for i=1:nfss
        slowpnts(i)=fix(sqrt(2)*max(s(i).x)/abs(min(diff(s(i).x)))+1);
    end
end
if(nargin<4 || isempty(method)); method='linear'; end

% expand step sizes
if(isscalar(azipnts)); azipnts=azipnts(ones(nfss,1),1); end
if(isscalar(slowpnts)); slowpnts=slowpnts(ones(nfss,1),1); end

% check optional inputs
if(~isreal(azipnts) || ~isequal(numel(s),numel(azipnts)) ...
        || any(azipnts<=0))
    error('seizmo:fsscart2pol:badInput',...
        'AZIPNTS must be a scalar or equal in size to S & positive!');
end
if(~isreal(slowpnts) || ~isequal(numel(s),numel(slowpnts)) ...
        || any(slowpnts<=0))
    error('seizmo:fsscart2pol:badInput',...
        'SLOWPNTS must be a scalar or equal in size to S & positive!');
end
if(~ischar(method) || size(method,1)~=1 || ndims(method)>2)
    error('seizmo:fsscart2pol:badInput',...
        'METHOD must be a string!');
end

% loop over each fss element
for i=1:nfss
    % skip if polar
    if(s(i).polar); continue; end
    
    % previous x/y grid
    [east,north]=meshgrid(s(i).x,s(i).y);
    
    % azi/slow grid
    smax=sqrt(2)*max(s(i).x);
    s(i).x=(0:azipnts(i)-1)/(azipnts(i)-1)*360;
    s(i).y=(0:slowpnts(i)-1).'/(slowpnts(i)-1)*smax;
    [baz,slow]=meshgrid(s(i).x,s(i).y);
    
    % convert to x/y (east/north)
    [easti,northi]=slowbaz2kxy(slow,baz);
    
    % interpolate
    nfreq=size(s(i).spectra,3);
    zi=nan(slowpnts(i),azipnts(i),nfreq);
    for j=1:nfreq
        zi(:,:,j)=interp2(east,north,s(i).spectra(:,:,j),...
            easti,northi,method);
    end
    
    % assign
    s(i).spectra=zi;
    s(i).polar=true;
    
    % remove slowness magnitudes with all nans
    % - this is for repeated calls between fsspol2cart & fsscart2pol
    nans=isnan(s(i).spectra);
    lastgoodslow=find(mean(nans(:,:),2)<1,1,'last');
    s(i).y(lastgoodslow+1:end)=[];
    s(i).spectra=s(i).spectra(1:lastgoodslow,:,:);
end

end
