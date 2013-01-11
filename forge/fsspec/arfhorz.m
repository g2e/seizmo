function [r,t]=arfhorz(varargin)
%ARFHORZ    Returns array response functions of horizontals
%
%    Usage:    [r,t]=arfhorz(stlalo,smax,spts,baz0,slow0,f0,radial0)
%              [r,t]=arfhorz(...,'polar',true|false,...)
%              [r,t]=arfhorz(...,'method',string,...)
%              [r,t]=arfhorz(...,'weights',w,...)
%              [r,t]=arfhorz(...,'avg',true|false,...)
%
%    Description:
%     [R,T]=ARFHORZ(STLALO,SMAX,SPTS,baz0,slow0,f0,radial0) computes the
%     radial and transverse array response functions (ARFs) for an array at
%     the locations in STLALO with horizontal-motion plane waves passing
%     through the array defined by backazimuths baz0 (in deg), horizontal
%     slownesses slow0 (in sec/deg), frequencies f0 (in Hz), and motion
%     radial0 (true=radial, false=trans).  The ARFs are computed in
%     cartesian slowness space where the range of the sampling is given by
%     SMAX (sec/deg) and extends from -SMAX to SMAX for both East/West and
%     North/South directions.  SPTS controls the number of slowness points
%     for both directions (SPTSxSPTS grid).  Latitudes & longitudes are
%     expected in degrees with the position arrays organized as [LAT LON].
%     The ARF is created using the 'center' method as this is faster (see
%     below to adjust the method).  The output struct S has the same fields
%     as FSSHORZ output but also includes a S.source field which has the
%     following layout:
%      .source.nsrc   - number of plane wave sources
%      .source.baz    - backazimuth (deg)
%      .source.slow   - horizontal slowness (sec/deg)
%      .source.freq   - frequency (Hz)
%      .source.radial - is radial motion (true) or transverse (false)
%     By default SPTS is 101, baz0 is 0, slow0 is 0, f0 is 1, & radial0 is
%     true.
%
%     [R,T]=ARFHORZ(...,'POLAR',TRUE|FALSE,...) decides if the spectra of
%     the array response is sampled regularly with cartesian or polar
%     coordinates.  Polar coords are useful for slicing the spectra by
%     azimuth (pie slice) or slowness (rings).  Cartesian coords (the
%     default) samples the slowness space regularly in the East/West &
%     North/South directions and so exhibits less distortion in plots of
%     the slowness space. If POLAR=TRUE, SPTS may be [spnts bazpnts] to
%     control the azimuthal resolution (default is 181 points).
%
%     [R,T]=ARFHORZ(...,'METHOD',STRING,...) sets the beamforming method.
%     STRING may be 'center', 'coarray', 'full', or [LAT LON].
%
%     [R,T]=ARFHORZ(...,'WEIGHTS',W,...) specifies the relative weights for
%     each station (must have the same number of rows as STLALO) or pair
%     (only if METHOD is 'coarray' or 'full').
%
%     [R,T]=ARFHORZ(...,'AVG',TRUE|FALSE,...) indicates if the spectra is
%     averaged across frequency during computation.  This can save a
%     significant amount of memory.  The default is false.
%
%    Notes:
%
%    Examples:
%     % Compare the vertical & horizontal 1Hz array response functions
%     % out to 5 sec/deg for an array in a SEIZMO dataset:
%     plotarf(arf(getheader(data,'st'),5));
%     [r,t]=arfhorz(getheader(data,'st'),5);
%     plotarf(r);
%     plotarf(t);
%
%     % Get multi-plane wave responses at 0.03Hz for a random array:
%     st=randlatlon(20)/45;
%     [r,t]=arfhorz(st,50,201,[0 0 45],[20 10 20],0.03,[true false true]);
%     plotarf(fssavg(r));
%     plotarf(fssavg(t));
%
%    See also: PLOTARF, ARF, FSSHORZ, FSSHORZXC, FSS, FSSXC, SNYQUIST

%     Version History:
%        Sep. 22, 2012 - initial version
%        Sep. 27, 2012 - pv pair inputs, doc update
%        Jan.  8, 2013 - avg option (tricky)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan.  8, 2013 at 14:40 GMT

% todo:

% defaults
if(nargin<4 || isempty(varargin{4})); baz0=0; else baz0=varargin{4}; end
if(nargin<5 || isempty(varargin{5})); slow0=0; else slow0=varargin{5}; end
if(nargin<6 || isempty(varargin{6})); f0=1; else f0=varargin{6}; end
if(nargin<7 || isempty(varargin{7})); radial0=true;
else radial0=varargin{7};
end

% check radial parameter
if(~islogical(radial0))
    error('seizmo:arfhorz:badInput',...
        'radial0 must be true (radial motion) or false (transverse)!');
elseif(~isequalsizeorscalar(radial0,slow0,baz0,f0))
    error('seizmo:arfhorz:badInput',...
        'baz0, slow0, f0, radial0 must be equal sized or scalar!');
end

% let arf compute the undistorted array response
r=arf(varargin{1:min(6,nargin)},varargin{8:end},'avg',false);
t=r;
r.orient='radial';
t.orient='transverse';

% expand needed planewave parameters
npw=size(r.spectra,3);
if(isscalar(baz0)); baz0=baz0(ones(npw,1)); end
if(isscalar(radial0)); radial0=radial0(ones(npw,1)); end

% build slowness baz grid
if(r.polar)
    baz=r.x(ones(numel(r.y),1),:);
else % cartesian
    baz=atan2(r.x(ones(numel(r.y),1),:),r.y(:,ones(numel(r.x),1)))*180/pi;
end

% loop over planewaves
for i=1:npw
    if(radial0(i))
        r.spectra(:,:,i)=r.spectra(:,:,i).*abs(cosd(baz-baz0(i))).^2;
        t.spectra(:,:,i)=t.spectra(:,:,i).*abs(sind(baz-baz0(i))).^2;
    else
        r.spectra(:,:,i)=r.spectra(:,:,i).*abs(sind(baz-baz0(i))).^2;
        t.spectra(:,:,i)=t.spectra(:,:,i).*abs(cosd(baz-baz0(i))).^2;
    end
end

% average
pv.avg=false;
for i=9:2:nargin
    switch lower(varargin{i})
        case {'average' 'av' 'avg' 'a'}
            pv.avg=varargin{i+1};
    end
end
if(pv.avg)
    r.spectra=mean(r.spectra,3);
    r.spectra=mean(t.spectra,3);
end

end
