function [mov]=geofssslowslide(s,varargin)
%GEOFSSSLOWSLIDE    Slides through slownesses of the freq-slow spectra
%
%    Usage:    mov=geofssslowslide(s)
%              mov=geofssslowslide(s,...)
%
%    Description:
%     MOV=GEOFSSSLOWSLIDE(S) slides through the slownesses of the
%     frequency-slowness-position spectra S, plotting each slowness
%     individually and saving it as a frame in the output movie MOV.  This
%     will average over frequency so you may want to pass S to GEOFSSSUB to
%     prefilter the spectra.
%
%     MOV=GEOFSSSLOWSLIDE(S,...) passes additional options to PLOTGEOFSS.
%     See that function for more details.
%
%    Notes:
%
%    Examples:
%     % Make a movie and save to an AVI file:
%     mov=geofssslowslide(s);
%     movie2avi(mov,'example.avi');
%     unixcompressavi('example.avi');
%
%    See also: GEOFSS, GEOFSSXC, GEOFSSSUB, PLOTGEOFSS, PLOTGEOFSSUPDATE,
%              GEOFSSFREQSLIDE, GEOFSSAVG, GEOFSSFRAMESLIDE

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Apr.  4, 2012 - minor doc update
%        Apr. 25, 2012 - allow non-volume in slowness domain
%        June  8, 2012 - adapted from geofkfreqslide
%        Sep. 29, 2012 - update for struct & geofsssub changes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 29, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check struct
error(chkgeofss(s));

% require scalar geofss struct
if(~isscalar(s))
    error('seizmo:geofssslowslide:badInput',...
        'S must be a scalar GEOFSS struct!');
end

% require slowness vector
if(size(s.spectra,2)~=numel(s.slow))
    error('seizmo:geofssslowslide:badInput',...
        'S must NOT be averaged with GEOFSSAVG!');
end

% make initial plot
ax=plotgeofss(geofssavg(geofsssub(s,[],[],[],[],[],1)),varargin{:});
fh=get(ax,'parent');
mov=getframe(fh);

% loop over remaining frequencies
for i=2:numel(s.slow)
    plotgeofssupdate(geofssavg(geofsssub(s,[],[],[],[],[],i)),ax);
    mov(i)=getframe(fh);
end

end
