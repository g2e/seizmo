function [mov]=geofssframeslide(s,varargin)
%GEOFSSFRAMESLIDE    Slides through times of the freq-slow spectra
%
%    Usage:    mov=geofssframeslide(s)
%              mov=geofssframeslide(s,...)
%
%    Description:
%     MOV=GEOFSSFRAMESLIDE(S) slides through the times of the
%     frequency-slowness-position spectra S, plotting each time
%     individually and saving it as a frame in the output movie MOV.  This
%     will average over frequency and slowness so you may want to pass S to
%     GEOFSSSUB to prefilter the spectra.
%
%     MOV=GEOFSSFRAMESLIDE(S,...) passes additional options to PLOTGEOFSS.
%     See that function for more details.
%
%    Notes:
%
%    Examples:
%     % Make a movie and save to an AVI file:
%     mov=geofssframeslide(s);
%     movie2avi(mov,'example.avi');
%     unixcompressavi('example.avi');
%
%    See also: GEOFSS, GEOFSSXC, GEOFSSSUB, PLOTGEOFSS, PLOTGEOFSSUPDATE,
%              GEOFSSFREQSLIDE, GEOFSSSLOWSLIDE, GEOFSSAVG

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Apr.  4, 2012 - minor doc update
%        Apr. 25, 2012 - allow non-volume in slowness domain
%        June  8, 2012 - adapted from geofkfreqslide
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  8, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check struct
error(chkgeofss(s));

% get frames
nf=numel(s);

% make initial plot
ax=plotgeofss(geofssavg(s(1)),varargin{:});
fh=get(ax,'parent');
mov=getframe(fh);

% loop over remaining frames
for i=2:nf
    plotgeofssupdate(geofssavg(s(i)),ax);
    mov(i)=getframe(fh);
end

end
