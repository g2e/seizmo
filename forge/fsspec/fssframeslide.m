function [mov]=fssframeslide(s,varargin)
%FSSFRAMESLIDE    Slides through times of the freq-slow spectra
%
%    Usage:    mov=fssframeslide(s)
%              mov=fssframeslide(s,...)
%
%    Description:
%     MOV=FSSFRAMESLIDE(S) slides through the times of the frequency-
%     slowness spectra S, plotting each time individually and saving it as
%     a frame in the output movie MOV.  This will average over frequency
%     so you may want to pass S to FSSSUB to prefilter the spectra.
%
%     MOV=FSSFRAMESLIDE(S,...) passes additional options to PLOTFSS.
%     See that function for more details.
%
%    Notes:
%
%    Examples:
%     % Make a movie and save to an AVI file:
%     mov=fssframeslide(s);
%     movie2avi(mov,'example.avi');
%     unixcompressavi('example.avi');
%
%    See also: FSS, FSSXC, FSSHORZ, FSSHORZXC, FSSSUB, PLOTFSS,
%              PLOTFSSUPDATE, FSSFREQSLIDE, FSSAVG

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Apr.  4, 2012 - minor doc update
%        Apr. 25, 2012 - allow non-volume in slowness domain
%        June  8, 2012 - adapted from fkfreqslide
%        Sep. 13, 2012 - adapted from geofssframeslide
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 13, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check struct
error(chkfss(s));

% get frames
nf=numel(s);

% make initial plot
ax=plotfss(fssavg(s(1)),varargin{:});
fh=get(ax,'parent');
mov=getframe(fh);

% loop over remaining frames
for i=2:nf
    plotfssupdate(fssavg(s(i)),ax);
    mov(i)=getframe(fh);
end

end
