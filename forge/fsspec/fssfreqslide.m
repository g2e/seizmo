function [mov]=fssfreqslide(s,varargin)
%FSSFREQSLIDE    Slides through frequencies of the freq-slow spectra
%
%    Usage:    mov=fssfreqslide(s)
%              mov=fssfreqslide(s,...)
%
%    Description:
%     MOV=FSSFREQSLIDE(S) slides through the frequencies of the frequency-
%     slowness spectra S, plotting each frequency individually and saving
%     it as a frame in the output movie MOV.
%
%     MOV=FSSFREQSLIDE(S,...) passes additional options to PLOTFSS.
%     See that function for more details.
%
%    Notes:
%
%    Examples:
%     % Make a movie and save to an AVI file:
%     mov=fssfreqslide(s);
%     movie2avi(mov,'example.avi');
%     unixcompressavi('example.avi');
%
%    See also: FSS, FSSXC, FSSHORZ, FSSHORZXC, FSSSUB, PLOTFSS,
%              PLOTFSSUPDATE, FSSFRAMESLIDE, FSSAVG

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Apr.  4, 2012 - minor doc update
%        Apr. 25, 2012 - allow non-volume in slowness domain
%        June  8, 2012 - adapted from geofkfreqslide
%        Sep. 13, 2012 - adapted from geofssfreqslide
%        Sep. 29, 2012 - allow scalar freq (non-avg) spectra
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 29, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check struct
error(chkfss(s));

% require scalar geofss struct
if(~isscalar(s))
    error('seizmo:fssfreqslide:badInput',...
        'S must be a scalar FSS struct!');
end

% require frequency vector
if(size(s.spectra,3)~=numel(s.freq))
    error('seizmo:fssfreqslide:badInput',...
        'S must not be averaged with FSSAVG!');
end

% make initial plot
ax=plotfss(fsssub(s,[],1),varargin{:});
fh=get(ax,'parent');
mov=getframe(fh);

% loop over remaining frequencies
for i=2:numel(s.freq)
    plotfssupdate(fsssub(s,[],i),ax);
    mov(i)=getframe(fh);
end

end
