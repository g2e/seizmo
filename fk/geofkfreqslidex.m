function [varargout]=geofkfreqslidex(vol,frng,srng,delay,varargin)
%GEOFKFREQSLIDEX    Marks peak while sliding through geofk frequencies
%
%    Usage:    geofkfreqslidex(vol)
%              geofkfreqslidex(vol,frng)
%              geofkfreqslidex(vol,frng,srng)
%              geofkfreqslidex(vol,frng,srng,delay)
%              geofkfreqslidex(vol,frng,srng,delay,dblim)
%              geofkfreqslidex(vol,frng,srng,delay,dblim,zerodb)
%              geofkfreqslidex(vol,frng,srng,delay,dblim,zerodb,...
%                              'mmap_opt1',mmap_val1,...)
%              mov=geofkfreqslidex(...);
%
%    Description:
%     GEOFKFREQSLIDEX(VOL) slides through a geofk volume by plotting each
%     of the frequencies in sequence in a single plot.  Note that the
%     volume is averaged across all slownesses (see the following Usage
%     forms to adjust this).  There is a 1/3 second delay between each
%     replotting by default (see the following Usage forms to adjust this).
%     This differs from GEOFKFREQSLIDE in that an "X" marker is placed at
%     the maximum db point for each frame, giving a visual reference to
%     follow.
%
%     GEOFKFREQSLIDEX(VOL,FRNG) sets the frequency range to slide through.
%     FRNG gives the frequency range to extract as [FREQLOW FREQHIGH] in
%     Hz.  The default is [] and will slide across all frequencies in VOL.
%
%     GEOFKFREQSLIDEX(VOL,FRNG,SRNG) defines the range of horizontal
%     slownesses to __AVERAGE__ across in seconds per degree as
%     [SLOWLOW SLOWHIGH].  The default is [] and will average across all
%     slownesses in the VOL.
%
%     GEOFKFREQSLIDEX(VOL,FRNG,SRNG,DELAY) specifies the delay between the
%     plotting of each frequency in seconds.  The default DELAY is 0.33s.
%
%     GEOFKFREQSLIDEX(VOL,FRNG,SRNG,DELAY,DBLIM) sets the dB limits
%     for coloring the beam data.  The default is [-12 0] for the default
%     ZERODB (see next Usage form).  If ZERODB IS 'min' or 'median', the
%     default DBLIM is [0 12].  DBLIM must be a real-valued 2-element
%     vector.
%
%     GEOFKFREQSLIDEX(VOL,FRNG,SRNG,DELAY,DBLIM,ZERODB) changes what
%     0dB corresponds to in the plot.  The allowed values are 'min', 'max',
%     'median', & 'abs'.  The default is 'max'.
%
%     GEOFKFREQSLIDEX(VOL,FRNG,SRNG,DELAY,DBLIM,ZERODB,...
%                     'MMAP_OPT1',MMAP_VAL1,...) passes additional options
%     on to MMAP to alter the map.
%
%     MOV=GEOFKFREQSLIDEX(...) creates a Matlab movie MOV.  This can played
%     using the function MOVIE.  See MOVIE2AVI for exporting as an AVI
%     file.
%
%    Notes:
%
%    Examples:
%
%    See also: GEOFKFREQSLIDE, GEOFKXCVOLUME, PLOTGEOFKMAP, UPDATEGEOFKMAP,
%              GEOFKSLOWSLIDE, GEOFKVOL2MAP, GEOFKXCHORZVOLUME

%     Version History:
%        July 18, 2010 - initial version
%        Apr.  4, 2012 - minor doc update
%        Mar.  5, 2014 - update doc for plotgeofkmap input changes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  5, 2014 at 23:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check geofk struct
error(chkgeofkstruct(vol));

% don't allow array/volume
if(~isscalar(vol) || ~all(vol.volume))
    error('seizmo:geofkfreqslidex:badInput',...
        'VOL must be a scalar geofk struct and a beam volume!');
end

% do we make the movie
makemovie=false;
if(nargout); makemovie=true; end

% check delay
if(nargin<2); frng=[]; end
if(nargin<3); srng=[]; end
if(nargin<4 || isempty(delay)); delay=0.33; end
sf=size(frng);
ss=size(srng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<=0)))
    error('seizmo:geofkfreqslidex:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
end
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:geofkfreqslidex:badInput',...
        'SRNG must be a 1x2 array of [SLOWLOW SLOWHIGH] in sec/deg!');
end
if(~isreal(delay) || ~isscalar(delay) || delay<0)
    error('seizmo:geofkfreqslidex:badDelay',...
        'DELAY must be a positive scalar in seconds!');
end

% extract subvolume
if(~isempty(frng) || ~isempty(srng))
    vol=geofksubvol(vol,frng,srng);
end

% average across slowness
nslow=size(vol.beam,2);
vol.beam=sum(vol.beam,2)/nslow;
vol.volume(1)=false;

% get frequencies
freqs=vol.freq;
nfreq=numel(freqs);

% make initial plot
map=geofkvol2map(vol,[freqs(1) freqs(1)]);
maxdb=geofkdbinfo(map);
ax=plotgeofkmap(map,varargin{:});
hold on
x(1)=m_scatter(maxdb.latlon(1,2),maxdb.latlon(1,1),200,'kx','linewidth',5);
x(2)=m_scatter(maxdb.latlon(1,2),maxdb.latlon(1,1),100,'wx','linewidth',2);
hold off
fh=get(ax,'parent');
if(makemovie); varargout{1}=getframe(fh); else drawnow; end

% loop over remaining frequencies
for i=2:nfreq
    pause(delay);
    map=geofkvol2map(vol,[freqs(i) freqs(i)]);
    maxdb=geofkdbinfo(map);
    updategeofkmap(map,ax);
    delete(x);
    hold on
    x(1)=m_scatter(maxdb.latlon(1,2),maxdb.latlon(1,1),200,'kx',...
        'linewidth',5);
    x(2)=m_scatter(maxdb.latlon(1,2),maxdb.latlon(1,1),100,'wx',...
        'linewidth',2);
    hold off
    if(makemovie); varargout{1}(i)=getframe(fh); else drawnow; end
end

end
