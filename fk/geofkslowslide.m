function [varargout]=geofkslowslide(vol,frng,srng,delay,varargin)
%GEOFKSLOWSLIDE    Slides through a geofk volume plotting each slowness
%
%    Usage:    geofkslowslide(vol)
%              geofkslowslide(vol,frng)
%              geofkslowslide(vol,frng,srng)
%              geofkslowslide(vol,frng,srng,delay)
%              geofkslowslide(vol,frng,srng,delay,projopt)
%              geofkslowslide(vol,frng,srng,delay,projopt,dblim)
%              geofkslowslide(vol,frng,srng,delay,projopt,dblim,zerodb)
%              geofkslowslide(vol,frng,srng,delay,projopt,dblim,zerodb,...
%                             fgcolor,bgcolor)
%              geofkslowslide(vol,frng,srng,delay,projopt,dblim,zerodb,...
%                             fgcolor,bgcolor,ax)
%              mov=geofkslowslide(...);
%
%    Description:
%     GEOFKSLOWSLIDE(VOL) slides through a geofk volume by plotting each of
%     the slownesses in sequence in a single plot.  Note that the volume is
%     averaged across all frequencies (see the following Usage forms to
%     adjust this).  There is a 1/3 second delay between each replotting by
%     default (see the following Usage forms to adjust this).
%
%     GEOFKSLOWSLIDE(VOL,FRNG) sets the frequency range to _AVERAGE_.
%     FRNG gives the frequency range as [FREQLOW FREQHIGH] in Hz.  The
%     default is [] and will average across all frequencies in VOL.
%
%     GEOFKSLOWSLIDE(VOL,FRNG,SRNG) defines the range of horizontal
%     slownesses to slide across in seconds per degree as
%     [SLOWLOW SLOWHIGH].  The default is [] and will slide across all
%     slownesses in the VOL.
%
%     GEOFKSLOWSLIDE(VOL,FRNG,SRNG,DELAY) specifies the delay between the
%     plotting of each slowness in seconds.  The default DELAY is 0.33s.
%
%     GEOFKSLOWSLIDE(VOL,FRNG,SRNG,DELAY,PROJOPT) allows passing options to
%     M_PROJ.  See M_PROJ('SET') for possible projections and see
%     M_PROJ('GET',PROJ) for a list of possible additional options specific
%     to that projection.
%
%     GEOFKSLOWSLIDE(VOL,FRNG,SRNG,DELAY,PROJOPT,DBLIM) sets the dB limits
%     for coloring the beam data.  The default is [-12 0] for the default
%     ZERODB (see next Usage form).  If ZERODB IS 'min' or 'median', the
%     default DBLIM is [0 12].  DBLIM must be a real-valued 2-element
%     vector.
%
%     GEOFKSLOWSLIDE(VOL,FRNG,SRNG,DELAY,PROJOPT,DBLIM,ZERODB) changes what
%     0dB corresponds to in the plot.  The allowed values are 'min', 'max',
%     'median', & 'abs'.  The default is 'max'.
%
%     GEOFKSLOWSLIDE(VOL,FRNG,SRNG,DELAY,PROJOPT,DBLIM,ZERODB,FGCOLOR,...
%     BGCOLOR) specifies the foreground and background colors of the plot. 
%     The default is 'w' for FGCOLOR and 'k' for BGCOLOR.  Note that if one
%     is specified and the other is not, an opposing color is found using
%     INVERTCOLOR.  The color scale is also changed so the noise clip is at
%     BGCOLOR.
%
%     GEOFKSLOWSLIDE(VOL,FRNG,SRNG,DELAY,PROJOPT,DBLIM,ZERODB,FGCOLOR,...
%     BGCOLOR,AX) sets the axes that the map is drawn in.  This is useful
%     for subplots, guis, etc.
%
%     MOV=GEOFKSLOWSLIDE(...) creates a Matlab movie MOV.  This can played
%     using the function MOVIE.  See MOVIE2AVI for exporting as an AVI
%     file.
%
%    Notes:
%
%    Examples:
%
%    See also: GEOFKXCVOLUME, PLOTGEOFKMAP, UPDATEGEOFKMAP, GEOFKFREQSLIDE,
%              GEOFKVOL2MAP, GEOFKXCHORZVOLUME

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        July  8, 2010 - doc update
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 10:25 GMT

% todo:

% check nargin
error(nargchk(1,10,nargin));

% check geofk struct
error(chkgeofkstruct(vol));

% don't allow array/volume
if(~isscalar(vol) || ~all(vol.volume))
    error('seizmo:geofkslowslide:badInput',...
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
    error('seizmo:geofkslowslide:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
end
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:geofkslowslide:badInput',...
        'SRNG must be a 1x2 array of [SLOWLOW SLOWHIGH] in sec/deg!');
end
if(~isreal(delay) || ~isscalar(delay) || delay<0)
    error('seizmo:geofkslowslide:badDelay',...
        'DELAY must be a positive scalar in seconds!');
end

% extract subvolume
vol=geofksubvol(vol,frng,srng);

% average across frequency
nfreq=size(vol.beam,3);
vol.beam=sum(vol.beam,3)/nfreq;
vol.volume(2)=false;

% get slownesses
slows=vol.horzslow;
nslow=numel(slows);

% make initial plot
ax=plotgeofkmap(geofkvol2map(vol,[],[slows(1) slows(1)]),varargin{:});
fh=get(ax,'parent');
if(makemovie); varargout{1}=getframe(fh); end

% loop over remaining frequencies
for i=2:nslow
    pause(delay);
    updategeofkmap(geofkvol2map(vol,[],[slows(i) slows(i)]),ax);
    if(makemovie); varargout{1}(i)=getframe(fh); end
end

end
