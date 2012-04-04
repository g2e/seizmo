function [varargout]=geofkarfslowslide(arf,srng,delay,varargin)
%GEOFKARFSLOWSLIDE    Slides through the slownesses of a geofkarf volume
%
%    Usage:    geofkarfslowslide(arf)
%              geofkarfslowslide(arf,srng)
%              geofkarfslowslide(arf,srng,delay)
%              geofkarfslowslide(arf,srng,delay,projopt)
%              geofkarfslowslide(arf,srng,delay,projopt,dblim)
%              geofkarfslowslide(arf,srng,delay,projopt,dblim,zerodb)
%              geofkarfslowslide(arf,srng,delay,projopt,dblim,zerodb,...
%                             fgcolor,bgcolor)
%              geofkarfslowslide(arf,srng,delay,projopt,dblim,zerodb,...
%                             fgcolor,bgcolor,ax)
%              mov=geofkarfslowslide(...);
%
%    Description:
%     GEOFKARFSLOWSLIDE(ARF) slides through a geofk ARF by plotting each of
%     the slownesses in sequence in a single plot.  There is a 1/3 second
%     delay between each replotting by default (see the following Usage
%     forms to adjust this).
%
%     GEOFKARFSLOWSLIDE(ARF,SRNG) defines the range of horizontal
%     slownesses to slide across in seconds per degree as
%     [SLOWLOW SLOWHIGH].  The default is [] and will slide across all
%     slownesses in the ARF.
%
%     GEOFKARFSLOWSLIDE(ARF,SRNG,DELAY) specifies the delay between the
%     plotting of each slowness in seconds.  The default DELAY is 0.33s.
%
%     GEOFKARFSLOWSLIDE(ARF,SRNG,DELAY,PROJOPT) allows passing options to
%     M_PROJ.  See M_PROJ('SET') for possible projections and see
%     M_PROJ('GET',PROJ) for a list of possible additional options specific
%     to that projection.
%
%     GEOFKARFSLOWSLIDE(ARF,SRNG,DELAY,PROJOPT,DBLIM) sets the dB limits
%     for coloring the beam data.  The default is [-12 0] for the default
%     ZERODB (see next Usage form).  If ZERODB IS 'min' or 'median', the
%     default DBLIM is [0 12].  DBLIM must be a real-valued 2-element
%     vector.
%
%     GEOFKARFSLOWSLIDE(ARF,SRNG,DELAY,PROJOPT,DBLIM,ZERODB) changes what
%     0dB corresponds to in the plot.  The allowed values are 'min', 'max',
%     'median', & 'abs'.  The default is 'max'.
%
%     GEOFKARFSLOWSLIDE(ARF,SRNG,DELAY,PROJOPT,DBLIM,ZERODB,FGCOLOR,...
%     BGCOLOR) specifies the foreground and background colors of the plot. 
%     The default is 'w' for FGCOLOR and 'k' for BGCOLOR.  Note that if one
%     is specified and the other is not, an opposing color is found using
%     INVERTCOLOR.  The color scale is also changed so the noise clip is at
%     BGCOLOR.
%
%     GEOFKARFSLOWSLIDE(ARF,SRNG,DELAY,PROJOPT,DBLIM,ZERODB,FGCOLOR,...
%     BGCOLOR,AX) sets the axes that the map is drawn in.  This is useful
%     for subplots, guis, etc.
%
%     MOV=GEOFKARFSLOWSLIDE(...) creates a Matlab movie MOV.  This can
%     played using the function MOVIE.  See MOVIE2AVI for exporting as an
%     AVI file.
%
%    Notes:
%
%    Examples:
%
%    See also: GEOFKARF, PLOTGEOFKARF, GEOFKARF2MAP, GEOFKSUBARF,
%              UPDATEGEOFKARF, CHKGEOFKARFSTRUCT

%     Version History:
%        July  8, 2010 - initial version
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 10:05 GMT

% todo:

% check nargin
error(nargchk(1,10,nargin));

% check geofk struct
error(chkgeofkarfstruct(arf));

% don't allow array/volume
if(~isscalar(arf) || any(~sum(reshape([arf.volume],[2 numel(arf)]))))
    error('seizmo:geofkarfslowslide:badInput',...
        'ARF must be a scalar geofkarf struct and a volume!');
end

% do we make the movie
makemovie=false;
if(nargout); makemovie=true; end

% check delay
if(nargin<2); srng=[]; end
if(nargin<3 || isempty(delay)); delay=0.33; end
ss=size(srng);
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:geofkarfslowslide:badInput',...
        'SRNG must be a 1x2 array of [SLOWLOW SLOWHIGH] in sec/deg!');
end
if(~isreal(delay) || ~isscalar(delay) || delay<0)
    error('seizmo:geofkarfslowslide:badDelay',...
        'DELAY must be a positive scalar in seconds!');
end

% extract subvolume
arf=geofksubarf(arf,srng);

% get slownesses
slows=arf.horzslow;
nslow=numel(slows);

% make initial plot
ax=plotgeofkarf(geofkarf2map(arf,[slows(1) slows(1)]),varargin{:});
fh=get(ax,'parent');
if(makemovie); varargout{1}=getframe(fh); end

% loop over remaining frequencies
for i=2:nslow
    pause(delay);
    updategeofkarf(geofkarf2map(arf,[slows(i) slows(i)]),ax);
    if(makemovie); varargout{1}(i)=getframe(fh); end
end

end
