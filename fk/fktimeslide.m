function [varargout]=fktimeslide(map,delay,varargin)
%FKTIMESLIDE    Slides through a sequence of fk maps plotting each one
%
%    Usage:    fktimeslide(map)
%              fktimeslide(map,delay)
%              fktimeslide(map,delay,dblim)
%              fktimeslide(map,delay,dblim,zerodb)
%              fktimeslide(map,delay,dblim,zerodb,fgcolor,bgcolor)
%              fktimeslide(map,delay,dblim,zerodb,fgcolor,bgcolor,ax)
%              mov=fktimeslide(...);
%
%    Description: FKTIMESLIDE(MAP) slides through an array of fk maps
%     produced by FKMAP or FKVOLUME/FK4D+FKVOL2MAP by plotting each
%     map sequentially in a single plot.  There is a 1/3 second delay
%     between each replotting by default (see next Usage form to adjust
%     this).
%
%     FKTIMESLIDE(MAP,DELAY) specifies the delay between the plotting of
%     each map in seconds.  The default DELAY is 0.33s.
%
%     FKTIMESLIDE(MAP,DELAY,DBLIM) sets the dB limits for coloring the
%     response info.  The default is [-12 0] for the default ZERODB (see
%     next Usage form).  If ZERODB IS 'min' or 'median', the default DBLIM
%     is [0 12].  DBLIM must be a real-valued 2-element vector.
%
%     FKTIMESLIDE(MAP,DELAY,DBLIM,ZERODB) changes what 0dB corresponds to
%     in the plot.  The allowed values are 'min', 'max', 'median', & 'abs'.
%     The default is 'max'.
%
%     FKTIMESLIDE(MAP,DELAY,DBLIM,ZERODB,FGCOLOR,BGCOLOR) specifies the
%     foreground and background colors of the plot.  The default is 'w' for
%     FGCOLOR and 'k' for BGCOLOR.  Note that if one is specified and the
%     other is not, an opposing color is found using INVERTCOLOR.  The
%     color scale is also changed so the noise clip is at BGCOLOR.
%
%     FKTIMESLIDE(MAP,DELAY,DBLIM,ZERODB,FGCOLOR,BGCOLOR,AX) sets the axes
%     that the map is drawn in.  This is useful for subplots, guis, etc.
%
%     MOV=FKTIMESLIDE(...) creates a Matlab movie MOV.  This can played
%     using the function MOVIE.  See MOVIE2AVI for exporting as an AVI
%     file.
%
%    Notes:
%
%    Examples:
%     Get a 4D fk dataset, average over frequencies, and play:
%      s4d=fk4d(data,1,75,50,201,[1/50 1/20]);
%      s3d=fkvol2map(s4d);
%      fktimeslide(s3d);
%
%    See also: FKVOLUME, PLOTFKMAP, UPDATEFKMAP, FKMAP, FK4D, FKVOL2MAP

%     Version History:
%        May  11, 2010 - initial version
%        May  26, 2010 - update for new plotfkmap args
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  26, 2010 at 10:00 GMT

% todo:

% check nargin
error(nargchk(1,7,nargin));

% check fk struct
error(chkfkstruct(map));

% don't allow array/volume
if(any([map.volume]))
    error('seizmo:fktimeslide:badInput',...
        'MAP must be an array of fk response maps!');
end

% do we make the movie
makemovie=false;
if(nargout); makemovie=true; end

% check delay
if(nargin<2 || isempty(delay)); delay=0.33; end
if(~isreal(delay) || ~isscalar(delay) || delay<0)
    error('seizmo:fktimeslide:badDelay',...
        'DELAY must be a positive scalar in seconds!');
end

% make initial plot
ax=plotfkmap(map(1),varargin{:});
fh=get(ax,'parent');
if(makemovie); varargout{1}=getframe(fh); end

% loop over remaining frequencies
for i=2:numel(map)
    pause(delay);
    updatefkmap(map(i),ax);
    if(makemovie); varargout{1}(i)=getframe(fh); end
end

end
