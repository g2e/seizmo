function [varargout]=geofkframeslide(map,delay,varargin)
%GEOFKFRAMESLIDE    Slides through a set of geofk maps plotting each one
%
%    Usage:    geofkframeslide(map)
%              geofkframeslide(map,delay)
%              geofkframeslide(map,delay,dblim)
%              geofkframeslide(map,delay,dblim,zerodb)
%              geofkframeslide(map,delay,dblim,zerodb,...
%                              'mmap_opt1',mmap_val1,...)
%              mov=geofkframeslide(...);
%
%    Description:
%     GEOFKFRAMESLIDE(MAP) slides through an array of geofk maps produced
%     by geofk functions + FKVOL2MAP by plotting each map sequentially in a
%     single plot.  There is a 1/3 second delay between each replotting by
%     default (see next Usage form to adjust this).
%
%     GEOFKFRAMESLIDE(MAP,DELAY) specifies the delay between the plotting
%     of each map in seconds.  The default DELAY is 0.33s.
%
%     GEOFKFRAMESLIDE(MAP,DELAY,DBLIM) sets the dB limits for
%     coloring the response info.  The default is [-12 0] for the default
%     ZERODB (see next Usage form).  If ZERODB IS 'min' or 'median', the
%     default DBLIM is [0 12].  DBLIM must be a real-valued 2-element
%     vector.
%
%     GEOFKFRAMESLIDE(MAP,DELAY,DBLIM,ZERODB) changes what 0dB
%     corresponds to in the plot.  The allowed values are 'min', 'max',
%     'median', & 'abs'.  The default is 'max'.
%
%     GEOFKFRAMESLIDE(MAP,DELAY,DBLIM,ZERODB,'MMAP_OPT1',MMAP_VAL1,...)
%     passes additional options on to MMAP to alter the map.
%
%     MOV=GEOFKFRAMESLIDE(...) creates a Matlab movie MOV.  This can played
%     using the function MOVIE.  See MOVIE2AVI for exporting as an AVI
%     file.
%
%    Notes:
%
%    Examples:
%
%    See also: GEOFKXCVOLUME, GEOFKXCHORZVOLUME, GEOFKVOL2MAP, FKFRAMESLIDE
%              GEOFKFREQSLIDE, GEOFKSLOWSLIDE, PLOTGEOFKMAP, UPDATEGEOFKMAP

%     Version History:
%        July  7, 2010 - initial version
%        Apr.  4, 2012 - minor doc update
%        Mar.  5, 2014 - update doc for plotgeofkmap input changes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  5, 2014 at 23:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check fk struct
error(chkgeofkstruct(map));

% don't allow array/volume
if(any([map.volume]))
    error('seizmo:geofkframeslide:badInput',...
        'MAP must be an array of geofk beam maps!');
end

% do we make the movie
makemovie=false;
if(nargout); makemovie=true; end

% check delay
if(nargin<2 || isempty(delay)); delay=0.33; end
if(~isreal(delay) || ~isscalar(delay) || delay<0)
    error('seizmo:geofkframeslide:badDelay',...
        'DELAY must be a positive scalar in seconds!');
end

% make initial plot
ax=plotgeofkmap(map(1),varargin{:});
fh=get(ax,'parent');
if(makemovie); varargout{1}=getframe(fh); else drawnow; end

% loop over remaining frequencies
for i=2:numel(map)
    pause(delay);
    updategeofkmap(map(i),ax);
    if(makemovie); varargout{1}(i)=getframe(fh); else drawnow; end
end

end
