function [varargout]=fkfreqslide(vol,delay,varargin)
%FKFREQSLIDE    Slides through a fk volume plotting each frequency
%
%    Usage:    fkfreqslide(vol)
%              fkfreqslide(vol,delay)
%              fkfreqslide(vol,delay,fgcolor,bgcolor)
%              fkfreqslide(vol,delay,fgcolor,bgcolor,ax)
%              mov=fkfreqslide(...);
%
%    Description: FKFREQSLIDE(VOL) slides through a fk volume produced by
%     FKVOLUME by plotting each of the frequencies in sequence in a single
%     plot.  There is a 1/3 second delay between each replotting by default
%     (see next Usage form to adjust this).
%
%     FKFREQSLIDE(VOL,DELAY) specifies the delay between the plotting of
%     each frequency in seconds.  The default DELAY is 0.33s.
%
%     FKFREQSLIDE(VOL,DELAY,FGCOLOR,BGCOLOR) specifies the foreground and
%     background colors of the plot.  The default is 'w' for FGCOLOR and
%     'k' for BGCOLOR.  Note that if one is specified and the other is not,
%     an opposing color is found using INVERTCOLOR.  The color scale is
%     also changed so the noise clip is at BGCOLOR.
%
%     FKFREQSLIDE(VOL,DELAY,FGCOLOR,BGCOLOR,AX) sets the axes that the map
%     is drawn in.  This is useful for subplots, guis, etc.
%
%     MOV=FKFREQSLIDE(...) creates a Matlab movie MOV.  This can played
%     using the function MOVIE.  See MOVIE2AVI for exporting as an AVI
%     file.
%
%    Notes:
%
%    Examples:
%     Show frequency-slowness volume for a dataset at 20-50s periods:
%      svol=fkvolume(data,50,201,[1/50 1/20]);
%      fkfreqslide(svol);
%
%    See also: FKVOLUME, PLOTFKMAP, UPDATEFKMAP, FKMAP, FK4D

%     Version History:
%        May  11, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  11, 2010 at 16:20 GMT

% todo:

% check nargin
error(nargchk(1,5,nargin));

% check fk struct
error(chkfkstruct(vol));

% don't allow array/volume
if(~isscalar(vol) || ~vol.volume)
    error('seizmo:fkfreqslide:badInput',...
        'VOL must be a scalar fk struct and a volume response!');
end

% do we make the movie
makemovie=false;
if(nargout); makemovie=true; end

% check delay
if(nargin<2 || isempty(delay)); delay=0.33; end
if(~isreal(delay) || ~isscalar(delay) || delay<0)
    error('seizmo:fkfreqslide:badDelay',...
        'DELAY must be a positive scalar in seconds!');
end

% get frequencies
freqs=vol.z;
nfreq=numel(freqs);

% make initial plot
ax=plotfkmap(fkvol2map(vol,[freqs(1) freqs(1)]),varargin{:});
fh=get(ax,'parent');
if(makemovie); varargout{1}=getframe(fh); end

% loop over remaining frequencies
for i=2:nfreq
    pause(delay);
    updatefkmap(fkvol2map(vol,[freqs(i) freqs(i)]),ax);
    if(makemovie); varargout{1}(i)=getframe(fh); end
end

end
