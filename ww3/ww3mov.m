function [varargout]=ww3mov(file,delay,varargin)
%WW3MOV    Create movie from a WaveWatch III hindcast GRiB1/2 file
%
%    Usage:    mov=ww3mov('file')
%              [mov1,...,movN]=ww3mov('file')
%              mov=ww3mov('file',delay)
%              mov=ww3mov('file',delay,rng)
%              mov=ww3mov('file',delay,rng,fgcolor,bgcolor)
%              mov=ww3mov('file',delay,rng,fgcolor,bgcolor,ax)
%
%    Description:
%     MOV=WW3MOV('FILE') creates a Matlab movie of the WaveWatch III
%     handcast data contained in the GRiB file FILE.  MOV is the movie
%     struct that can then be converted to an AVI file using MOVIE2AVI.
%     There is a 1/3 second delay between each frame by default (see next
%     Usage form to adjust this).  If FILE is omitted or is empty then a
%     GUI is presented for GRiB file selection.  If no output is assigned
%     then WW3MOV will "play" the data.
%
%     [MOV1,...,MOVN]=WW3MOV('FILE') returns the movies for each data type
%     in FILE (eg for wind data there is a movie for each component).
%
%     MOV=WW3MOV('FILE',DELAY) specifies the delay between the plotting of
%     each time step in seconds.  The default DELAY is 0.33s.
%
%     MOV=WW3MOV('FILE',DELAY,RNG) sets the limits for coloring the data.
%     The default is [0 15] which works well for significant wave heights.
%
%     MOV=WW3MOV('FILE',DELAY,RNG,FGCOLOR,BGCOLOR) specifies foreground and
%     background colors of the movie.  The default is 'w' for FGCOLOR & 'k'
%     for BGCOLOR.  Note that if one is specified and the other is not, an
%     opposing color is found using INVERTCOLOR.  The color scale is also
%     changed so the noise clip is at BGCOLOR.
%
%     MOV=WW3MOV('FILE',DELAY,RNG,FGCOLOR,BGCOLOR,AX) sets the axes drawn
%     in.  This is useful for subplots, guis, etc.  The default creates a
%     new figure.
%
%    Notes:
%     - Requires that the njtbx toolbox is installed!
%
%    Examples:
%     % Calling WW3MOV with no args lets you graphically choose a file:
%     mov=ww3mov;
%
%    See also: PLOTWW3, MOVIE2AVI, UNIXCOMPRESSAVI, WW3STRUCT

%     Version History:
%        June 15, 2010 - initial version
%        July  2, 2010 - adjusted oneliner description
%        Aug. 30, 2010 - fix documentation
%        Feb. 15, 2012 - use ww3struct, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 15, 2012 at 00:40 GMT

% todo:

% check nargin
error(nargchk(0,6,nargin));

% attempt reading in first record of file
% - this does the gui & checks file is valid
s=ww3struct(file,1);

% get number of time steps
h=mDataset(fullfile(s.path,s.name));
nrecs=numel(h{'time'}(:));
close(h);

% check delay
if(nargin<2 || isempty(delay)); delay=0.33; end
if(~isreal(delay) || ~isscalar(delay) || delay<0)
    error('seizmo:ww3mov:badDelay',...
        'DELAY must be a positive scalar in seconds!');
end

% only make movie if output
makemovie=false;
if(nargout); makemovie=true; end

% make initial plot
ax=plotww3(s,varargin{:});
varargin{4}=ax;
fh=get(ax,'parent');
if(iscell(fh)); fh=cell2mat(fh); end
for j=1:numel(fh)
    if(makemovie); varargout{j}=getframe(fh(j)); end
end

% now loop over records and plot them
for i=2:nrecs
    pause(delay);
    s=ww3struct(file,i);
    if(any(~ishghandle(ax,'axes')))
        error('seizmo:ww3mov:userClose',...
            'Axes disappeared! Did someone turn off the lights?');
    end
    plotww3(s,varargin{:});
    drawnow;
    for j=1:numel(fh)
        if(makemovie); varargout{j}(i)=getframe(fh(j)); end
    end
end

end
