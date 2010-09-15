function [varargout]=ww3mov(file,delay,varargin)
%WW3MOV    Create movie from a WaveWatch III file
%
%    Usage:    mov=ww3mov(gribfile)
%              mov=ww3mov(gribfile,delay)
%              mov=ww3mov(gribfile,delay,rng)
%              mov=ww3mov(gribfile,delay,rng,fgcolor,bgcolor)
%              mov=ww3mov(gribfile,delay,rng,fgcolor,bgcolor,h)
%
%    Description: MOV=WW3MOV(GRIBFILE) creates a Matlab movie of the
%     WaveWatch III data contained in GRIBFILE.  Note that this requires
%     that READ_GRIB is available and working in Matlab.  MOV is the movie
%     struct that can then be converted to an AVI file using MOVIE2AVI.
%     There is a 1/3 second delay between each frame by default (see next
%     Usage form to adjust this).
%
%     MOV=WW3MOV(WW3,DELAY) specifies the delay between the plotting of
%     each frequency in seconds.  The default DELAY is 0.33s.
%
%     MOV=WW3MOV(WW3,DELAY,RNG) sets the limits for coloring the data. The
%     default is [0 15] which works well for significant wave height.
%
%     MOV=WW3MOV(WW3,DELAY,RNG,FGCOLOR,BGCOLOR) specifies foreground and
%     background colors of the movie.  The default is 'w' for FGCOLOR & 'k'
%     for BGCOLOR.  Note that if one is specified and the other is not, an
%     opposing color is found using INVERTCOLOR.  The color scale is also
%     changed so the noise clip is at BGCOLOR.
%
%     MOV=WW3MOV(WW3,DELAY,RNG,FGCOLOR,BGCOLOR,H) sets the axes to draw in.
%     This is useful for subplots, guis, etc.  The default is a new plot.
%
%    Notes:
%
%    Examples:
%     Calling WW3MOV with no args lets you graphically choose a file:
%      mov=ww3mov();
%
%    See also: READ_GRIB, PLOTWW3, MOVIE2AVI, UNIXCOMPRESSAVI, WW3MAT,
%              WW3TS

%     Version History:
%        June 15, 2010 - initial version
%        July  2, 2010 - adjusted oneliner description
%        Aug. 30, 2010 - fix documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 30, 2010 at 00:40 GMT

% todo:

% check nargin
error(nargchk(0,6,nargin));

% graphical selection
if(nargin<1 || isempty(file))
    [file,path]=uigetfile(...
        {'*.grb;*.grib' 'GRIB Files (*.grb,*.grib)';
        '*.*' 'All Files (*.*)'},...
        'Select GRIB File');
    if(isequal(0,file))
        error('seizmo:ww3mov:noFileSelected','No input file selected!');
    end
    file=strcat(path,filesep,file);
else
    % check file
    if(~ischar(file))
        error('seizmo:ww3mov:fileNotString',...
            'FILE must be a string!');
    end
    if(~exist(file,'file'))
        error('seizmo:ww3mov:fileDoesNotExist',...
            'File: %s\nDoes Not Exist!',file);
    elseif(exist(file,'dir'))
        error('seizmo:ww3mov:dirConflict',...
            'File: %s\nIs A Directory!',file);
    end
end

% read in grib file headers
ww3=read_grib(file,-1,'screendiag',0,'dataflag',0);
nrecs=numel(ww3);

% error if nothing
if(~nrecs)
    error('seizmo:ww3mov:badGRIB',...
        'WW3 GRIB file has no records or is not a GRIB file!');
end

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
ww3=read_grib(file,1,'screendiag',0);
ax=plotww3(ww3,varargin{:});
varargin{4}=ax;
fh=get(ax,'parent');
if(makemovie); varargout{1}=getframe(fh); end

% now loop over records and plot them
for i=2:nrecs
    pause(delay);
    ww3=read_grib(file,i,'screendiag',0);
    plotww3(ww3,varargin{:});
    if(makemovie); varargout{1}(i)=getframe(fh); end
end

end
