function [varargout]=plotww3avg(file,varargin)
%PLOTWW3AVG    Plots an average WaveWatch III data map
%
%    Usage:    plotww3avg(files)
%              plotww3avg(ww3,rng)
%              plotww3avg(ww3,rng,fgcolor,bgcolor)
%              plotww3avg(ww3,rng,fgcolor,bgcolor,ax)
%              ax=plotww3avg(...)
%
%    Description:
%     PLOTWW3AVG(FILES) takes the average of all the WaveWatch III data
%     contained in the GRIB files given by FILES and plots it.  This is for
%     making monthly, seasonal, or longer averages.
%
%     PLOTWW3AVG(WW3,RNG) sets the limits for coloring the data. The
%     default is [0 15] which works well for significant wave height.
%
%     PLOTWW3AVG(WW3,RNG,FGCOLOR,BGCOLOR) specifies foreground and
%     background colors of the plot.  The default is 'w' for FGCOLOR & 'k'
%     for BGCOLOR.  Note that if one is specified and the other is not, an
%     opposing color is found using INVERTCOLOR.  The color scale is also
%     changed so the noise clip is at BGCOLOR.
%
%     PLOTWW3AVG(WW3,RNG,FGCOLOR,BGCOLOR,AX) sets the axes to draw in. This
%     is useful for subplots, guis, etc.
%
%     AX=PLOTWW3AVG(...) returns the axes drawn in.  This is useful for
%     after-touching.
%
%    Notes:
%
%    Examples:
%     % Calling with no files will present a prompt to select files
%     % graphically:
%     ax=plotww3avg();
%
%    See also: READ_GRIB, WW3MOV, PLOTWW3, WW3MAT

%     Version History:
%        Aug. 30, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 30, 2010 at 11:00 GMT

% todo:

% check nargin
error(nargchk(0,3,nargin));

% graphical selection
if(nargin<1 || isempty(file))
    [file,path]=uigetfile(...
        {'*.grb;*.grib' 'GRIB Files (*.grb,*.grib)';
        '*.*' 'All Files (*.*)'},...
        'Select GRIB File',...
        'multiselect','on');
    if(isequal(0,file) || isequal(0,path))
        error('seizmo:plotww3avg:noFileSelected',...
            'No input file selected!');
    end
    file=strcat(path,filesep,file);
    if(ischar(file)); file=cellstr(file); end
    nfiles=numel(file);
else
    % check file
    if(~ischar(file) || ~iscellstr(file))
        error('seizmo:plotww3avg:fileNotString',...
            'FILES must be a string or a cell array of strings!');
    end
    if(ischar(file)); file=cellstr(file); end
    nfiles=numel(file);
    for i=1:nfiles
        if(~exist(file{i},'file'))
            error('seizmo:plotww3avg:fileDoesNotExist',...
                'File: %s\nDoes Not Exist!',file{i});
        elseif(exist(file{i},'dir'))
            error('seizmo:plotww3avg:dirConflict',...
                'File: %s\nIs A Directory!',file{i});
        end
    end
end

% read in first record of first file and set to zero
ww3=read_grib(file{1},1,'screendiag',0);
try
    ww3.fltarray=ww3.fltarray.*0;
catch
    error('seizmo:plotww3avg:badGRIB',...
        'WW3 GRIB file has no records or is not a GRIB file!');
end

% loop over all files
nrecs=zeros(nfiles,1);
mint=inf; maxt=-inf; cnt=0;
for i=1:nfiles
    % read in grib file headers
    tmp=read_grib(file{i},-1,'screendiag',0,'dataflag',0);
    nrecs(i)=numel(tmp);
    
    % error if nothing
    if(~nrecs(i))
        error('seizmo:plotww3avg:badGRIB',...
            'WW3 GRIB file has no records or is not a GRIB file!');
    end
    
    % now loop over records and add them
    for j=1:nrecs(i)
        % increment counter (for averaging)
        cnt=cnt+1;
        
        % read in data
        tmp=read_grib(file{i},j,'screendiag',0);
        
        % add in
        ww3.fltarray=ww3.fltarray+tmp.fltarray;
        time=[tmp.pds.year tmp.pds.month tmp.pds.day ...
            tmp.pds.hour tmp.pds.min 0];
        time=datenum(time);
        if(time<mint)
            smint=tmp.stime;
            mint=time;
        end
        if(time>maxt)
            smaxt=tmp.stime;
            maxt=time;
        end
    end
end

% now plot average
ww3.fltarray=ww3.fltarray./cnt;
ax=plotww3(ww3,varargin{:});

% fix title
set(get(ax,'title'),'string',...
    {'NOAA WaveWatch III Hindcast' ww3.description [smint ' to ' smaxt]});

% handle output
if(nargout); varargout{1}=ax; end

end
