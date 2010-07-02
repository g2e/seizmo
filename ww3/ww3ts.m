function [x,t,lat,lon]=ww3ts(file,latrng,lonrng)
%WW3TS    Extract a timeseries of values from WaveWatch III files
%
%    Usage:    [x,t,lat,lon]=ww3ts(files,latrng,lonrng)
%
%    Description: [X,T,LAT,LON]=WW3TS(FILES,LATRNG,LONRNG) extracts a time
%     series of values for points in the WaveWatch III files FILES within
%     the range specified by LATRNG & LONRNG.  LATRNG & LONRNG are [] by
%     default which does not exclude any points.  
%
%    Notes:
%     - LATRNG should be in -90 to 90
%     - LONRNG should be in 0 to 360
%
%    Examples:
%     Calling with no files will present a prompt to select files
%     graphically:
%      [x,t,lat,lon]=ww3ts([],[0 4],[6 10]);
%
%    See also: READ_GRIB, WW3MOV, PLOTWW3, WW3MAT

%     Version History:
%        June 30, 2010 - initial version
%        July  1, 2010 - fixed latlon bug
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July  1, 2010 at 11:00 GMT

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
        error('seizmo:ww3ts:noFileSelected','No input file selected!');
    end
    file=strcat(path,filesep,file);
    if(ischar(file)); file=cellstr(file); end
    nfiles=numel(file);
else
    % check file
    if(~ischar(file) || ~iscellstr(file))
        error('seizmo:ww3ts:fileNotString',...
            'FILES must be a string or a cell array of strings!');
    end
    if(ischar(file)); file=cellstr(file); end
    nfiles=numel(file);
    for i=1:nfiles
        if(~exist(file{i},'file'))
            error('seizmo:ww3ts:fileDoesNotExist',...
                'File: %s\nDoes Not Exist!',file{i});
        elseif(exist(file{i},'dir'))
            error('seizmo:ww3ts:dirConflict',...
                'File: %s\nIs A Directory!',file{i});
        end
    end
end

% check lat/lon ranges
if(nargin<2 || isempty(latrng)); latrng=[-90 90]; end
if(nargin<3 || isempty(lonrng)); lonrng=[0 360]; end
if(~isreal(latrng) || ~isequal(size(latrng),[1 2]) || ~issorted(latrng))
    error('seizmo:ww3ts:badLatRange',...
        'LATRNG must be a real-valued vector of [MINLAT MAXLAT]!');
end
if(~isreal(lonrng) || ~isequal(size(lonrng),[1 2]) || ~issorted(lonrng))
    error('seizmo:ww3ts:badLonRange',...
        'LONRNG must be a real-valued vector of [MINLON MAXLON]!');
end

% loop over all files
nrecs=zeros(nfiles,1);
cnt=1;
for i=1:nfiles
    % read in grib file headers
    ww3=read_grib(file{i},-1,'screendiag',0,'dataflag',0);
    nrecs(i)=numel(ww3);
    
    % error if nothing
    if(~nrecs(i))
        error('seizmo:ww3ts:badGRIB',...
            'WW3 GRIB file has no records or is not a GRIB file!');
    end
    
    % now loop over records and extract
    for j=1:nrecs(i)
        % read in data
        [map,time,lat,lon]=ww3mat(read_grib(file{i},j,...
            'screendiag',0));
        
        % extract lat/lon range
        latidx=(lat>=latrng(1) & lat<=latrng(2));
        lonidx=(lon>=lonrng(1) & lon<=lonrng(2));
        map=map(latidx,lonidx);
        lat=lat(latidx);
        lon=lon(lonidx);
        
        % push into column vector array
        x(cnt,:)=map(:);
        t(cnt,:)=time;
        cnt=cnt+1;
    end
end

% lat/lon for each point
[lon,lat]=meshgrid(lon,lat);
lat=lat(:);
lon=lon(:);

% sort by times
t=datenum(t);
[t,idx]=sort(t);
x=x(idx,:);

end
