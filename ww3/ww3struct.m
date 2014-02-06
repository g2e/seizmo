function [s]=ww3struct(file,rec,stime,etime,latrng,lonrng)
%WW3STRUCT    Puts WaveWatch III hindcast GRiB1/2 files into a structure
%
%    Usage:    s=ww3struct('file')
%              s=ww3struct({'file1' ... 'fileN'})
%              s=ww3struct()
%              s=ww3struct('file',rec)
%              s=ww3struct('file',rec,stime,etime,latrng,lonrng)
%
%    Description:
%     S=WW3STRUCT('FILE') reads the WW3 hindcast GRiB or GRiB2 file FILE,
%     returning the contents as a struct S with the following layout:
%       .path        - path to the WW3 hindcast file
%       .name        - name of the WW3 hindcast file
%       .description - description of the data
%       .units       - units of the data
%       .data        - data values as LATxLONxTIME
%       .lat         - latitudes (degrees)
%       .lon         - longitudes (degrees)
%       .time        - times (in serial number format aka day numer)
%       .latstep     - latitude grid step size (degrees)
%       .lonstep     - longitude grid step size (degrees)
%       .timestep    - time step size (in fraction of a day)
%     Note that the fields .description, .units & .data are returned as
%     cell arrays so that WW3 hindcast files with multiple variables are
%     supported (this allows reading wind u/v data).
%
%     S=WW3STRUCT({'FILE1' ... 'FILEN'}) processes multiple WW3 hindcast
%     files.  The output structure S is Nx1 where N is the number of valid
%     hindcast files.
%
%     S=WW3STRUCT() presents a GUI for the user to select WW3 hindcast
%     file(s).
%
%     S=WW3STRUCT('FILE',REC) only reads the record numbers specified in
%     REC.  This is useful for working time step by time step.
%
%     S=WW3STRUCT('FILE',REC,STIME,ETIME,LATRNG,LONRNG) limits output to
%     the ranges specified.  The defaults are [], which do not limit the
%     output.
%
%    Notes:
%     - Requires that the njtbx toolbox is installed!
%     - STIME & ETIME must be in one of the following formats:
%        1x1 -> SERIAL (See datenum)
%        1x2 -> [YR DOY]
%        1x3 -> [YR MO CDAY]
%        1x5 -> [YR DOY HR MIN SEC]
%        1x6 -> [YR MO CDAY HR MIN SEC]
%
%    Examples:
%     % Read the first record of a NOAA WW3 grib file and plot it up:
%     s=ww3struct('nww3.hs.200607.grb',1);
%     plotww3(s);
%
%    See also: WW3REC, WW3CAT, PLOTWW3, PLOTWW3TS, WW3MOV, WW3MAP,
%              WW3MAPMOV, WW3UV2SA, WW3BAZ2AZ

%     Version History:
%        June 30, 2010 - initial version
%        July  2, 2010 - added WW3TS to see also
%        Dec.  5, 2011 - doc update
%        Feb. 14, 2012 - njtbx instead of read_grib, range limits, file
%                        input, multi-file support, multi-field support,
%                        renamed from ww3mat to ww3struct
%        Feb. 27, 2012 - minor doc update
%        May   5, 2012 - minor doc update
%        May  11, 2012 - skips .data access if there is a null range
%        Jan. 15, 2014 - updated See also list
%        Jan. 26, 2014 - fixed error msg for no njtbx installed
%        Feb.  5, 2014 - minor doc update, bugfix: pixel longitude
%                        correction was not needed (decided from ww3map)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  5, 2014 at 00:40 GMT

% todo:

% check nargin
error(nargchk(0,6,nargin));

% defaults
if(nargin<1); file=[]; end
if(nargin<2); rec=[]; end
if(nargin<3); stime=[]; end
if(nargin<4); etime=[]; end
if(nargin<5); latrng=[]; end
if(nargin<6); lonrng=[]; end

% check that njtbx is installed (somewhat)
if(~exist('mDataset','file') || ~exist('nj_time','file'))
    error('seizmo:ww3struct:noNJTBX',...
        'NJTBX is not installed!  Run function WEBINSTALL_NJTBX.');
end

% set filterspec appropriately
global SEIZMO
SEIZMO.ONEFILELIST.FILTERSPEC=...
    {'*.grb;*.GRB;*.grb2;*.GRB2' 'GRiB Files (*.grb,*.GRB,*.grb2,*.GRB2)'};
 
% parse file input
file=onefilelist(file);
nfiles=numel(file);

% error if no files
if(nfiles<1); error('seizmo:ww3struct:noFiles','No files to read!'); end

% check the limits
stsz=size(stime); etsz=size(etime);
if(~isreal(rec) || ~isreal(stime) || ~isreal(etime) ...
        || ~isreal(latrng) || ~isreal(lonrng))
    error('seizmo:ww3struct:badInput',...
        'Range inputs must be real-valued!');
elseif(~isempty(rec) && (~isvector(rec) || any(fix(rec)~=rec)))
    error('seizmo:ww3struct:badInput',...
        'REC must be a vector of integers!');
elseif(numel(stsz)>2 || stsz(1)>1 || ~any(stsz(2)==[0 1 2 3 5 6]))
    error('seizmo:ww3struct:badInput',...
        'STIME must be a 1x? numeric vector giving a time!');
elseif(numel(etsz)>2 || etsz(1)>1 || ~any(etsz(2)==[0 1 2 3 5 6]))
    error('seizmo:ww3struct:badInput',...
        'ETIME must be a 1x? numeric vector giving a time!');
elseif(~isempty(latrng) && (~isequal(size(latrng),[1 2]) ...
        || latrng(1)>latrng(2)))
    error('seizmo:ww3struct:badInput',...
        'LATRNG must be a 1x2 vector as [MINLAT MAXLAT]!');
elseif(~isempty(lonrng) && (~isequal(size(lonrng),[1 2]) ...
        || lonrng(1)>lonrng(2)))
    error('seizmo:ww3struct:badInput',...
        'LONRNG must be a 1x2 vector as [MINLON MAXLON]!');
end

% convert times to serial if needed
if(numel(stime)>1); stime=gregorian2serial(stime); end
if(numel(etime)>1); etime=gregorian2serial(etime); end

% pre-allocate structure
s(nfiles,1)=struct('path',[],'name',[],'description',[],'units',[],...
    'data',[],'lat',[],'lon',[],'time',[]);

% loop over each file
bad=false(nfiles,1);
for i=1:nfiles
    % load the file
    gh=mDataset([file(i).path file(i).name]);
    if(isempty(gh)); bad(i)=true; continue; end
    
    % get field names
    vf={'latLonCoordSys' 'time' 'lat' 'lon'};
    f=setdiff(getVars(gh),vf);
    
    % add supporting info
    % - we assume all fields have the same time/lat/lon
    s(i).path=file(i).path;
    s(i).name=file(i).name;
    s(i).lat=gh{'lat'}(:);
    s(i).lon=gh{'lon'}(:);
    s(i).time=nj_time(gh,'time');
    s(i).latstep=s(i).lat(2)-s(i).lat(1);
    s(i).lonstep=s(i).lon(2)-s(i).lon(1);
    s(i).timestep=s(i).time(2)-s(i).time(1);
    
    % handle no record limiting
    if(isempty(rec))
        newrec=1:numel(s(i).time);
    else
        newrec=rec;
    end
    
    % limit output
    if(~isempty(latrng))
        oklat=s(i).lat>=latrng(1) & s(i).lat<=latrng(2);
        s(i).lat=s(i).lat(oklat);
    else
        oklat=true(size(s(i).lat));
    end
    if(~isempty(lonrng))
        oklon=(s(i).lon>=lonrng(1) & s(i).lon<=lonrng(2)) ...
            | (s(i).lon+360>=lonrng(1) & s(i).lon+360<=lonrng(2)) ...
            | (s(i).lon-360>=lonrng(1) & s(i).lon-360<=lonrng(2));
        s(i).lon=s(i).lon(oklon);
    else
        oklon=true(size(s(i).lon));
    end
    oktime=false(size(s(i).time));
    oktime(newrec)=true;
    if(~isempty(stime))
        oktime=oktime & s(i).time>=stime;
    end
    if(~isempty(etime))
        oktime=oktime & s(i).time<=etime;
    end
    s(i).time=s(i).time(oktime);
    
    % convert to indices
    oklat=find(oklat); % logical indexing has
    oklon=find(oklon); % issues so use linear
    oktime=find(oktime);
    
    % size of output
    ntime=numel(s(i).time);
    nlat=numel(s(i).lat);
    nlon=numel(s(i).lon);
    
    % loop over fields
    for j=1:numel(f)
        % field specific info
        s(i).description{j}=strrep(f{j},'_',' ');
        s(i).units{j}=gh{f{j}}.units;
        
        % preallocate data output and loop over each time
        % - this tries to avoid java memory issues
        s(i).data{j}=nan(nlat,nlon,ntime);
        if(nlat*nlon*ntime>0)
            for k=1:ntime
                s(i).data{j}(:,:,k)=gh{f{j}}(...
                    oktime(k),oklat,oklon).data; %#ok<FNDSB>
            end
        end
    end
    
    % close file
    close(gh);
end

end
