function [storms]=read_hurdat(file,allstorms,latmult,lonmult)
%READ_HURDAT    Reads hurricane storm track data in the HURDAT format
%
%    Usage:    storms=read_hurdat(file)
%              storms=read_hurdat(file,allstormsflag)
%              storms=read_hurdat(file,[],latmult,lonmult)
%
%    Description:
%     STORMS=READ_HURDAT(FILE) reads storm track data stored in a HURDAT
%     formatted ascii file.  There are minor variations in the format so
%     you will probably need to use the flags or multipliers to properly
%     read the data in.  The output STORMS is a scalar struct with the
%     following fields:
%       STORMS.name             -- name of storm
%             .time             -- times of measurements (datenum format)
%             .stage            -- storm stage code (1 char code)
%             .lat              -- latitude (degrees North)
%             .lon              -- longitude (degrees East)
%             .wind             -- wind in knots
%             .pressure         -- pressure in mbar
%             .impact_us        -- impacted US coast?
%             .impact_category  -- hurricane category at landfall
%             .max_stage        -- maximum storm stage code (2 letter code)
%     Where .time, .stage, .lat, .lon, .wind, & .pressure are the storm
%     track data for each storm.  Note that for the ith storm the data is
%     in the ith element of each field.  See SSIDX & SSCAT for help with
%     working with 
%
%     STORMS=READ_HURDAT(FILE,ALLSTORMSFLAG) is useful for handling IBTrACS
%     files with all tropical storms.  These are formatted slightly
%     different to handle latitude and longitude in all hemispheres.  Set
%     ALLSTORMSFLAG to TRUE to handle this modification.  The default is
%     FALSE.
%
%     STORMS=READ_HURDAT(FILE,[],LATMULT,LONMULT) is useful for putting
%     "basin" latitudes & longitudes in their true place.  The default
%     multipliers (1) are correct for north latitudes and east longitudes.
%     Use -1 to map southern hemisphere & west hemisphere data in their
%     proper locations.
%
%    Notes:
%     - See the following for a detailed explaination of the format:
%        http://www.aoml.noaa.gov/hrd/data_sub/hurdat.html
%     - The following document explains the IBTrACS modifications:
%        ftp://eclipse.ncdc.noaa.gov/pub/ibtracs/v03r04/wmo/...
%        hurdat_format/README.hurdat_format
%
%    Examples:
%     % Read Atlantic data:
%     storms=read_hurdat('hurdat_atlantic_1851-2011.txt',[],1,-1);
%
%    See also: READ_GISS_STORMDB, MAPSTORMS

%     Version History:
%        Feb. 16, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 16, 2013 at 13:30 GMT

% todo:

% check nargin
error(nargchk(0,4,nargin));

% defaults
if(~nargin); file=[]; end
if(nargin<2 || isempty(allstorms)); allstorms=false; end
if(nargin<3 || isempty(latmult)); latmult=1; end
if(nargin<4 || isempty(lonmult)); lonmult=1; end

% filterspec
filterspec={'*.*' 'All Files (*.*)'};

% use readtxt to get text of file
txt=readtxt(file,filterspec);

% check optional inputs
if(~isscalar(allstorms) || ~islogical(allstorms))
    error('seizmo:read_hurdat:badInput',...
        'ALLSTORMS flag must be TRUE or FALSE!');
elseif(~isscalar(latmult) || abs(latmult)~=1)
    error('seizmo:read_hurdat:badInput',...
        'LATMULT must be 1 or -1!');
elseif(~isscalar(lonmult) || abs(lonmult)~=1)
    error('seizmo:read_hurdat:badInput',...
        'LONMULT must be 1 or -1!');
end

% delete carriage return characters
txt(txt==13)=[];

% split lines
lines=getwords(txt,sprintf('\n'));
nlines=numel(lines);

% now lets count the storms and how long they are
[ln,nd]=deal(nan(1e6,1)); % allocating by guess
line=1; nstorm=0; bad=false(nlines,1);
while(line<nlines)
    if(isempty(lines{line})); bad(line)=true; line=line+1; continue; end
    nstorm=nstorm+1;
    ln(nstorm)=line; % starting line of each storm
    nd(nstorm)=str2double(lines{line}(20:21)); % number of days for storm
    line=line+nd(nstorm)+2; % jump past the storm track data & trailer
end

 % deallocate unfilled
ln(nstorm+1:end)=[];
nd(nstorm+1:end)=[];

% allocate output
storms=struct('name',[],'time',[],'stage',[],'lat',[],...
    'lon',[],'wind',[],'pressure',[],'impact_us',[],...
    'impact_category',[],'max_stage',[]);

% split up
head=char(lines(ln));
tail=char(lines(ln+nd+1));
bad([ln; ln+nd+1])=true;
data=char(lines(~bad));
clear lines;

% header info
storms.name=cellstr(head(:,36:46));
storms.time=mat2cell(datenum(...
    [str2num(head(:,13:16)) str2num(head(:,7:8)) ...
    str2num(head(:,10:11)) zeros(nstorm,3)]),ones(nstorm,1));
storms.impact_us=strcmp(cellstr(head(:,53)),'1');
storms.impact_category=str2num(head(:,59));

% trailer info
storms.max_stage=cellstr(tail(:,7:8));
%lcat=tail(:,12:4:end); lcat=reshape(str2num(lcat(:)),nstorm,[]);
%tail(:,[1:8 12:4:end])=[];
%lpl=mat2cell(tail,ones(size(tail,1),1),3*ones(1,size(tail,2)/3));

% storm track data
tmp=data(:,12:17:end-1)';
storms.stage=mat2cell(tmp(:),4*nd);
tmp=str2double(mat2cell(data(:,sort([13:17:79 14:17:79 15:17:79])),...
    ones(sum(nd),1),3*ones(1,4)))'/10;
if(allstorms); tmp=tmp.*(1-2*(data(:,20:17:79)=='s')'); end
storms.lat=mat2cell(latmult*tmp(:),4*nd);
tmp=str2double(mat2cell(data(:,...
    sort([16:17:79 17:17:79 18:17:79 19:17:79])),...
    ones(sum(nd),1),4*ones(1,4)))'/10;
storms.lon=mat2cell(lonmult*tmp(:),4*nd);
tmp=str2double(mat2cell(data(:,sort([21:17:79 22:17:79 23:17:79])),...
    ones(sum(nd),1),3*ones(1,4)))'/10;
storms.wind=mat2cell(tmp(:),4*nd);
tmp=str2double(mat2cell(data(:,...
    sort([25:17:79 26:17:79 27:17:79 28:17:79])),...
    ones(sum(nd),1),4*ones(1,4)))'/10;
storms.pressure=mat2cell(tmp(:),4*nd);

% storm track times & clear null storm track info
for i=1:nstorm
    storms.time{i}=storms.time{i}+(((1:nd(i)*4)-1)/4)';
    tmp=storms.lat{i}==0 & storms.lon{i}==0 & storms.wind{i}==0;
    storms.time{i}(tmp)=[];
    storms.stage{i}(tmp)=[];
    storms.lat{i}(tmp)=[];
    storms.lon{i}(tmp)=[];
    storms.wind{i}(tmp)=[];
    storms.pressure{i}(tmp)=[];
end

end
