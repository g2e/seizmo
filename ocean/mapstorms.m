function [varargout]=mapstorms(storms,lineopt,varargin)
%MAPSTORMS    Maps storm track data
%
%    Usage:    mapstorms(storms)
%              mapstorms(storms,lineopt)
%              mapstorms(storms,lineopt,'opt1',val1,'opt2',val2,...)
%              ax=mapstorms(...)
%
%    Description:
%     MAPSTORMS(STORMS) plots the storms tracks in struct STORMS.
%
%     MAPSTORMS(STORMS,LINEOPT) passes a cell array of options to M_LINE
%     which will then pass them to LINE.  See the Examples section for
%     more.
%
%     MAPSTORMS(STORMS,LINEOPT,'OPT1',VAL1,'OPT2',VAL2,...) allows passing
%     options for mapping.  See MMAP for details.
%
%     AX=MAPSTORMS(...) returns the axis handle AX of the plot.
%
%    Notes:
%
%    Examples:
%     % First get the HURDAT "all storms" data from here:
%     %  http://www.ncdc.noaa.gov/oa/ibtracs/index.php?name=wmo-data
%     % Now read in the file and map the tracks:
%     storms=read_hurdat('AllStorms.ibtracs_hurdat.v03r04.hdat',true);
%     mapstorms(storms);
%
%     % Map tracks in red & with linewidth 2:
%     mapstorms(storms,{'color','r','linewidth',2});
%
%    See also: READ_HURDAT, READ_GISS_STORMDB, MMAP, MAPFEATURE

%     Version History:
%        Feb. 16, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 16, 2013 at 13:30 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin))
if(nargin>2 && mod(nargin,2))
    error('seizmo:mapstorms:badNumInputs',...
        'Unpaired Option/Value!');
end

% check storms struct
error(chkstorms(storms));

% check lineopt
if(nargin<2 || isempty(lineopt)); lineopt={}; end
if(~iscell(lineopt))
    error('seizmo:mapstorms:badNumInputs',...
        'LINEOPT must be a cell array!');
end

% call mmap
ax=mmap(varargin{:});

% plot the storms
% - simple plotting for now
%   - would like to allow coloring by wind/pressure/category/stage
for i=1:numel(storms.lat)
    m_line(unwrap(storms.lon{i}*(pi/180))/(pi/180),storms.lat{i},...
        'tag','stormtracks','parent',ax,lineopt{:});
end

% optional output
if(nargout); varargout{1}=ax; end

end
