function [events]=readsodeventcsv(varargin)
%READSODEVENTCSV    Read in a SOD event .csv formatted file as a structure
%
%    Usage:    events=readsodeventcsv(file)
%
%    Description:
%     EVENTS=READSODEVENTCSV(FILE) reads in a SOD (Standing Order for Data)
%     generated event CSV file as a scalar structure array.  The CSV (comma
%     separated values) file is expected to have 1 header line of comma
%     separated field names.  These define the fields of the structure
%     which are returned as NLx1 arrays where NL is the number of lines in
%     the file excluding the header line.  The function SSIDX is useful for
%     accessing specific entries of scalar structs (see the Examples
%     section below).  Calling READSODEVENTCSV without a FILE argument or
%     with FILE set to '' will present a graphical file selection menu.
%
%    Notes:
%     - Converts latitude, longitude, depth and magnitude to numeric form.
%     - Converts time from string to [yr mon cday hr min secs].
%     - Fields of a standard SOD Event CSV file:
%           time
%           latitude
%           longitude
%           depth
%           depthUnits
%           magnitude
%           magnitudeType
%           catalog
%           contributor
%           name
%           flinnEngdahlRegion
%           flinnEngdahlRegionType
%
%     - SOD (Standing Order for Data) is written/maintained by
%       Philip Crotwell.  Website: http://www.seis.sc.edu/sod/
%
%    Examples:
%     % READSODEVENTCSV uses READCSV but also converts several
%     % fields into a more useful format.  To show the differences:
%     ev1=readcsv(file)
%     ev2=readsodeventcsv(file)
%
%    See also: WRITESODEVENTCSV, READCSV, WRITECSV, SETEVENT, SSIDX, SSCAT

%     Version History:
%        Sep. 16, 2009 - initial version
%        Jan. 26, 2010 - allow no input (select file graphically)
%        Mar. 30, 2010 - check fields are available to modify, doc update
%        Feb. 11, 2011 - mass nargchk fix
%        Jan. 28, 2012 - doc update, pass char to strnlen
%        Feb. 28, 2012 - update for READCSV scalar struct output
%        Jan. 17, 2013 - fix example, add ssidx/sscat to see also
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 17, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% send input on to readcsv
events=readcsv(varargin{:});

% require certain fields are present
req={'time' 'latitude' 'longitude' 'depth' 'magnitude'};
fields=fieldnames(events);
if(~all(ismember(req,fields)))
    error('seizmo:readsodeventcsv:missingFields',...
        ['CSV file must have the following fields:\n' ...
        sprintf('%s ',req{:})]);
end

% convert some fields to numbers
f={'latitude' 'longitude' 'depth' 'magnitude'};
for i=1:numel(f)
    events.(f{i})=str2double(events.(f{i}));
end

% clean up depths
events.depth=round(100*events.depth)/100;

% convert time to something useful
% - note that we strip off the Z first b/c there is
%   a bug in dtstr2dtvecmx (datevec) that occasionally
%   returns an incorrect time when the Z is included
events.time=datevec(strnlen(char(events.time),23),...
    'yyyy-mm-ddTHH:MM:SS.FFF');

end
