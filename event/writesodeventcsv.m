function []=writesodeventcsv(file,events,varargin)
%WRITESODEVENTCSV    Write out a SOD event .csv file from a structure
%
%    Usage:    writesodeventcsv(file,events)
%              writesodeventcsv(file,events,overwrite)
%
%    Description:
%     WRITESODEVENTCSV(FILE,EVENTS) writes out a CSV (comma-separated
%     values) file like event .csv files produced with SOD (Standing Order
%     for Data).  FILE is the path/filename of the .csv file.  EVENTS is a
%     scalar struct (probably created with READSODEVENTCSV) with the fields
%     indicated in the Notes section below.  If FILE is an empty string a
%     graphical file creation menu is presented.
%
%     WRITESODEVENTCSV(FILE,EVENTS,OVERWRITE) quietly overwrites a pre-
%     existing CSV file without confirmation when OVERWRITE is set to
%     TRUE.  By default OVERWRITE is FALSE.  OVERWRITE is ignored in the
%     graphical file creation menu.
%
%    Notes:
%     - Converts latitude, longitude, depth and magnitude from numeric form
%       to strings (for function WRITECSV).
%     - Converts time from [yr mon cday hr min secs] to string.
%     - Fields of a standard SOD Event CSV file (in this order too):
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
%     % Using READSODEVENTCSV and WRITESODEVENTCSV together allows for
%     % adjusting SOD Event CSV files via matlab.  For instance to create
%     % a file with only events with magnitude > 5.5:
%     events=readsodeventcsv('my.csv')
%     events=ssidx(events,events.magnitude>5.5);
%     writesodeventcsv('my_gt_5p5.csv',events);
%
%    See also: READSODEVENTCSV, READCSV, WRITECSV, SSIDX, SSCAT

%     Version History:
%        Sep. 16, 2009 - initial version
%        Sep. 22, 2009 - overwrite confirmation skip option added
%        Jan. 26, 2010 - minor doc update, graphical file creation
%        Mar. 30, 2010 - check fields are available to modify, doc update
%        Feb. 11, 2011 - mass nargchk fix
%        Feb. 29, 2012 - update for scalar struct changes
%        Jan. 17, 2013 - fix example, add ssidx/sscat to see also
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 17, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% require certain fields are present
req={'time' 'latitude' 'longitude' 'depth' 'magnitude'};
fields=fieldnames(events);
if(~all(ismember(req,fields)))
    error('seizmo:readsodeventcsv:missingFields',...
        ['EVENTS struct must have the following fields:\n' ...
        sprintf('%s ',req{:})]);
end

% convert some fields from numbers to strings
f={'latitude' 'longitude' 'depth' 'magnitude'};
for i=1:numel(f)
    events.(f{i})=strtrim(cellstr(num2str(events.(f{i}))));
end

% convert time from numeric matrix to string
events.time=cellstr(reshape(sprintf(...
    '%04d-%02d-%02dT%02d:%02d:%06.3fZ',...
    fixtimes(events.time).'),24,[]).');

% write out
writecsv(file,events,varargin{:});

end
