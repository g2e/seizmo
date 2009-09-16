function [events]=readsodeventcsv(file)
%READSODEVENTCSV    Read in a SOD event .csv formatted file as a structure
%
%    Usage:    events=readsodeventcsv(file)
%
%    Description: EVENTS=READSODEVENTCSV(FILE) reads in a SOD (Standing
%     Order for Data) generated event CSV file as a structure array.  The
%     CSV (comma-separated values) file is expected to have 1 header line
%     of comma-separated field names.  Each line in FILE is placed as a
%     separated index in EVENTS, such that all values from line 4 (counting
%     the header line) can be accessed in EVENTS(3).
%
%    Notes:
%     - converts latitude, longitude, depth and magnitude to numeric form
%     - converts time from string to [yr mon cday hr min secs]
%     - SOD (Standing Order for Data) was written by Philip Crotwell
%     - fields of a standard SOD Event CSV file (in this order too):
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
%    Examples:
%     READSODEVENTCSV uses readcsv but also converts several fields into
%     a more useful form.  To show the differences:
%      ev1=readcsv(file);
%      ev2=readsodeventcsv(file);
%      ev1(1)
%      ev2(1)
%
%    See also: writesodeventcsv, readcsv, writecsv

%     Version History:
%        Sep. 16, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 16, 2009 at 22:30 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end;

% send on to readcsv
events=readcsv(file);

% convert some fields to numbers
f={'latitude' 'longitude' 'magnitude'};
for i=1:numel(f)
    tmp=num2cell(str2double({events.(f{i})}));
    [events.(f{i})]=deal(tmp{:});
end

% clean up depths
tmp=num2cell(round(100*str2double({events.depth}))/100);
[events.depth]=deal(tmp{:});

% convert time to something useful
% - note that we strip off the Z first b/c there is
%   a bug in dtstr2dtvecmx (datevec) that occasionally
%   returns an incorrect time when the Z is included
tmp=mat2cell(datevec(strnlen({events.time}.',23),...
    'yyyy-mm-ddTHH:MM:SS.FFF'),ones(numel(events),1));
[events.time]=deal(tmp{:});

end
