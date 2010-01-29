function []=writesodeventcsv(file,events,varargin)
%WRITESODEVENTCSV    Write out a SOD event .csv file from a structure
%
%    Usage:    writesodeventcsv(file,events)
%              writesodeventcsv(file,events,overwrite)
%
%    Description: WRITESODEVENTCSV(FILE,EVENTS) writes out a CSV (comma-
%     separated values) file similar to that of event .csv files produced
%     with SOD (Standing Order for Data).  FILE is the path/filename of the
%     .csv file.  EVENTS is a struct with the fields indicated in the Notes
%     section below.  If FILE is an empty string a graphical file creation
%     menu is presented.
%
%     WRITESODEVENTCSV(FILE,EVENTS,OVERWRITE) quietly overwrites pre-
%     existing CSV files without confirmation when OVERWRITE is set to
%     TRUE.  By default OVERWRITE is FALSE.  OVERWRITE is ignored in the
%     graphical file creation menu.
%
%    Notes:
%     - converts latitude, longitude, depth and magnitude from numeric form
%       to strings (for function WRITECSV)
%     - converts time from [yr mon cday hr min secs] to string
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
%     Using READSODEVENTCSV and WRITESODEVENTCSV together allows for
%     adjusting SOD Event CSV files via matlab.  For instance to create a
%     file with only events with magnitude > 5.5:
%      events=readsodeventcsv('my.csv')
%      events=events([events.magnitude]>5.5));
%      writesodeventcsv('my_gt_5p5.csv',events);
%
%    See also: READSODEVENTCSV, READCSV, WRITECSV

%     Version History:
%        Sep. 16, 2009 - initial version
%        Sep. 22, 2009 - overwrite confirmation skip option added
%        Jan. 26, 2010 - minor doc update, graphical file creation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2010 at 22:35 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
if(~isempty(msg)); error(msg); end;

% convert some fields from numbers to strings
f={'latitude' 'longitude' 'depth' 'magnitude'};
for i=1:numel(f)
    tmp=strtrim(cellstr(num2str([events.(f{i})].')));
    [events.(f{i})]=deal(tmp{:});
end

% convert time from numeric matrix to string
% - clean up times
% - need to get millisecond values
nev=numel(events);
tmp=[events.time];
tmp=reshape(tmp,6,nev);
tmp=fixtimes(tmp.').';
tmp(6,:)=round(1000*tmp(6,:));
tmp(7,:)=mod(tmp(6,:),1000);
tmp(6,:)=fix(tmp(6,:)/1000);
tmp=cellstr(reshape(...
    sprintf('%04d-%02d-%02dT%02d:%02d:%02d.%03dZ',tmp),24,nev).');
[events.time]=deal(tmp{:});

% write out
writecsv(file,events,varargin{:});

end
