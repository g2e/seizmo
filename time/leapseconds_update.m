function []=leapseconds_update(url)
%LEAPSECONDS_UPDATE    Updates the leapseconds file
%
%    Usage:    leapseconds_update
%              leapseconds_update(url)
%
%    Description:
%     LEAPSECONDS_UPDATE saves leap second information from a trusted
%     online source (http://maia.usno.navy.mil/ser7/leapsec.dat) to a file
%     named 'leapsec.dat' in the same directory as this function.
%
%     LEAPSECONDS_UPDATE(URL) allows using a different location to download
%     the leapsecond info.  Please note that the file must have the same
%     formatting as the info at the default url or it will break the
%     function LEAPSECONDS which all UTC timing depends on.
%
%    Notes:
%     - Requires ability to write to SEIZMO's installation directory.
%     - Since leapseconds were introduced in 1972, UTC times before that
%       are not properly accounted for with leap seconds to maintain timing
%       near UT1 (GMT).  There was actually a different method implemented
%       but that is not a matter for this function.  The data for that is
%       found here: http://maia.usno.navy.mil/ser7/tai-utc.dat
%
%    Examples:
%     % You should update the leap second file every time a new leap second
%     % is added to UTC (update before the leap second actually occurs):
%     leapseconds_update
%
%    See also: LEAPSECONDS, UTCOFFSET, LOD, FIXTIMES, TIMEDIFF, UTC2TAI,
%              TAI2UTC, ISLEAPYEAR

%     Version History:
%        Nov.  1, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% default url
if(~nargin || isempty(url))
    url='http://maia.usno.navy.mil/ser7/leapsec.dat';
end

% check url
if(~ischar(url) || ~isvector(url))
    error('seizmo:leapseconds_update:badURL',...
        'URL must be a string!');
end

% read info
[txt,ok]=urlread(url);

% check existence & non-empty
if(~ok || isempty(txt))
    error('seizmo:leapseconds_update:badURL',...
        'No leap second info found at: %s',url);
end

% file check
words=getwords(txt);
if(fix(numel(txt)/81)~=numel(txt)/81 ...
        || fix(numel(words)/15)~=numel(words)/15 ...
        || any(str2double(words(1:15:end))<1972) ...
        || any(fix(str2double(words(1:15:end)))...
            ~=str2double(words(1:15:end))) ...
        || any(~strcmp(unique(words(2:15:end)),{'JAN' 'JUL'})) ...
        || ~isscalar(unique(words(3:15:end))) ...
        || ~strcmp(unique(words(3:15:end)),'1') ...
        || ~strcmp(unique(words(4:15:end)),'=JD') ...
        || ~strcmp(unique(words(6:15:end)),'TAI-UTC='))
    error('seizmo:leapseconds_update:badURL',...
        'Malformed leap second info found at: %s',url);
end

% open for writing
path=fileparts(mfilename('fullpath'));
fid=fopen([path filesep 'leapsec.dat'],'wt');
if(fid<0)
    error('seizmo:leapseconds_update:cannotOpenFile',...
        'File: %s\nNot Openable!',[path filesep 'leapsec.dat']);
end

% write & close
fwrite(fid,txt,'char');
fclose(fid);

end
