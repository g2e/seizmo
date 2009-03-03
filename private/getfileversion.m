function [filetype,version,endian]=getfileversion(filename)
%GETFILEVERSION    Get filetype, version and byte-order of SEIZMO datafile
%
%    Description: [FILETYPE,VERSION,ENDIAN]=GETFILEVERSION(FILENAME) gets
%     the filetype FILETYPE, version VERSION, and byte-order ENDIAN of a
%     SEIZMO compatible file FILENAME.  If a datafile cannot be validated
%     (occurs when the file is not a SEIZMO datafile or cannot be opened or
%     read) a warning is given & FILETYPE, VERSION, ENDIAN are set empty.
%
%    Notes:
%     - Currently this is solely based on the sac header version field
%       validity (a 32bit signed integer occupying bytes 305 to 308).
%
%    Tested on: Matlab r2007b
%
%    Usage:    [filetype,version,endian]=getfileversion('filename')
%
%    Examples:
%     Figure out a file's version so that we can pull up the definition:
%      [filetype,version,endian]=getfileversion('myfile')
%      definition=seizmodef(version)
%
%    See also:  readheader, writeheader, seizmodef, validseizmo

%     Version History:
%        Jan. 27, 2008 - initial version
%        Feb. 23, 2008 - minor doc update
%        Feb. 28, 2008 - uses vvseis now, fix warnings
%        Mar.  4, 2008 - doc update, fix warnings
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - minor doc update, input checks
%        Oct. 17, 2008 - renamed from GV to getfileversion, filetype outputj
%        Nov. 15, 2008 - update for new naming schema, better separation of
%                        methods for getting versions
%        Mar.  3, 2009 - renamed from getfileversion to GETFILEVERSION,
%                        minor code cleaning
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  3, 2009 at 19:15 GMT

% todo:

% check input
error(nargchk(1,1,nargin))
if(~ischar(filename))
    error('seizmo:getfileversion:badInput','FILENAME must be a string!');
end

% preset filetype/version/endian to invalid
filetype=[]; version=[]; endian=[];

% open file for reading
fid=fopen(filename);

% check for invalid fid (for directories etc)
if(fid<0)
    warning('seizmo:getfileversion:badFID',...
        'File not openable, %s !',filename);
    return;
end

% try different methods, catching errors
try
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHOD ONE -- SAC VERSION FIELD %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [filetype,version,endian]=getsacversion(fid);
    fclose(fid);
    return;
catch
    % other methods would go here
end

% all methods failed - clean up and give some info
fclose(fid);
report=lasterror;
warning(report.identifier,report.message);

end


function [filetype,version,endian]=getsacversion(fid)
% gets version for sac files, throws error on non-sac files

try
    % seek to version field
    fseek(fid,304,'bof');
catch
    % seek failed
    error('seizmo:getfileversion:fileTooShort',...
        'File too short, %s !',filename);
end

% at end of file
if(feof(fid))
    % seeked to eof...
    error('seizmo:getfileversion:fileTooShort',...
        'File too short, %s !',filename);
end

try
    % read in version as little-endian
    ver=fread(fid,1,'int32','ieee-le');
catch
    % read version failed - close file and warn
    error('seizmo:getfileversion:readVerFail',...
        'Unable to read header version of file, %s !',filename);
end

% check for empty
if(isempty(ver))
    % read returned nothing...
    error('seizmo:getfileversion:readVerFail',...
        'Unable to read header version of file, %s !',filename);
end

% check if valid SAC or SEIZMO file
sac=validseizmo('SAC Binary');
seizmo=validseizmo('SEIZMO Binary');
if(any(sac==ver))
    filetype='SAC Binary';
    version=ver;
    endian='ieee-le';
    return;
elseif(any(seizmo==ver))
    filetype='SEIZMO Binary';
    version=ver;
    endian='ieee-le';
    return;
end

% no good - seek back
fseek(fid,-4,'cof');

try
    % read in version as big-endian
    ver=fread(fid,1,'int32','ieee-be');
catch
    % read version failed - close file and warn
    error('seizmo:getfileversion:readVerFail',...
        'Unable to read header version of file, %s !',filename);
end

% check for empty
if(isempty(ver))
    % read returned nothing...
    error('seizmo:getfileversion:readVerFail',...
        'Unable to read header version of file, %s !',filename);
end

% check if valid SAC or SEIZMO file
if(any(sac==ver))
    filetype='SAC Binary';
    version=ver;
    endian='ieee-be';
    return;
elseif(any(seizmo==ver))
    filetype='SEIZMO Binary';
    version=ver;
    endian='ieee-be';
    return;
end

% unknown version
error('seizmo:getfileversion:versionUnknown',...
    'Unknown/Unsupported version!');

end
