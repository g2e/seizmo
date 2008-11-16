function [filetype,version,endian]=getversion(filename)
%GETVERSION    Get filetype, version and byte-order of SEIZMO datafile
%
%    Description: [FILETYPE,VERSION,ENDIAN]=GETVERSION(FILENAME) determines
%     the filetype FILETYPE, version VERSION and byte-order ENDIAN of a
%     SEIZMO compatible file FILENAME.  If a datafile cannot be validated
%     (occurs when the file is not a SEIZMO datafile or cannot be opened or
%     read) a warning is given & FILETYPE, VERSION, ENDIAN are set empty.
%
%    Notes:
%     - Currently this is solely based on the header version field validity 
%       which is a 32bit signed integer occupying bytes 305 to 308.
%
%    Tested on: Matlab r2007b
%
%    Usage:    [filetype,version,endian]=getversion('filename')
%
%    Examples:
%     Figure out a file's version so that we can pull up the definition:
%      [filetype,version,endian]=getversion('myfile')
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
%        Oct. 17, 2008 - renamed from GV to GETVERSION, filetype outputj
%        Nov. 15, 2008 - update for new naming schema, better separation of
%                        methods for getting versions
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 15, 2008 at 20:45 GMT

% todo:

% check input
error(nargchk(1,1,nargin))
if(~ischar(filename))
    error('seizmo:getversion:badInput','FILENAME must be a string!');
end

% preset filetype/version/endian to invalid
filetype=[]; version=[]; endian=[];

% open file for reading
fid=fopen(filename);

% check for invalid fid (for directories etc)
if(fid<0)
    warning('seizmo:getversion:badFID',...
        'File not openable, %s !',filename);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% METHOD ONE -- SAC VERSION FIELD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    [filetype,version,endian]=getsacversion(fid);
    fclose(fid);
    return;
catch
    % go on to the next method below
end

% all methods failed - clean up and give some info
fclose(fid);
report=lasterror;
warning(report.identifier,report.message);

end


function [filetype,version,endian]=getsacversion(fid)
% gets version for sac files

try
    % seek to version field
    fseek(fid,304,'bof');
catch
    % seek failed
    error('seizmo:getversion:fileTooShort',...
        'File too short, %s !',filename);
end

% at end of file
if(feof(fid))
    % seeked to eof...
    error('seizmo:getversion:fileTooShort',...
        'File too short, %s !',filename);
end

try
    % read in version as little-endian
    ver=fread(fid,1,'int32','ieee-le');
catch
    % read version failed - close file and warn
    error('seizmo:getversion:readVerFail',...
        'Unable to read header version of file, %s !',filename);
end

% check for empty
if(isempty(ver))
    % read returned nothing...
    error('seizmo:getversion:readVerFail',...
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
    error('seizmo:getversion:readVerFail',...
        'Unable to read header version of file, %s !',filename);
end

% check for empty
if(isempty(ver))
    % read returned nothing...
    error('seizmo:getversion:readVerFail',...
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
error('seizmo:getversion:versionUnknown',...
    'Unknown/Unsupported version!');

end
