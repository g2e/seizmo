function [filetype,version,endian]=getversion(filename)
%GETVERSION    Get filetype, version and byte-order of SAClab datafile
%
%    Description: [FILETYPE,VERSION,ENDIAN]=GETVERSION(FILENAME) determines
%     the filetype FILETYPE, version VERSION and byte-order ENDIAN of a
%     SAClab compatible file FILENAME.  If a datafile cannot be validated
%     (occurs when the file is not a SAClab datafile or cannot be opened or
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
%      definition=seisdef(version)
%
%    See also:  rh, wh, seisdef, vvseis

%     Version History:
%        Jan. 27, 2008 - initial version
%        Feb. 23, 2008 - minor doc update
%        Feb. 28, 2008 - uses vvseis now, fix warnings
%        Mar.  4, 2008 - doc update, fix warnings
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - minor doc update, input checks
%        Oct. 17, 2008 - renamed from GV to GETVERSION, filetype output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 17, 2008 at 00:55 GMT

% todo:

% check input
error(nargchk(1,1,nargin))
if(~ischar(filename))
    error('SAClab:gv:badInput','FILENAME must be a string!');
end

% preset filetype/version/endian to invalid
filetype=[]; version=[]; endian=[];

% open file for reading
fid=fopen(filename);

% get valid versions for SAClab Binary
valid=vvseis('SAClab Binary');

% check for invalid fid (for directories etc)
if(fid<0)
    warning('SAClab:gv:badFID','File not openable, %s !',filename);
    return;
end

try
    % seek to version field
    fseek(fid,304,'bof');
catch
    % seek failed
    fclose(fid);
    warning('SAClab:gv:fileTooShort','File too short, %s !',filename);
    return;
end

% at end of file
if(feof(fid))
    % seeked to eof...
    fclose(fid);
    warning('SAClab:gv:fileTooShort','File too short, %s !',filename);
    return;
end

try
    % read in version as little-endian
    ver=fread(fid,1,'int32','ieee-le');
catch
    % read version failed - close file and warn
    fclose(fid);
    warning('SAClab:gv:readVerFail',...
        'Unable to read header version of file, %s !',filename);
    return;
end

% check if valid
if(isempty(ver))
    % read returned nothing...
    fclose(fid);
    warning('SAClab:gv:readVerFail',...
        'Unable to read header version of file, %s !',filename);
    return;
elseif(~any(valid==ver))
    % no good - seek back
    fseek(fid,-4,'cof');
    
    try
        % read in version as big-endian
        ver=fread(fid,1,'int32','ieee-be');
    catch
        % read version failed - close file and warn
        fclose(fid);
        warning('SAClab:gv:readVerFail',...
            'Unable to read header version of file, %s !',filename);
        return;
    end
    
    % check if valid
    if(~any(valid==ver))
        % no good again - close file and warn
        fclose(fid);
        warning('SAClab:gv:versionUnknown',...
            'Unknown header version for file, %s !',filename);
        return;
    else
        filetype='SAClab Binary';
        version=ver;
        endian='ieee-be';
    end
else
    filetype='SAClab Binary';
    version=ver;
    endian='ieee-le';
end
fclose(fid);

end
