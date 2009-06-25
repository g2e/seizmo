function [filetype,version,endian]=getfileversion(filename,verbose)
%GETFILEVERSION    Get filetype, version and byte-order of SEIZMO datafile
%
%    Usage:    [filetype,version,endian]=getfileversion('filename')
%              [filetype,version,endian]=getfileversion('filename',verbose)
%
%    Description: [FILETYPE,VERSION,ENDIAN]=GETFILEVERSION(FILENAME) gets
%     the filetype FILETYPE, version VERSION, and byte-order ENDIAN of a
%     SEIZMO compatible file FILENAME.  If a datafile cannot be validated
%     (occurs when the file is not a SEIZMO datafile or cannot be opened or
%     read) a warning is given & FILETYPE, VERSION, ENDIAN are set empty.
%
%     [FILETYPE,VERSION,ENDIAN]=GETFILEVERSION(FILENAME,VERBOSE) outputs
%     all error messages encountered during reading attempts as warnings.
%     This is useful for when a file is not reading in when it should be.
%
%    Notes:
%     - Currently this is solely based on the sac header version field
%       validity (a 32bit signed integer occupying bytes 305 to 308).
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
%        Oct. 17, 2008 - renamed from GV to GETVERSION, filetype output
%        Nov. 15, 2008 - update for new naming schema, better separation of
%                        methods for getting versions
%        Mar.  3, 2009 - renamed from GETVERSION to GETFILEVERSION,
%                        minor code cleaning
%        Apr.  7, 2009 - verbose option, output cleanup, better subfunction
%                        handling, total separation of type methods
%        Apr. 23, 2009 - fix for array of function handles (octave needs
%                        comma separated list), move usage up
%        June 12, 2009 - minor error msg change, add testing table
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2009 at 17:30 GMT

% todo:

% check input
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end

if(~ischar(filename))
    error('seizmo:getfileversion:badInput','FILENAME must be a string!');
end
if(nargin==1 || isempty(verbose))
    verbose=false;
end
if(~islogical(verbose))
    error('seizmo:getfileversion:badInput','VERBOSE must be logical!');
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

% methods
funcs={@getsacbinaryversion, @getseizmobinaryversion};

% try different methods, catching errors
for i=1:numel(funcs)
    try
        [filetype,version,endian]=funcs{i}(fid);
        fclose(fid);
        return;
    catch
        % report error if verbose
        if(verbose)
            report=lasterror;
            warning(report.identifier,report.message);
        end
    end
end

% all methods failed - clean up and give some info
fclose(fid);
warning('seizmo:getfileversion:typeUnknown',...
    'File: %s\nUnknown Filetype! Skipping!',filename)

end


function [filetype,version,endian]=getsacbinaryversion(fid)
% gets version for sac binary files, throws error on non-sac files

try
    % seek to version field
    fseek(fid,304,'bof');
catch
    % seek failed
    error('seizmo:getsacbinaryversion:fileTooShort',...
        'SAC Binary File too short, %s !',filename);
end

% at end of file
if(feof(fid))
    % seeked to eof...
    error('seizmo:getsacbinaryversion:fileTooShort',...
        'SAC Binary File too short, %s !',filename);
end

try
    % read in version as little-endian
    ver=fread(fid,1,'int32','ieee-le');
catch
    % read version failed - close file and warn
    error('seizmo:getsacbinaryversion:readVerFail',...
        'Unable to read SAC Binary header version of file, %s !',filename);
end

% check for empty
if(isempty(ver))
    % read returned nothing...
    error('seizmo:getsacbinaryversion:readVerFail',...
        'Unable to read SAC Binary header version of file, %s !',filename);
end

% check if valid SAC file
sac=validseizmo('SAC Binary');
if(any(sac==ver))
    filetype='SAC Binary';
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
    error('seizmo:getsacbinaryversion:readVerFail',...
        'Unable to read SAC Binary header version of file, %s !',filename);
end

% check for empty
if(isempty(ver))
    % read returned nothing...
    error('seizmo:getsacbinaryversion:readVerFail',...
        'Unable to read SAC Binary header version of file, %s !',filename);
end

% check if valid SAC or SEIZMO file
if(any(sac==ver))
    filetype='SAC Binary';
    version=ver;
    endian='ieee-be';
    return;
end

% unknown version
error('seizmo:getsacbinaryversion:versionUnknown',...
    'File: %s\nUnknown/Unsupported SAC Binary Version!',filename);

end


function [filetype,version,endian]=getseizmobinaryversion(fid)
% gets version for sac binary files, throws error on non-sac files

try
    % seek to version field
    fseek(fid,304,'bof');
catch
    % seek failed
    error('seizmo:getseizmobinaryversion:fileTooShort',...
        'SEIZMO Binary File too short, %s !',filename);
end

% at end of file
if(feof(fid))
    % seeked to eof...
    error('seizmo:getseizmobinaryversion:fileTooShort',...
        'SEIZMO Binary File too short, %s !',filename);
end

try
    % read in version as little-endian
    ver=fread(fid,1,'int32','ieee-le');
catch
    % read version failed - close file and warn
    error('seizmo:getseizmobinaryversion:readVerFail',...
        'Unable to read SEIZMO Binary header version of file, %s !',...
        filename);
end

% check for empty
if(isempty(ver))
    % read returned nothing...
    error('seizmo:getseizmobinaryversion:readVerFail',...
        'Unable to read SEIZMO Binary header version of file, %s !',...
        filename);
end

% check if valid SEIZMO file
seizmo=validseizmo('SEIZMO Binary');
if(any(seizmo==ver))
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
    error('seizmo:getseizmobinaryversion:readVerFail',...
        'Unable to read SEIZMO Binary header version of file, %s !',...
        filename);
end

% check for empty
if(isempty(ver))
    % read returned nothing...
    error('seizmo:getseizmobinaryversion:readVerFail',...
        'Unable to read SEIZMO Binary header version of file, %s !',...
        filename);
end

% check if valid SEIZMO file
if(any(seizmo==ver))
    filetype='SEIZMO Binary';
    version=ver;
    endian='ieee-be';
    return;
end

% unknown version
error('seizmo:getseizmobinaryversion:versionUnknown',...
    'File: %s\nUnknown/Unsupported SEIZMO Binary Version!',filename);

end
