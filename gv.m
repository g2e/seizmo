function [version,endian]=gv(filename)
%GV    Get version and byte-order of SAClab datafile
%
%    Description: Automatically determines version and byte-order of a
%     SAClab binary datafile.  Currently based on the header version field
%     validity - a 32bit signed integer occupying bytes 305 to 308.  If 
%     the datafile cannot be validated (usually occurs when the file is not
%     a SAClab datafile or cannot be opened) a warning is issued and the 
%     version is set to 0.
%
%    Usage:    [version,endian]=gv('filename')
%
%    Examples:
%
%    See also:  rh, wh, seishi, vvseis

% get valid versions
valid=vvseis();

% preset version/endian to invalid
version=0; endian='';

% open file for reading
fid=fopen(filename);

% check for invalid fid (for directories etc)
if(fid<0)
    warning('SAClab:gv:badFID','File not openable, %s',filename);
    return;
end

% seek to version field
try
    fseek(fid,304,'bof');
catch
    % seek failed
    fclose(fid);
    warning('SAClab:gv:fileTooShort','File too short, %s',filename);
    return;
end

% at end of file
if(feof(fid))
    % seeked to eof...
    fclose(fid);
    warning('SAClab:gv:fileTooShort','File too short, %s',filename);
    return;
end

% read in version as little-endian
endian='ieee-le';
try
    version=fread(fid,1,'int32',endian);
catch
    % read version failed - close file and warn
    fclose(fid);
    warning('SAClab:gv:readVerFail',...
        'Unable to read header version of file, %s',filename);
    version=0;
    return;
end

% check if valid
if(isempty(version))
    % read returned nothing...
    fclose(fid);
    warning('SAClab:gv:readVerFail',...
        'Unable to read header version of file, %s',filename);
    version=0;
    return;
elseif(~any(valid==version))
    % no good - seek back and read as big-endian
    fseek(fid,-4,'cof');
    endian='ieee-be';
    try
        version=fread(fid,1,'int32',endian);
    catch
        % read version failed - close file and warn
        fclose(fid);
        warning('SAClab:gv:readVerFail',...
            'Unable to read header version of file, %s',filename);
        version=0;
        return;
    end
    
    % check if valid
    if(~any(valid==version))
        % no good again - close file and warn
        fclose(fid);
        warning('SAClab:gv:versionUnknown',...
            'Unknown header version for file, %s',filename);
        version=0;
        return;
    end
end
fclose(fid);

end
