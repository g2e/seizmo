function [version,endian]=gv(filename)
%GV    Get version and byte-order of SAC file
%
%    Description: Automatically determines version and byte-order of a SAC
%     file based on the validity of the header version field - a 32bit 
%     signed integer occupying bytes 305 to 308).
%
%    Usage:    [version,endian]=gv('filename')
%
%    See also:  ch, lh, gh, rh, wh, rpdw, rdata, rsac, bsac, wsac, sachi 
%               doubleit, fixdelta, glgc, genum, genumdesc, sacsize

% get valid versions
valid=sachi();

% preset version/endian to invalid
version=0; endian='';

% open file for reading
fid=fopen(filename);

% check for invalid fid (for directories etc)
if(fid<0)
    warning('SAClab:badFID','File not openable, %s',filename);
    return;
end

% seek to version field
try
    fseek(fid,304,'bof');
catch
    % seek failed
    fclose(fid);
    warning('SAClab:fileTooShort','File too short, %s',filename);
    return;
end

% at end of file
if(feof(fid))
    % seeked to eof...
    fclose(fid);
    warning('SAClab:fileTooShort','File too short, %s',filename);
    return;
end

% read in version as little-endian
endian='ieee-le';
try
    version=fread(fid,1,'int32',endian);
catch
    % read version failed - close file and warn
    fclose(fid);
    warning('SAClab:readVerFail',...
        'Unable to read header version of file, %s',filename);
    version=0;
    return;
end

% check if valid
if(isempty(version))
    % read returned nothing...
    fclose(fid);
    warning('SAClab:readVerFail',...
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
        warning('SAClab:readVerFail',...
            'Unable to read header version of file, %s',filename);
        version=0;
        return;
    end
    
    % check if valid
    if(~any(valid==version))
        % no good again - close file and warn
        fclose(fid);
        warning('SAClab:versionUnknown',...
            'Unknown header version for file, %s',filename);
        version=0;
        return;
    end
end
fclose(fid);

end
