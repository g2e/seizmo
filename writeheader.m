function []=writeheader(data)
%WRITEHEADER    Write SEIZMO data header info to datafiles
%
%    Description: WRITEHEADER(DATA) writes SEIZMO data headers as datafiles
%     on disk.  Primarily this is for updating the headers of existing
%     datafiles to match the SEIZMO DATA structure.
%
%    Warning:  
%     If you want to modify the filetype, version or byte-order of
%     datafiles with SEIZMO, use READSEIZMO and WRITESEIZMO to read/write
%     the entire datafile so that the data can also be adjusted to the new
%     format.  Using READHEADER and WRITEHEADER in conjunction with a
%     filetype, version or byte-order change is NOT recommended as it only
%     changes the header and will likely corrupt your datafiles!
%
%    Tested on: Matlab r2007b
%
%    Header changes: NONE
%
%    Usage:    writeheader(data)
%
%    Examples:
%     Read in some datafile's headers, modify them, and write out changes:
%      writeheader(changeheader(readheader('A.SAC'),'kuser0','QCed'))
%
%    See also:  writeseizmo, readseizmo, bseizmo, readdatawindow, readdata,
%               readheader, seizmodef, getfileversion, changeheader,
%               getheader, listheader, getlgc, getenumid, getenumdesc

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 11, 2008 - new SACHP support, better checks
%        Mar.  4, 2008 - code cleaning, more checks, doc update
%        June 12, 2008 - doc update
%        Sep. 15, 2008 - history fix, doc update
%        Sep. 22, 2008 - blank name and endian handling
%        Sep. 26, 2008 - VINFO & NATIVEENDIAN added
%        Oct. 15, 2008 - new SEISCHK support (blank name and endian not
%                        allowed) => NATIVEENDIAN dropped
%        Oct. 17, 2008 - supports new struct layout, added CHKHDR support
%        Oct. 27, 2008 - update for struct changes, remove CHKHDR so user
%                        can write whatever they want to disk
%        Nov. 17, 2008 - update for new name schema (now WRITEHEADER)
%        Dec. 13, 2008 - added mkdir call to make sure path exists
%        Mar.  3, 2009 - update for GETFILEVERSION
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  3, 2009 at 19:20 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% headers setup (checks struct too)
[h,vi]=versioninfo(data);

% loop over records
for i=1:numel(data)
    % construct fullname
    name=fullfile(data(i).path,data(i).name);
    
    % make sure directory exists
    [ok,msg,msgid]=mkdir(data(i).path);
    if(~ok)
        warning(msgid,msg);
        error('seizmo:writeheader:pathBad',...
            ['Record: %d, File: %s\n' ...
            'Cannot write record to path!'],i,name);
    end
    
    % open existing file for writing
    if(exist(name,'file'))
        % get version/byte-order of datafile on disk
        [filetype,fileversion,fileendian]=getfileversion(name);
        
        % non-zero version ==> file is SEIZMO compatible
        if(~isempty(filetype))
            % check for filetype/version/byte-order change
            if(filetype~=data(i).filetype)
                warning('seizmo:writeheader:filetypeMismatch',...
                    ['Filetype of existing file does '...
                    'NOT match the output filetype!\n'...
                    'Record: %d, File: %s\n' ...
                    'Data corruption is likely to occur!'],i,name);
            end
            if(fileversion~=data(i).version)
                warning('seizmo:writeheader:versionMismatch',...
                    ['Version of existing file does ' ...
                    'NOT match output filetype version!\n' ...
                    'Record: %d, File: %s\n' ...
                    'Data corruption is likely to occur!'],i,name);
            end
            if(fileendian~=data(i).byteorder)
                warning('seizmo:writeheader:endianMismatch',...
                    ['Byte-order of existing file does '...
                    'NOT match output header byte-order!\n'...
                    'Record: %d, File: %s\n' ...
                    'Data corruption is likely to occur!'],i,name);
            end
            
            % open file for modification
            fid=fopen(name,'r+',data(i).byteorder);
        % file exists but is not SEIZMO datafile
        else
            % SEIZMO is gonna trash your file!
            warning('seizmo:writeheader:badFile',...
                ['Existing file is not a SEIZMO datafile.\n' ...
                'Record: %d, File: %s\n' ...
                'Attempting to overwrite file...'],i,name);
            
            % overwrite file
            fid=fopen(name,'w',data(i).byteorder);
        end
    % file doesn't exist ==> make new file
    else
        % new file
        fid=fopen(name,'w',data(i).byteorder);
    end
    
    % check fid
    if(fid<0)
        % unopenable file for writing (permissions/directory?)
        error('seizmo:writeheader:badFID',...
            ['File not openable for writing!\n'...
            '(Permissions problem / conflict with directory?)\n'...
            'Record: %d, File: %s\n'],i,name);
    end
    
    % fill header with dummy bytes (so we can seek around)
    fseek(fid,0,'bof');
    count=fwrite(fid,zeros(h(vi(i)).data.startbyte,1),'char');
    
    % verify write
    if(count<h(vi(i)).data.startbyte)
        % write failed
        fclose(fid);
        error('seizmo:writeheader:writeFailed',...
            'Writing failed!\nRecord: %d, File: %s',i,name);
    end
    
    % write header
    n=h(vi(i)).types;
    for m=1:length(n)
        for k=1:length(h(vi(i)).(n{m}))
            fseek(fid,h(vi(i)).(n{m})(k).startbyte,'bof');
            fwrite(fid,data(i).head(h(vi(i)).(n{m})(k).minpos:...
                h(vi(i)).(n{m})(k).maxpos),h(vi(i)).(n{m})(k).store);
        end
    end
    
    % close file
    fclose(fid);
end

end
