function []=wh(data)
%WH    Write SEIZMO data header info to datafiles
%
%    Description: WH(DATA) writes SEIZMO data headers as datafiles on disk.  
%     Primarily this is for updating the headers of existing datafiles to 
%     match the SEIZMO DATA structure.  Struct fields 'location' and 'name'
%     are used to create or write to a file on disk.  Byte-order is set
%     using the 'endian' field.
%
%    Warning:  
%     If you want to modify the filetype, version or byte-order of
%     datafiles with SEIZMO, use RSEIS and WSEIS to read/write the entire
%     datafile so that the data can also be adjusted to the new format.
%     Using RH and WH in conjunction with a filetype, version or byte-order
%     change is NOT recommended as it only changes the header and will
%     likely corrupt your datafiles!
%
%    Tested on: Matlab r2007b
%
%    Header changes: NONE
%
%    Usage:    wh(data)
%
%    Examples:
%     Read in some datafile's headers, modify them, and write out changes:
%      wh(ch(rh('A.SAC'),'kuser0','QCed'))
%
%    See also:  wseis, rseis, bseis, rpdw, rdata, rh, seisdef, gv

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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 27, 2008 at 04:00 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% headers setup (checks struct too)
[h,vi]=vinfo(data);

% loop over records
for i=1:numel(data)
    % construct fullname
    name=fullfile(data(i).location,data(i).name);
    
    % open existing file for writing
    if(exist(name,'file'))
        % get version/byte-order of datafile on disk
        [filetype,fileversion,fileendian]=getversion(name);
        
        % non-zero version ==> file is SEIZMO compatible
        if(~isempty(filetype))
            % check for filetype/version/byte-order change
            if(filetype~=data(i).filetype)
                warning('seizmo:wh:filetypeMismatch',...
                    ['Filetype of existing file %s does NOT '...
                    'match the output filetype!\n'...
                    'Data corruption is likely to occur!'],name);
            end
            if(fileversion~=data(i).version)
                warning('seizmo:wh:versionMismatch',...
                    ['Version of existing file %s does ' ...
                    'NOT match output filetype version!\n' ...
                    'Data corruption is likely to occur!'],name);
            end
            if(fileendian~=data(i).endian)
                warning('seizmo:wh:endianMismatch',...
                    ['Byte-order of existing file %s does '...
                    'NOT match output header byte-order!\n'...
                    'Data corruption is likely to occur!'],name);
            end
            
            % open file for modification
            fid=fopen(name,'r+',data(i).endian);
        % file exists but is not SEIZMO datafile
        else
            % SEIZMO is gonna trash your file!
            warning('seizmo:wh:badFile',...
                ['Existing file %s is not a SEIZMO datafile.\n' ...
                'Attempting to overwrite file...'],name);
            
            % overwrite file
            fid=fopen(name,'w',data(i).endian);
        end
    % file doesn't exist ==> make new file
    else
        % new file
        fid=fopen(name,'w',data(i).endian);
    end
    
    % check fid
    if(fid<0)
        % unopenable file for writing (permissions/directory?)
        error('seizmo:wh:badFID',...
            ['File not openable for writing!\n'...
            '(Permissions problem / conflict with directory?)\n'...
            'File: %s\n'],name);
    end
    
    % fill header with dummy bytes (so we can seek around)
    fseek(fid,0,'bof');
    count=fwrite(fid,zeros(h(vi(i)).data.startbyte,1),'char');
    
    % verify write
    if(count<h(vi(i)).data.startbyte)
        % write failed
        fclose(fid);
        error('seizmo:wh:writeFailed',...
            'Writing failed for file %s !\n',name);
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
