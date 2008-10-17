function []=wh(data)
%WH    Write SAClab data header info to datafiles
%
%    Description: WH(DATA) writes SAClab data headers as datafiles on disk.  
%     Primarily this is for updating the headers of existing datafiles to 
%     match the SAClab DATA structure.  For files that do not already exist
%     the 'name' field is used as the output filename.  To write files to a
%     specific directory, either add the path to the 'name' field or browse 
%     there using the 'cd' command before using WH.  Output file byte-order
%     can be set using the 'endian' field.
%
%    Warning:  If you want to modify the version and/or byte-order of a
%     datafile with SAClab, use RSEIS and WSEIS to read/write the datafile
%     so that the data can also be adjusted to the new format.  Using RH 
%     and WH in conjunction with a version and/or byte-order change is NOT 
%     recommended as it only changes the header and will likely corrupt
%     your datafiles!
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage:    wh(data)
%
%    Examples:
%     Read in some datafile's headers, modify them, and write out changes:
%      wh(ch(rh('*.SAC'),'kuser0','QCed'))
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 15, 2008 at 13:30 GMT

% todo:
% - filetype support
% - clean up warnings

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data))

% check header
data=chkhdr(data);

% headers setup
[h,vi]=vinfo(data);

% loop over records
for i=1:numel(data)
    % construct fullname
    name=fullfile(data(i).dir,data(i).name);
    
    % open existing file for writing
    if(exist(name,'file'))
        % get version/byte-order of datafile on disk
        [filetype,fileversion,fileendian]=getversion(name);
        
        % non-zero version ==> file is SAClab compatible
        if(fileversion)
            % check for version/byte-order change
            if(fileversion~=data(i).version)
                warning('SAClab:wh:versionMismatch',...
                    ['Version of existing file %s does ' ...
                    'not match output header version!\n' ...
                    'Accurate reading of data in existing file '...
                    'may not be possible now!'],name);
            end
            if(fileendian~=data(i).endian)
                warning('SAClab:wh:endianMismatch',...
                    ['Byte-order of existing file %s does '...
                    'not match output header byte-order!\n'...
                    'Accurate reading of data in existing file '...
                    'may not be possible now!'],name);
            end
            
            % open file for modification
            fid=fopen(name,'r+',data(i).endian);
        % file exists but is not SAClab datafile
        else
            % SAClab is gonna trash your file!
            warning('SAClab:wh:badFile',...
                ['Existing file %s is not a SAClab datafile.\n' ...
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
        error('SAClab:wh:badFID',...
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
        error('SAClab:wh:writeFailed',...
            'Writing failed for file %s !\n',name);
    end
    
    % write header
    n=h(vi(i)).types;
    for m=1:length(n)
        for k=1:length(h(vi(i)).(n{m}))
            fseek(fid,h(vi(i)).(n{m})(k).startbyte,'bof');
            fwrite(fid,data(i).head(h(vi(i)).(n{m})(k).minpos:h(vi(i)).(n{m})(k).maxpos),h(vi(i)).(n{m})(k).store);
        end
    end
    
    % close file
    fclose(fid);
end

end
