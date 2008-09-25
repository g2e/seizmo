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
%     recommended as it will likely corrupt your datafiles!
%
%    System requirements: Matlab 7
%
%    Input/Output requirements: Data structure must have 'name' and
%     'endian' fields.
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 22, 2008 at 07:15 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'name','endian'))

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=seisdef(vers(nver));
for i=1:nver-1
    h(i)=seisdef(vers(i));
end

% get this platform's native byte-order
[platform,maxint,nativeendian]=computer;
clear platform maxint
if(strcmpi(nativeendian,'L')); nativeendian='ieee-le';
else nativeendian='ieee-be'; end

% loop over records
for i=1:length(data)
    % logical index of header info
    v=(data(i).version==vers);
    
    % check for empty filename
    if(isempty(data(i).name))
        warning('SAClab:wh:namelessData',...
            ['Record %d has no associated filename!\n'...
            '==> Set name as rec%d.sac !'],i,i);
        data(i).name=['SAClab.' num2str(i) '.sac'];
    end
    
    % open existing file for writing
    if(exist(data(i).name,'file'))
        % get version/byte-order of datafile on disk
        [version,endian]=gv(data(i).name);
        
        % non-zero version ==> file is SAClab compatible
        if(version)
            % check if current data has no byte-order 
            % if so ==> set to match datafile on disk
            if(isempty(data(i).endian))
                warning('SAClab:wh:endianlessData',...
                    ['Record %d has no associated byteorder!\n'...
                    '==> Set byteorder as %s to match existing file!'],...
                    i,endian);
                data(i).endian=endian;
            end
            
            % check for version/byte-order change
            if(version~=data(i).version)
                warning('SAClab:wh:versionMismatch',...
                    ['Version of existing file %s does ' ...
                    'not match output header version!\n' ...
                    'Accurate reading of data in existing file '...
                    'may not be possible now!'],data(i).name);
            end
            if(endian~=data(i).endian)
                warning('SAClab:wh:endianMismatch',...
                    ['Byte-order of existing file %s does '...
                    'not match output header byte-order!\n'...
                    'Accurate reading of data in existing file '...
                    'may not be possible now!'],data(i).name);
            end
            
            % open file for modification
            fid=fopen(data(i).name,'r+',data(i).endian);
        % file exists but is not SAClab datafile
        else
            % check if current data has no byte-order 
            % if so ==> set to platforms native byteorder
            if(isempty(data(i).endian))
                warning('SAClab:wh:endianlessData',...
                    ['Record %d has no associated byteorder!\n'...
                    '==> Set byteorder as %s to match platform!'],...
                    i,nativeendian);
                data(i).endian=nativeendian;
            end
            
            % SAClab is gonna trash your file!
            warning('SAClab:wh:badFile',...
                ['Existing file %s is not a SAClab datafile.\n' ...
                'Attempting to overwrite file...'],data(i).name);
            
            % overwrite file
            fid=fopen(data(i).name,'w',data(i).endian);
        end
    % file doesn't exist ==> make new file
    else
        % check if current data has no byte-order 
        % if so ==> set to platforms native byteorder
        if(isempty(data(i).endian))
            warning('SAClab:wh:endianlessData',...
                ['Record %d has no associated byteorder!\n'...
                '==> Set byteorder as %s to match platform!'],...
                i,nativeendian);
            data(i).endian=nativeendian;
        end
        
        % new file
        fid=fopen(data(i).name,'w',data(i).endian);
    end
    
    % check fid
    if(fid<0)
        % unopenable file for writing (permissions/directory?)
        error('SAClab:wh:badFID',...
            ['File not openable for writing '...
            '(permissions/conflict with directory?):\n%s\n'],...
            data(i).name);
    end
    
    % fill header with dummy bytes (so we can seek around)
    fseek(fid,0,'bof');
    count=fwrite(fid,zeros(h(v).data.startbyte,1),'char');
    
    % verify write
    if(count<h(v).data.startbyte)
        % write failed
        fclose(fid);
        error('SAClab:wh:writeFailed',...
            'Writing failed for file %s !\n',data(i).name);
    end
    
    % write header
    n=h(v).types;
    for m=1:length(n)
        for k=1:length(h(v).(n{m}))
            fseek(fid,h(v).(n{m})(k).startbyte,'bof');
            fwrite(fid,data(i).head(h(v).(n{m})(k).minpos:h(v).(n{m})(k).maxpos),h(v).(n{m})(k).store);
        end
    end
    
    % close file
    fclose(fid);
end

end
