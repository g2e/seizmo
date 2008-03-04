function []=wh(data)
%WH    Write SAClab data headers to binary seismic datafiles
%
%    Description: Writes SAClab data headers to binary seismic datafiles.  
%     Basically this is for updating the headers of existing datafiles to 
%     match the SAClab data structure.  For files that do not already exist
%     the 'name' field is used for output file naming.  To write files to a
%     specific directory either add the path to the 'name' field or browse 
%     Matlab there using the 'cd' command.  Output file byte-order can be 
%     set using the 'endian' field.
%
%    Warning:  Using wh after modifying the version/byte-order in SAClab
%              could cause the file to become corrupted (eg. header has
%              one byte-order while the data has the other).  Please use
%              wseis if you want to modify the version/byte-order of a file
%              in order to avoid corruption issues.
%
%    Usage:    wh(SAClab_struct)
%
%    Examples:
%     Read in headers, modify them, and write out changes
%
%         data=rh(filelist)
%         data=ch(data,'kuser0','QCed')
%         wh(data)
%
%    See also:  wseis, rseis, bseis, rpdw, rdata, rh, seishi, gv

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'name','endian'))

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=seishi(vers(nver));
for i=1:nver-1
    h(i)=seishi(vers(i));
end

% loop over records
for i=1:length(data)
    % logical index of header info
    v=(data(i).version==vers);
    
    % open existing file for writing
    if(exist(data(i).name,'file'))
        % check for version/byte-order change
        [version,endian]=gv(data(i).name);
        if(version)
            if(version~=data(i).version)
                warning('SAClab:wh:versionMismatch',...
                    ['Version of existing file %s does ' ...
                    'not match output header version.\n' ...
                    'Accurate reading of data in existing file '...
                    'may not be possible now!'],data(i).name);
            end
            if(endian~=data(i).endian)
                warning('SAClab:wh:endianMismatch',...
                    ['Byte-order of existing file %s does '...
                    'not match output header byte-order.\n'...
                    'Accurate reading of data in existing file '...
                    'may not be possible now!'],data(i).name);
            end
            
            % open file for modification
            fid=fopen(data(i).name,'r+',data(i).endian);
        % file exists but is not SAClab datafile
        else
            warning('SAClab:wh:badFile',...
                ['Existing file %s is not a SAClab datafile.\n' ...
                'Attempting to overwrite file...'],data(i).name);
            % overwrite file
            fid=fopen(data(i).name,'w',data(i).endian);
        end
    else
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
            'Writing failed for file %s\n',data(i).name);
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
