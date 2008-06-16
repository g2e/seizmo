function []=wseis(data)
%WSEIS    Write SAClab data as binary datafiles
%
%    Description: Writes SAClab data as binary datafiles.  Uses the 'name'
%     field for naming output files.  To write files to a specific 
%     directory either add the path to the name field or browse Matlab 
%     there using the 'cd' command.  Output file's endianness can be set 
%     using the 'endian' field.
%
%    Usage:    wseis(SAClab_struct)
%
%    Examples:
%     Read in some files, clean them up and write over:
%      wseis(taper(rtrend(rseis('*'))))
%
%     To write out all records in big endian:
%      data.endian=deal('ieee-be');
%      wseis(data)
%
%     Alter the first records path/filename and write only to it:
%      data(1).name='../../myfile';
%      wseis(data(1))
%
%    See also:  rseis, bseis, seisdef, gv, rpdw, rdata, rh, wh

%     Version History:
%        ????????????? - Initial Version
%        June 12, 2008 - Added examples and updated documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2008 at 05:15 GMT

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'name','endian'))

% header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
iftype=genumdesc(data,'iftype');
warning('off','SAClab:gh:fieldInvalid')
[npts,ncmp]=gh(data,'npts','ncmp');
warning('on','SAClab:gh:fieldInvalid')

% clean up and check ncmp
ncmp(isnan(ncmp))=1;
if(any(ncmp<1 | fix(ncmp)~=ncmp))
    error('SAClab:rdata:badNumCmp',...
        'field ncmp must be a positive integer')
end

% estimated filesize from header
est_bytes=seissize(data);

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=seisdef(vers(nver));
for i=1:nver-1
    h(i)=seisdef(vers(i));
end

% loop over records
for i=1:length(data)
    % logical index of header info
    v=(data(i).version==vers);
    
    % open file for writing
    fid=fopen(data(i).name,'w',data(i).endian);
    
    % check fid
    if(fid<0)
        % unopenable file for writing (permissions/directory?)
        error('SAClab:wseis:badFID',...
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
        error('SAClab:wseis:writeFailed',...
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
    
    % skip if npts==0 (dataless)
    if(npts(i)<1); fclose(fid); continue; end
    
    % act by file type (new filetypes will have to be added here)
    fseek(fid,h(v).data.startbyte,'bof');
    if(strcmp(iftype(i),'Time Series File'))
        % time series file - amplitude and time
        for k=1:ncmp(i)
            fwrite(fid,data(i).x(:,k),h(v).data.store);
        end
        
        % timing of amp data if uneven
        if(strcmp(leven(i),'false'))
            fwrite(fid,data(i).t(:,1),h(v).data.store);
        end
    elseif(strcmp(iftype(i),'Spectral File-Real/Imag'))
        % spectral file - real and imaginary
        for k=1:ncmp(i)
            fwrite(fid,data(i).x(:,2*k-1),h(v).data.store);
            fwrite(fid,data(i).x(:,2*k),h(v).data.store);
        end
    elseif(strcmp(iftype(i),'Spectral File-Ampl/Phase'))
        % spectral file - amplitude and phase
        for k=1:ncmp(i)
            fwrite(fid,data(i).x(:,2*k-1),h(v).data.store);
            fwrite(fid,data(i).x(:,2*k),h(v).data.store);
        end
    elseif(strcmp(iftype(i),'General X vs Y file'))
        % general x vs y data (x is 'dependent')
        for k=1:ncmp(i)
            fwrite(fid,data(i).x(:,k),h(v).data.store);
        end
        
        % independent data (if uneven)
        if(strcmp(leven(i),'false'))
            fwrite(fid,data(i).t(:,1),h(v).data.store);
        end
    elseif(strcmp(iftype(i),'General XYZ (3-D) file'))
        % general xyz (3D) grid - nodes are evenly spaced
        for k=1:ncmp(i)
            fwrite(fid,data(i).x(:,k),h(v).data.store);
        end
    else
        % unknown filetype
        fclose(fid)
        error('SAClab:wseis:iftypeBad',...
            'File: %s\nBad filetype: %s',data(i).name,iftype(i));
    end
    
    % verify written file's size
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    if(bytes~=est_bytes(i))
        % write failed/incomplete
        fclose(fid);
        error('SAClab:wseis:badFileSize',...
            ['Output file''s written size does not match expected size\n'...
            'File: %s'],data(i).name);
    end
    
    % close file
    fclose(fid);
end

end
