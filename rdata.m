function [data,destroy]=rdata(data,trim)
%RH    Read binary seismic datafile data
%
%    Description: Reads all binary seismic datafile data into a seislab
%     structure using the info from the input seislab structure. Optional
%     logical parameter trim sets rdata to delete entries in data that had
%     errors (default is true). Optional output destroy is the logical
%     matrix indicating which records to remove.
%
%     Fields for timeseries files:
%      data.x(:,1) - amplitudes
%      data.t - times (if uneven spacing)
%
%     Fields for spectral amp/phase files:
%      data.x(:,1) - spectral amplitudes
%      data.x(:,2) - spectral phase
%
%     Fields for spectral real/imag files:
%      data.x(:,1) - spectral real
%      data.x(:,2) - spectral imaginary
%
%     Fields for general xy files:
%      data.x(:,1) - dependent component
%      data.t - independent component (if uneven spacing)
%
%     Fields for xyz grid files:
%      data.x(:,1) - matrix data (nodes evenly spaced; advances l2r,b2t)
%     
%    Notes:
%     - multi-component files will replicate the number of columns by the
%       number of components.  So a three component spectral file will have
%       six columns total.
%
%    Usage:    [data,destroy]=rdata(data,trim)
%
%    See also: rh, rpdw, rseis, wseis, bseis, seishi, gv, seissize

% check number of inputs
error(nargchk(1,2,nargin))

% check data structure
error(seischk(data,'name','endian'))

% default trim
if(nargin==1 || isempty(trim) || ...
        ~islogical(trim) || ~isscalar(trim))
    trim=true;
end

% number of records
nrecs=length(data);

% estimated filesize from header
est_bytes=seissize(data);

% header info
leven=glgc(data,'leven');
iftype=genumdesc(data,'iftype');
warning('off','seislab:gh:fieldInvalid')
[npts,ncmp]=gh(data,'npts','ncmp');
warning('on','seislab:gh:fieldInvalid')

% clean up and check ncmp
ncmp(isnan(ncmp))=1;
if(any(ncmp<1 | fix(ncmp)~=ncmp))
    error('seislab:rdata:badNumCmp',...
        'field ncmp must be a positive integer')
end

% check leven
t=strcmp(leven,'true');
f=strcmp(leven,'false');
if(~all(t | f))
    error('sieslab:rdata:levenBad',...
        'logical field leven needs to be set'); 
end

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=seishi(vers(nver));
for i=1:nver-1
    h(i)=seishi(vers(i));
end

% read loop
destroy=false(nrecs,1);
for i=1:nrecs
    % logical index of header info
    v=(data(i).version==vers);
    
    % open file for reading
    fid=fopen(data(i).name,'r',data(i).endian);
    
    % fid check
    if(fid<0)
        % non-existent file or directory
        warning('seislab:rdata:badFID',...
            'File not openable, %s',data(i).name);
        destroy(i)=1;
        continue;
    end
    
    % file size
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    
    % byte size check
    if(bytes>est_bytes(i))
        % size big enough but inconsistent - read anyways
        warning('seislab:rdata:badFileSize',...
            ['Filesize does not match header info (gt)\n'...
            'File: %s'],data(i).name);
    elseif(bytes<est_bytes(i))
        % size to small - skip
        fclose(fid);
        warning('seislab:rdata:badFileSize',...
            ['Filesize does not match header info (lt)\n'...
            'File: %s'],data(i).name);
        destroy(i)=1;
        continue;
    end
    
    % preallocate data record with NaNs, deallocate timing
    data(i).x=nan(npts(i),ncmp(i),h(v).data.store);
    data(i).t=[];
    
    % skip if npts==0
    if(npts(i)<1); fclose(fid); continue; end
    
    % act by file type (any new filetypes will have to be added here)
    fseek(fid,h(v).data.startbyte,'bof');
    if(strcmp(iftype(i),'Time Series File'))
        % time series file - amplitude and time
        for k=1:ncmp(i)
            data(i).x(:,k)=fread(fid,npts(i),['*' h(v).data.store]);
        end
        
        % timing of amp data if uneven
        if(strcmp(leven(i),'false'))
            data(i).t(:,1)=fread(fid,npts(i),['*' h(v).data.store]);
        end
    elseif(strcmp(iftype(i),'Spectral File-Real/Imag'))
        % preallocate data record with NaNs
        data(i).x=nan(npts(i),2*ncmp(i),h(v).data.store);
        
        % spectral file - real and imaginary
        for k=1:ncmp(i)
            data(i).x(:,2*k-1)=fread(fid,npts(i),['*' h(v).data.store]);
            data(i).x(:,2*k)=fread(fid,npts(i),['*' h(v).data.store]);
        end
    elseif(strcmp(iftype(i),'Spectral File-Ampl/Phase'))
        % preallocate data record with NaNs
        data(i).x=nan(npts(i),2*ncmp(i),h(v).data.store);
        
        % spectral file - amplitude and phase
        for k=1:ncmp(i)
            data(i).x(:,2*k-1)=fread(fid,npts(i),['*' h(v).data.store]);
            data(i).x(:,2*k)=fread(fid,npts(i),['*' h(v).data.store]);
        end
    elseif(strcmp(iftype(i),'General X vs Y file'))
        % general x vs y data (x is 'dependent')
        for k=1:ncmp(i)
            data(i).x(:,k)=fread(fid,npts(i),['*' h(v).data.store]);
        end
        
        % independent data (if uneven)
        if(strcmp(leven(i),'false'))
            data(i).t(:,1)=fread(fid,npts(i),['*' h(v).data.store]);
        end
    elseif(strcmp(iftype(i),'General XYZ (3-D) file'))
        % general xyz (3D) grid - nodes are evenly spaced
        for k=1:ncmp(i)
            data(i).x(:,k)=fread(fid,npts(i),['*' h(v).data.store]);
        end
    else
        % unknown filetype
        fclose(fid)
        warning('seislab:rdata:iftypeBad',...
            'File: %s\nBad filetype: %s',data(i).name,iftype(i));
        destroy(i)=1;
        continue;
    end
    
    % closing file
    fclose(fid);
end

% remove unread entries
if(trim); data(destroy)=[]; end

end
