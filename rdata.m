function [data]=rdata(data)
%RH    Read SAC binary file data
%
%    Description: Reads all SAC (seismic analysis code) binary file data
%     into a SAClab structure using the input SAClab structure header info.
%
%     Fields for SAC timeseries files:
%      data.x(:,1) - amplitudes
%      data.t - times (if uneven spacing)
%
%     Fields for SAC spectral amp/phase files:
%      data.x(:,1) - spectral amplitudes
%      data.x(:,2) - spectral phase
%
%     Fields for SAC spectral real/imag files:
%      data.x(:,1) - spectral real
%      data.x(:,2) - spectral imaginary
%
%     Fields for SAC general xy files:
%      data.x(:,1) - dependent component
%      data.t - independent component (if uneven spacing)
%
%     Fields for SAC xyz grid files:
%      data.x(:,1) - matrix data (nodes evenly spaced; advances l2r,b2t)
%     
%     Fields for SAC n-component timeseries files:
%      data.x(:,1) - 1st component
%      data.x(:,2) - 2nd component
%      data.x(:,n) - nth component
%      data.t - times (if uneven spacing)
%
%    Usage:    data=rdata(data)
%
%    By: Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: rh, rpdw, rsac, wsac, bsac, sachi, gh, lh, ch, wh, gv
%              sacsize

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'name') || ~isfield(data,'endian'))
    error('data structure does not have proper fields')
end

% number of records
nrecs=length(data);

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% initializing count of invalid SAC files
j=0;

% read loop
for i=1:nrecs
    % open file for reading
    fid=fopen(data(i-j).name,'r',data(i-j).endian);
    
    % invalid fid check (non-existent file or directory)
    if(fid<0)
        warning('SAClab:badFID','File not openable, %s',data(i-j).name);
        
        % delete and jump to next file
        data(i-j)=[];
        j=j+1;
        continue;
    end
    
    % file size
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    
    % grab some header fields (no warnings about ncmp)
    warning('off','SAClab:fieldInvalid');
    [iftype,npts,leven,ncmp]=gh(data(i-j),'iftype','npts','leven','ncmp');
    warning('on','SAClab:fieldInvalid');
    
    % byte size check
    total=sacsize(data(i-j).version,iftype,npts,leven,ncmp);
    if(bytes~=total)
        % inconsistent size
        fclose(fid);
        warning('SAClab:fileBadSize','File size wrong, %s',data(i-j).name);
        
        % delete and jump to next file
        data(i-j)=[];
        j=j+1;
        continue;
    end
    
    % logical index of header info
    v=data(i-j).version==vers;
    
    % act by file type (any new filetypes will have to be added here)
    fseek(fid,h(v).data.startbyte,'bof');
    if(iftype==h(v).enum(1).val.itime)
        % time series file - amplitude and time
        data(i-j).x(:,1)=fread(fid,npts,['*' h(v).data.store]);
        
        % timing of amp data if uneven
        if(leven==h(v).false)
            data(i-j).t(:,1)=fread(fid,npts,['*' h(v).data.store]);
        end
    elseif(iftype==h(v).enum(1).val.irlim)
        % spectral file - real and imaginary
        data(i-j).x(:,1)=fread(fid,npts,['*' h(v).data.store]);
        data(i-j).x(:,2)=fread(fid,npts,['*' h(v).data.store]);
    elseif(iftype==h(v).enum(1).val.iamph)
        % spectral file - amplitude and phase
        data(i-j).x(:,1)=fread(fid,npts,['*' h(v).data.store]);
        data(i-j).x(:,2)=fread(fid,npts,['*' h(v).data.store]);
    elseif(iftype==h(v).enum(1).val.ixy)
        % general x vs y data (x is 'dependent')
        data(i-j).x(:,1)=fread(fid,npts,['*' h(v).data.store]);
        
        % independent data (if uneven)
        if(leven==h(v).false)
            data(i-j).t(:,1)=fread(fid,npts,['*' h(v).data.store]);
        end
    elseif(iftype==h(v).enum(1).val.ixyz)
        % general xyz (3D) grid - nodes are evenly spaced
        data(i-j).x(:,1)=fread(fid,npts,['*' h(v).data.store]);
    elseif(isfield(h(v).enum(1).val,'incmp') ...
            && iftype==h(v).enum(1).val.incmp)
        % multi-component time series file
        for k=1:ncmp
            data(i-j).x(:,k)=fread(fid,npts,['*' h(v).data.store]);
        end
        
        % independent data (if uneven)
        if(leven==h(v).false)
            data(i-j).t(:,1)=fread(fid,npts,['*' h(v).data.store]);
        end
    else
        % unknown filetype
        fclose(fid)
        warning('SAClab:iftypeBad','Bad filetype code: %d',iftype);
        
        % delete and jump to next file
        data(i-j)=[];
        j=j+1;
        continue;
    end
    
    % closing file
    fclose(fid);
end

end