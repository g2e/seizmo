function []=wsac(data)
%WSAC    Write SAC binary files
%
%    Description: Writes SAC (seismic analysis code) binary format files
%     from a SAClab data structure.  Uses the 'name' field for naming
%     output files.  To write files to a specific directory either
%     add the path to the name field or browse Matlab there using the 'cd' 
%     command.  Output file's endianness can be set using the 'endian' 
%     field.
%
%    Usage:    wsac(saclab_struct)
%
%    See also:  rsac, bsac, sachi, gv, rpdw, rh, wh, lh, ch, gh

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'name') || ~isfield(data,'endian'))
    error('data structure does not have proper fields')
end

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% loop over records
for i=1:length(data)
    % logical index of header info
    v=(data(i).version==vers);
    
    % open file for writing
    fid=fopen(data(i).name,'w',data(i).endian);
    
    % fill header with dummy bytes (so we can seek around)
    fseek(fid,0,'bof');
    fwrite(fid,zeros(h(v).data.startbyte,1),'char');
    
    % write header
    n=h(v).types;
    for m=1:length(n)
        for k=1:length(h(v).(n{m}))
            fseek(fid,h(v).(n{m})(k).startbyte,'bof');
            fwrite(fid,data(i).head(h(v).(n{m})(k).minpos:h(v).(n{m})(k).maxpos),h(v).(n{m})(k).store);
        end
    end
    
    % get filetype and leven to write data
    [iftype,leven]=gh(data(i),'iftype','leven');
    
    % act by file type (another reason for enum(1) lookup only)
    fseek(fid,h(v).data.startbyte,'bof');
    if(iftype==h(v).enum(1).val.itime)
        % time series file - amplitude and time
        fwrite(fid,data(i).x(:,1),h(v).data.store);
        
        % timing of amp data if uneven
        if(leven==h(v).false)
            fwrite(fid,data(i).t(:,1),h(v).data.store);
        end
    elseif(iftype==h(v).enum(1).val.irlim)
        % spectral file - real and imaginary
        fwrite(fid,data(i).x(:,1),h(v).data.store);
        fwrite(fid,data(i).x(:,2),h(v).data.store);
    elseif(iftype==h(v).enum(1).val.iamph)
        % spectral file - amplitude and phase
        fwrite(fid,data(i).x(:,1),h(v).data.store);
        fwrite(fid,data(i).x(:,2),h(v).data.store);
    elseif(iftype==h(v).enum(1).val.ixy)
        % general x vs y data (y is 'dependent')
        fwrite(fid,data(i).x(:,1),h(v).data.store);
        
        % x data (if uneven)
        if(leven==h(v).false)
            fwrite(fid,data(i).t(:,1),h(v).data.store);
        end
    elseif(iftype==h(v).enum(1).val.ixyz)
        % general xyz (3D) - x,y are evenly spaced
        fwrite(fid,data(i).x(:,1),h(v).data.store);
    elseif(isfield(h(v).enum(1).val,'incmp') ...
            && iftype==h(v).enum(1).val.incmp)
        % multi-component time series file
        ncmp=gh(data(i),'ncmp');
        for k=1:ncmp
            fwrite(fid,data(i).x(:,k),h(v).data.store);
        end
        
        % independent data (if uneven)
        if(leven==h(v).false)
            fwrite(fid,data(i).t(:,1),h(v).data.store);
        end
    else
        % unknown filetype - close file and throw error
        fclose(fid);
        error('Unknown filetype code: %d',iftype);
    end
    
    % close file
    fclose(fid);
end

end
