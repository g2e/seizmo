function []=wh(data)
%WH    Write SAC binary file headers
%
%    Description: Writes out SAC (seismic analysis code) binary format file
%     headers from a SAClab data structure.  Basically this updates the
%     headers of existing files to match the SAClab data structure.  For
%     files that do not already exist the 'name' field is used for output
%     file naming.  To write files to a specific directory either add the 
%     path to the 'name' field or browse Matlab there using the 'cd' 
%     command.  Output file byte-order can be set using the 'endian' field.
%
%    Warning:  Using wh after modifying the version/byte-order in SAClab
%              could cause the file to become corrupted (eg. header has
%              one byte-order while the data has the other).  Please use
%              wsac for modifying version/byte-order to avoid such issues.
%
%    Usage:    wh(saclab_struct)
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also:  wsac, rsac, bsac, rpdw, rdata, rh, lh, ch, gh, sachp, gv

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
    v=data(i).version==vers;
    
    % open existing file for writing
    if(exist(data(i).name,'file'))
        % check for version/byte-order change
        [version,endian]=gv(data(i).name);
        if(version)
            if(version~=data(i).version)
                warning('SAClab:versionMismatch',['Version of existing'...
                    'file %s does not match output header version.  '...
                    'Accurate reading of data in existing file '...
                    'may not be possible now!'],data(i).name);
            end
            if(endian~=data(i).endian)
                warning('SAClab:endianMismatch',['Byte-order of existing'...
                    'file %s does not match output header byte-order.  '...
                    'Accurate reading of data in existing file may not '...
                    'be possible now!'],data(i).name);
            end
        end
        
        % open file for modification
        fid=fopen(data(i).name,'r+',data(i).endian);
    else
        % new file
        fid=fopen(data(i).name,'w',data(i).endian);
    end
    
    % fill header with zeros (so we can seek around)
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
    
    % close file
    fclose(fid);
end

end
