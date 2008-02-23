function [data]=rpdw(data,ref1,offset1,ref2,offset2,fill,filler)
%RPDW    Read partial data window from SAC binary file
%
%    Description: Reads a partial data window from a SAC (seismic analysis
%     code) binary file using the reference and offset parameters given.
%     This provides a mechanism similar to the SAC 'cut' command to limit
%     memory usage.  Input parameters are exactly the same as the SAClab 
%     'cutim' command - for details on cutting read there.
%
%    Usage: [data]=rpdw(data,ref1,offset1,ref2,offset2,fill,filler)
%
%    See also: cutim, rh, wh, ch, lh, gh, rsac, wsac, bsac, sachi, gv,
%              rdata

% input check
error(nargchk(1,7,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'name') || ~isfield(data,'endian'))
    error('data structure does not have proper fields')
end

% defaults
if(nargin<7); filler=0; end
if(nargin<6); fill=0; end
if(nargin<5); offset2=0; end
if(nargin<4); ref2='e'; end
if(nargin<3); offset1=0; end
if(nargin<2); ref1='b'; end

% empty cut parameter defaults
if(isempty(ref1)); if(isempty(offset1)); ref1='b'; else ref1='z'; end; end
if(isempty(ref2)); if(isempty(offset2)); ref2='e'; else ref2='z'; end; end
if(isempty(offset1)); offset1=0; end
if(isempty(offset2)); offset2=0; end
if(isempty(fill)); fill=0; end
if(isempty(filler)); filler=0; end

% number of seismograms
nrecs=length(data);

% cut parameter checks
if(~ischar(ref1) || ~ischar(ref2))
    error('references must be strings')
elseif(~isnumeric(offset1) || ~isnumeric(offset2) || ...
        ~isvector(offset1) || ~isvector(offset2))
    error('offsets must be a numeric vector')
elseif(~any(length(offset1)==[1 nrecs]) || ~any(length(offset2)==[1 nrecs]))
    error('offsets dimensions not correct')
end

% grab header info
[b,npts,delta,leven,iftype]=gh(data,'b','npts','delta','leven','iftype');
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% expand scalar offsets
if(length(offset1)==1)
    offset1=offset1(ones(nrecs,1));
end
if(length(offset2)==1)
    offset2=offset2(ones(nrecs,1));
end

% window begin point
if(strcmpi(ref1,'z'))
    bt(1:nrecs,1)=offset1;
elseif(strcmpi(ref1,'n'))
    bp=round(offset1);
else
    bt=gh(data,ref1)+offset1;
end

% window end time
if(strcmpi(ref2,'z'))
    et(1:nrecs,1)=offset2;
elseif(strcmpi(ref2,'n') || strcmpi(ref2,'l'))
    ep=round(offset2);
else
    et=gh(data,ref2)+offset2;
end

% loop through each file
destroy=false(nrecs,1);
for i=1:nrecs
    % header version index
    v=data(i).version==vers;
    
    % check for unsupported filetypes
    if(iftype(i)==h(v).enum(1).val.ixyz)
        destroy(i)=1;
        warning('SAClab:illegalFiletype','Illegal operation on xyz file');
        continue;
    elseif(iftype(i)==h(v).enum(1).val.irlim || iftype(i)==h(v).enum(1).val.iamph)
        destroy(i)=1;
        warning('SAClab:illegalFiletype','Illegal operation on spectral file');
        continue;
    end
    
    % check if evenly spaced
    if(leven(i)==h(v).true)
        % calculate begin and end points
        if(~strcmpi(ref1,'n'))
            bp(i)=round((bt(i)-b(i))/delta(i))+1;
        end
        if(strcmpi(ref2,'l'))
            ep2=bp(i)+ep(i)-1;
        elseif(strcmpi(ref2,'n'))
            ep2=ep(i);
        else
            ep2=round((et(i)-b(i))/delta(i))+1;
        end
        
        % boundary conditions
        nbp=max([bp(i) 1]);
        nep=min([ep2 npts(i)]);
        
        % open file
        fid=fopen(data(i).name,'r',data(i).endian);
        
        % check that it opened
        if(fid<0)
            warning('SAClab:fileNoExist','File does not exist: %s',data(i).name);
            destroy(i)=1;
            continue;
        end
        
        % move to first byte of window
        try
            fseek(fid,h(v).data.startbyte+h(v).data.bytesize*(nbp-1),'bof');
        catch
            fclose(fid);
            warning('SAClab:seekFailed','Seek to data failed: %s',data(i).name);
            destroy(i)=1;
            continue;
        end
        
        % read in data
        data(i).x(:,1)=fread(fid,nep-nbp+1,['*' h(v).data.store]);
        
        % verify read
        if(length(data(i).x(:,1))~=nep-nbp+1)
            fclose(fid);
            warning('SAClab:readFailed','Reading data failed: %s',data(i).name);
            destroy(i)=1;
            continue;
        end
        
        % incmp filetype reading continues here
        if(isfield(h(v).enum(1).val,'incmp') && iftype==h(v).enum(1).val.incmp)
            % read and check ncmp
            ncmp=gh(data(i),'ncmp');
            if(ncmp<1 || ncmp-fix(ncmp)~=0)
                warning('SAClab:ncmpInvalid','Number of components invalid: %s',data(i).name);
            end
            
            % loop through components
            for j=2:ncmp
                % move to first byte of window
                try
                    fseek(fid,h(v).data.startbyte+h(v).data.bytesize*(j*npts(i)+nbp-1),'bof');
                catch
                    fclose(fid);
                    warning('SAClab:seekFailed','Seek to data failed: %s',data(i).name);
                    destroy(i)=1;
                    continue;
                end
                
                % read in data
                temp=fread(fid,nep-nbp+1,['*' h(v).data.store]);
                
                % verify read
                if(length(temp)~=nep-nbp+1)
                    fclose(fid);
                    warning('SAClab:readFailed','Reading data failed: %s',data(i).name);
                    destroy(i)=1;
                    continue;
                end
                
                % assign data
                data(i).x(:,j)=temp;
            end
        end
        
        % close file
        fclose(fid);
        
        % fill or no
        if(fill)
            % add filler
            cmp=size(data(i).x,2);
            data(i).x=[ones(1-bp(i),cmp)*filler; data(i).x; ones(ep2-npts(i),cmp)*filler];
            
            % empty window - add to destroy list
            if(isempty(data(i).x))
                destroy(i)=1;
                continue;
            end
            
            % fix header
            data(i)=ch(data(i),'b',b(i)+(bp(i)-1)*delta(i),...
                'e',b(i)+(ep2-1)*delta(i),'npts',size(data(i).x,1),...
                'depmen',norm(mean(data(i).x)),'depmin',-norm(min(data(i).x)),...
                'depmax',norm(max(data(i).x)));
        else
            % empty window - add to destroy list
            if(isempty(data(i).x))
                destroy(i)=1;
                continue;
            end
                    
            % fix header
            data(i)=ch(data(i),'b',b(i)+(nbp-1)*delta(i),...
                'e',b(i)+(nep-1)*delta(i),'npts',size(data(i).x,1),...
                'depmen',norm(mean(data(i).x)),'depmin',-norm(min(data(i).x)),...
                'depmax',norm(max(data(i).x)));
        end
    % unevenly spaced
    else
        data(i)=rdata(data(i));
        data(i)=cutim(data(i),ref1,offset1,ref2,offset2,fill,filler);
    end
end

% destroy empty records
data(destroy)=[];

end
