function [data,destroy]=rpdw(data,varargin)
%RPDW    Read partial data window from binary seismic datafiles in SAClab
%
%    Description: Reads a partial data window from binary seismic datafiles
%     into a SAClab structure using the windowing parameters supplied.  
%     This provides a mechanism similar to the SAC 'cut' command to limit 
%     memory/cpu usage related to reading in large datafiles.  Input 
%     parameters are exactly the same as the SAClab 'cutim' command - for 
%     details on cutting read there.
%
%    Usage: [data]=rpdw(data,...,variable,list,of,cut,parameters,...)
%
%    Examples:
%     remember to read in the header info first:
%      files=dir('*.sac')  % reads in files ending in .sac from current dir
%      data=rh(files.name) % reads in headers of sac files
%
%     read in only the first 300 samples:
%      data=rpdw(data,'x',1,300)
%     
%     read in a 90 second window around t1 arrival, padding with zeros if
%     necessary:
%      data=rpdw(data,'t1',-30,60,'fill',true)
%
%     read in the 123rd sample:
%      data=rpdw(data,'x',123,123)
%     or
%      data=rpdw(data,'x',123,'n',1)
%
%    See also: cutim, rh, rdata, rseis, wh, wseis

% input check
error(nargchk(1,11,nargin))

% check data structure
error(seischk(data,'name','endian'))

% parse cut parameters
[ref1,ref2,offset1,offset2,fill,filler,trim]=cutparam(varargin{:});

% number of records
nrecs=length(data);

% cut parameter checks
if(~ischar(ref1) || ~ischar(ref2))
    error('SAClab:rpdw:badInput','ref must be a string')
elseif(~isvector(offset1) || ~isvector(offset2))
    error('SAClab:rpdw:badInput','offset must be a numeric vector')
elseif(~any(length(offset1)==[1 nrecs]) || ...
        ~any(length(offset2)==[1 nrecs]))
    error('SAClab:rpdw:badInputSize','offset dimensions not correct')
end

% grab header info
iftype=genumdesc(data,'iftype');
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
warning('off','SAClab:gh:fieldInvalid')
[b,npts,delta,ncmp]=gh(data,'b','npts','delta','ncmp');
warning('on','SAClab:gh:fieldInvalid')

% clean up and check ncmp
ncmp(isnan(ncmp))=1;
if(any(ncmp<1 | fix(ncmp)~=ncmp))
    error('SAClab:rpdw:badNumCmp',...
        'field ncmp must be a positive integer')
end

% indices based on leven
tru=strcmp(leven,'true');
fals=strcmp(leven,'false');
trui=find(tru);
falsi=find(fals);

% column vector/expand scalar offsets
offset1=offset1(:);
offset2=offset2(:);
if(length(offset1)==1)
    offset1=offset1(ones(nrecs,1));
end
if(length(offset2)==1)
    offset2=offset2(ones(nrecs,1));
end

% window begin point
bt=[]; bp=[];
if(strcmpi(ref1,'z'))
    bt=offset1;
elseif(strcmpi(ref1,'x'))
    bp=round(offset1);
else
    bt=gh(data,ref1)+offset1;
end

% window end time
et=[]; ep=[];
if(strcmpi(ref2,'z'))
    et=offset2;
elseif(strcmpi(ref2,'x') || strcmpi(ref2,'n'))
    ep=round(offset2);
else
    et=gh(data,ref2)+offset2;
end

% check for nans
if(any(isnan(bt)) || any(isnan(bp)) || any(isnan(et)) || any(isnan(ep)))
    error('SAClab:rpdw:nanHeaderField','header field returned NaN')
end

% grab header setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=seishi(vers(nver));
for i=1:nver-1
    h(i)=seishi(vers(i));
end

% allocate bad records matrix
destroy=false(nrecs,1);

% let rdata/cutim handle unevenly sampled records minus file deletion
if(~isempty(falsi))
    [data(falsi),destroy(falsi)]=rdata(data(falsi),false);
    if(~isempty(falsi(~destroy(falsi))))
        [data(falsi(~destroy(falsi))),destroy(falsi(~destroy(falsi)))]=...
            cutim(data(falsi(~destroy(falsi))),ref1,offset1,ref2,offset2,...
            'fill',fill,'filler',filler,'removedataless',false);
    end
end

% loop through each file
for i=trui.'
    % header version index
    v=(data(i).version==vers);
    
    % check for unsupported filetypes
    if(strcmp(iftype(i),'General XYZ (3-D) file'))
        destroy(i)=true;
        warning('SAClab:rpdw:illegalFiletype',...
            'illegal operation on xyz file');
        continue;
    elseif(any(strcmp(iftype(i),{'Spectral File-Real/Imag'...
            'Spectral File-Ampl/Phase'})))
        destroy(i)=true;
        warning('SAClab:rpdw:illegalFiletype',...
            'illegal operation on spectral file');
        continue;
    end
    
    % closest points to begin and end
    if(~strcmpi(ref1,'x'))
        bp(i)=round((bt(i)-b(i))/delta(i))+1;
    end
    if(strcmpi(ref2,'n'))
        ep2=bp(i)+ep(i)-1;
    elseif(strcmpi(ref2,'x'))
        ep2=ep(i);
    else
        ep2=round((et(i)-b(i))/delta(i))+1;
    end
    
    % boundary conditions
    nbp=max([bp(i) 1]);
    nep=min([ep2 npts(i)]);
    nnp=nep-nbp+1;
    
    % open file
    fid=fopen(data(i).name,'r',data(i).endian);
    
    % check that it opened
    if(fid<0)
        warning('SAClab:rpdw:badFID',...
            'File not openable, %s',data(i).name);
        destroy(i)=true;
        continue;
    end
    
    % preallocate data record with NaNs, deallocate timing
    data(i).x=nan(nnp,ncmp(i),h(v).data.store);
    data(i).t=[];
    
    % skip if new npts==0
    if(nnp<1)
        data(i)=ch(data(i),'b',0,'e',0,'npts',0,'delta',0,...
            'depmen',0,'depmin',0,'depmax',0);
        destroy(i)=true;
        continue;
    end
    
    % loop through each component
    for j=1:ncmp(i)
        % move to first byte of window and read
        try
            fseek(fid,h(v).data.startbyte+...
                h(v).data.bytesize*((j-1)*npts(i)+nbp-1),'bof');
            data(i).x(:,j)=fread(fid,nnp,['*' h(v).data.store]);
        catch
            warning('SAClab:rpdw:readFailed',...
                'Read in of data failed: %s',data(i).name);
            destroy(i)=true;
            break;
        end
    end
    
    % close file
    fclose(fid);
    
    % skip to next if read failed
    if(destroy(i)); continue; end
    
    % fill or no
    if(fill)
        % add filler
        data(i).x=[ones(1-bp(i),ncmp(i))*filler; ...
            data(i).x; ...
            ones(ep2-npts(i),ncmp(i))*filler];
        
        % empty window - add to destroy list
        if(isempty(data(i).x))
            data(i)=ch(data(i),'b',0,'e',0,'npts',0,'delta',0,...
                'depmen',0,'depmin',0,'depmax',0);
            destroy(i)=true;
            continue;
        end
        
        % fix header
        data(i)=ch(data(i),'b',b(i)+(bp(i)-1)*delta(i),...
            'e',b(i)+(ep2-1)*delta(i),'npts',size(data(i).x,1),...
            'depmen',mean(data(i).x(:)),...
            'depmin',min(data(i).x(:)),...
            'depmax',max(data(i).x(:)));
    else
        % empty window - add to destroy list
        if(isempty(data(i).x))
            data(i)=ch(data(i),'b',0,'e',0,'npts',0,'delta',0,...
                'depmen',0,'depmin',0,'depmax',0);
            destroy(i)=true;
            continue;
        end
        
        % fix header
        data(i)=ch(data(i),'b',b(i)+(nbp-1)*delta(i),...
            'e',b(i)+(nep-1)*delta(i),'npts',size(data(i).x,1),...
            'depmen',mean(data(i).x(:)),...
            'depmin',min(data(i).x(:)),...
            'depmax',max(data(i).x(:)));
    end
end

% destroy empty records
if(trim); data(destroy)=[]; end

end
