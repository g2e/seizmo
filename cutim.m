function [data,destroy]=cutim(data,varargin)
%CUTIM    Cut SAClab data records in memory
%
%    Description: Cuts data windows from SAClab data records using the
%     given cut parameters.  Windows that extend outside the timing of the 
%     data are pushed to the data bounds unless 'fill' is set to true/1. In
%     this case - and only for evenly spaced files - the window is padded
%     with filler values (default is zero) to meet the window parameters.
%     If a record ends up having no data or is not a timeseries filetype,
%     it is deleted from the data structure unless the 'removedataless'
%     parameter is set to false.  
%
%     Defining the window is done by supplying time references and/or
%     offsets in a manner just like SAC's cut command.
%
%     ref and offset values:
%        ref can be 'z','x','some other header timing field', and for ref2 
%        only: 'n'.  'z' corresponds to the zero time of the record. 'x' 
%        indicates the following offset is the sample number (first sample 
%        is 1).  ref2 can be 'n' which gives the window length in samples.
%        All header field references (except 'x' and 'n') assume offsets
%        are in seconds.
%
%        If no ref value is given and an offset is, the ref is assumed as 
%        'z' (zero).  Offsets default to 0 when not supplied.  If neither 
%        the second ref or offset are given, then they default to the end -
%        'e',0.  If both pairs are missing, then they default to 'b',0 and
%        'e',0 which under normal situations will not cut the records.
%
%    Usage: [data]=cutim(data,...,variable,list,of,cut,parameters,...)
%
%    Examples:
%     cut a 400 sample window starting from the 33rd sample
%       data=cutim(data,'x',33,'n',400);
%
%     cut out a 90s window around t1
%       data=cutim(data,'t1',-30,'t1',60);
%
%     cut one hour starting at the origin time, 
%     and pad incomplete records with zeros
%       data=cutim(data,'o',0,3600,'fill',1);
%
%     cut first 300 seconds of records
%       data=cutim(data,'b',0,300);
%
%     cut from 300 to 500 seconds relative to reference time
%       data=cutim(data,300,500);
%     or explicitly
%       data=cutim(data,'z',300,'z',500);
%
%    See also: rpdw

% input check
error(nargchk(1,11,nargin))

% check data structure
error(seischk(data,'x'))

% parse cut parameters
[ref1,ref2,offset1,offset2,fill,filler,trim]=cutparam(varargin{:});

% number of records
nrecs=length(data);

% cut parameter checks
if(~ischar(ref1) || ~ischar(ref2))
    error('SAClab:cutim:badInput','ref must be a string')
elseif(~isvector(offset1) || ~isvector(offset2))
    error('SAClab:cutim:badInput','offset must be a numeric vector')
elseif(~any(length(offset1)==[1 nrecs]) || ~any(length(offset2)==[1 nrecs]))
    error('SAClab:cutim:badInputSize','offset dimensions not correct')
end

% grab spacing info
iftype=genumdesc(data,'iftype');
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
warning('off','SAClab:gh:fieldInvalid')
[b,npts,delta,ncmp]=gh(data,'b','npts','delta','ncmp');
warning('on','SAClab:gh:fieldInvalid')

% clean up and check ncmp
ncmp(isnan(ncmp))=1;
if(any(ncmp<1 | fix(ncmp)~=ncmp))
    error('SAClab:cutim:badNumCmp',...
        'field ncmp must be a positive integer')
end

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
    error('SAClab:cutim:nanHeaderField','header field returned NaN')
end

% loop through each file
destroy=false(nrecs,1);
for i=1:nrecs
    % check for unsupported filetypes
    if(strcmp(iftype(i),'General XYZ (3-D) file'))
        destroy(i)=true;
        warning('SAClab:cutim:illegalFiletype',...
            'illegal operation on xyz file');
        continue;
    elseif(any(strcmp(iftype(i),{'Spectral File-Real/Imag'...
            'Spectral File-Ampl/Phase'})))
        destroy(i)=true;
        warning('SAClab:cutim:illegalFiletype',...
            'illegal operation on spectral file');
        continue;
    end
    
    % evenly spaced
    if(strcmp(leven(i),'true'))
        % calculate begin and end points
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
        
        % cut
        data(i).x=data(i).x(nbp:nep,:);
        
        % fill or no
        if(fill)
            % add filler
            data(i).x=[ones(1-bp(i),ncmp(i))*filler; 
                        data(i).x; 
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
    % unevenly spaced
    else
        % get begin point
        if(~strcmpi(ref1,'x'))
            % save to temp variable (avoids corruption in no result case)
            temp=find(data(i).t>=bt(i),1);
            if(isempty(temp))
                % empty window - add to destroy list
                data(i).x=data(i).x([],:); data(i).t=[];
                data(i)=ch(data(i),'b',0,'e',0,'npts',0,'delta',0,...
                    'depmen',0,'depmin',0,'depmax',0,'odelta',0);
                destroy(i)=true;
                continue;
            elseif(temp>1)
                temp2=temp-1;
                
                % figure out which is closer
                if(abs(bt(i)-data(i).t(temp))<=abs(bt(i)-data(i).t(temp2)))
                    bp(i)=temp;
                else
                    bp(i)=temp2;
                end
            else
                bp(i)=temp;
            end
        end
        
        % get end point
        if(strcmpi(ref2,'n'))
            ep2=bp(i)+ep(i)-1;
        elseif(strcmpi(ref2,'x'))
            ep2=ep(i);
        else
            % save to temp variable (avoids corruption in no result case)
            temp=find(data(i).t<=et(i),1,'last');
            if(isempty(temp))
                % empty window - add to destroy list
                data(i).x=data(i).x([],:); data(i).t=[];
                data(i)=ch(data(i),'b',0,'e',0,'npts',0,'delta',0,...
                    'depmen',0,'depmin',0,'depmax',0,'odelta',0);
                destroy(i)=true;
                continue;
            elseif(temp<npts(i))
                temp2=temp+1;
                
                % figure out which is closer
                if(abs(et(i)-data(i).t(temp))<=abs(et(i)-data(i).t(temp2)))
                    ep2=temp;
                else
                    ep2=temp2;
                end
            else
                ep2=temp;
            end
        end
        
        % boundary conditions
        nbp=max([bp(i) 1]);
        nep=min([ep2 npts(i)]);
        
        % cut
        data(i).x=data(i).x(nbp:nep,:);
        data(i).t=data(i).t(nbp:nep,1);
        
        % check if bp>ep
        if(nbp>nep)
            data(i)=ch(data(i),'b',0,'e',0,'npts',0,'delta',0,...
                'depmen',0,'depmin',0,'depmax',0,'odelta',0);
            destroy(i)=true;
            continue;
        end
        
        % fix header
        data(i)=ch(data(i),'b',data(i).t(1),...
            'e',data(i).t(end),'npts',size(data(i).x,1),...
            'depmen',mean(data(i).x(:)),...
            'depmin',min(data(i).x(:)),...
            'depmax',max(data(i).x(:)),...
            'delta',(data(i).t(end)-data(i).t(1))/(npts(i)-1),...
            'odelta',data(i).t(min([end 2]))-data(i).t(1));
    end
end

% destroy empty/wrong filetype records
if(trim); data(destroy)=[]; end

end
