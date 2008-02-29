function [data]=cutim(data,ref1,offset1,ref2,offset2,fill,filler)
%CUTIM    Cut seislab data records in memory
%
%    Description: Cuts data windows from seislab data records using the
%     given cut parameters.  Windows that extend outside the timing of the 
%     data are pushed to the data bounds unless 'fill' is set to 1. In this
%     case - and only for evenly spaced files - the window is padded with
%     filler values (default is zero) to meet the window parameters. If a
%     record ends up having no data it is deleted from the data structure.
%
%     ref values can be one of the following:
%        'z','n','some other header timing field', and for ref2 only: 'l'
%
%        'z' corresponds to the zero time of the record. 'n' indicates the 
%        offset is the sample number (first sample is 1).  ref2 can also
%        be 'l' which gives the window length in samples.  All header field
%        references (except 'n' and 'l') take offsets in seconds only.
%
%        If no ref value is given and an offset is, the ref is assumed as 
%        'z' (zero).  Offsets default to 0 when not supplied.  If neither 
%        a ref or offset are given for a pair, then they default to the
%        beginning and end - 'b',0 or 'e',0.  If both pairs are missing, no
%        cut is done.
%
%    Usage: [data]=cutim(data,ref1,offset1,ref2,offset2,fill,filler)
%
%    Examples:
%     cut a 400 sample window starting from the 33rd sample
%       data=cutim(data,'n',33,'l',400);
%
%     cut out a 90s window around t1
%       data=cutim(data,'t1',-30,'t1',60);
%
%     cut one hour starting at the origin time, 
%     and pad incomplete records with zeros
%       data=cutim(data,'o',0,'o',3600,1);
%
%     cut first 300 seconds of records
%       data=cutim(data,[],[],'b',300);
%     or more explicitly
%       data=cutim(data,'b',0,'b',300);
%
%     cut from 300 to 500 seconds relative to reference time
%       data=cutim(data,[],300,[],500);
%     or explicitly
%       data=cutim(data,'z',300,'z',500);
%
%    See also: rpdw

% input check
error(nargchk(1,7,nargin))

% check data structure
error(seischk(data,'x'))

% defaults
if(nargin==1); return; end
if(nargin<7); filler=0; end
if(nargin<6); fill=0; end
if(nargin<5); offset2=0; end
if(nargin<4); ref2='e'; end
if(nargin<3); offset1=0; end

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
    error('seislab:cutim:badInput','ref must be a string')
elseif(~isnumeric(offset1) || ~isnumeric(offset2) || ...
        ~isvector(offset1) || ~isvector(offset2))
    error('seislab:cutim:badInput','offset must be a numeric vector')
elseif(~any(length(offset1)==[1 nrecs]) || ~any(length(offset2)==[1 nrecs]))
    error('seislab:cutim:badInputSize','offset dimensions not correct')
end

% grab spacing info
iftype=genumdesc(data,'iftype');
leven=glgc(data,'leven');
[b,npts,delta]=gh(data,'b','npts','delta');

% check sample spacing logical
t=strcmp(leven,'true');
f=strcmp(leven,'false');
if(~all(t | f))
    error('sieslab:cutim:levenBad',...
        'logical field leven needs to be set'); 
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
    % check for unsupported filetypes
    if(strcmp(iftype(i),'General XYZ (3-D) file'))
        destroy(i)=1;
        warning('seislab:cutim:illegalFiletype',...
            'illegal operation on xyz file');
        continue;
    elseif(any(strcmp(iftype(i),{'Spectral File-Real/Imag'...
            'Spectral File-Ampl/Phase'})))
        destroy(i)=1;
        warning('seislab:cutim:illegalFiletype',...
            'illegal operation on spectral file');
        continue;
    end
    
    % evenly spaced
    if(strcmp(leven(i),'true'))
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
        
        % cut
        data(i).x=data(i).x(nbp:nep,:);
        
        % fill or no
        if(fill)
            % add filler
            cmp=size(data(i).x,2);
            data(i).x=[ones(1-bp(i),cmp)*filler; 
                        data(i).x; 
                        ones(ep2-npts(i),cmp)*filler];
            
            % empty window - add to destroy list
            if(isempty(data(i).x))
                destroy(i)=1;
                continue;
            end
            
            % fix header
            data(i)=ch(data(i),'b',b(i)+(bp(i)-1)*delta(i),...
                'e',b(i)+(ep2-1)*delta(i),'npts',size(data(i).x,1),...
                'depmen',norm(mean(data(i).x)),...
                'depmin',-norm(min(data(i).x)),...
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
                'depmen',norm(mean(data(i).x)),...
                'depmin',-norm(min(data(i).x)),...
                'depmax',norm(max(data(i).x)));
        end
    % unevenly spaced
    else
        % get begin point
        if(~strcmpi(ref1,'n'))
            % save to temp variable (avoids corruption in no result case)
            temp=find(data(i).t>=bt(i),1);
            if(isempty(temp))
                % empty window - add to destroy list
                destroy(i)=1;
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
        if(strcmpi(ref2,'l'))
            ep2=bp(i)+ep(i)-1;
        elseif(strcmpi(ref2,'n'))
            ep2=ep(i);
        else
            % save to temp variable (avoids corruption in no result case)
            temp=find(data(i).t<=et(i),1,'last');
            if(isempty(temp))
                % empty window - add to destroy list
                destroy(i)=1;
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
        
        % check if bp>ep
        if(nbp>nep)
            destroy(i)=1;
            continue;
        end
        
        % cut
        data(i).x=data(i).x(nbp:nep,:);
        data(i).t=data(i).t(nbp:nep,1);
        
        % fix header
        data(i)=ch(data(i),'b',data(i).t(1),...
            'e',data(i).t(end),'npts',size(data(i).x,1),...
            'depmen',norm(mean(data(i).x)),...
            'depmin',-norm(min(data(i).x)),...
            'depmax',norm(max(data(i).x)),...
            'delta',(data(i).t(end)-data(i).t(1))/(npts(i)-1),...
            'odelta',data(i).t(min([end 2]))-data(i).t(1));
    end
end

% destroy empty records
data(destroy)=[];

end
