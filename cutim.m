function [data,failed]=cutim(data,varargin)
%CUTIM    Cut SAClab data records in memory
%
%    Description: CUTIM(DATA,PDW) cuts the records in DATA to the window
%     limits defined by PDW and outputs the updated dataset.  Header fields
%     are updated to match the windowed data.  Works with unevenly sampled
%     data.
%
%     PDW is a set of several arguments of one of the following forms:
%                  (1) REF1,OFFSET1,REF2,OFFSET2
%                  (2) REF1,REF2,OFFSET2
%                  (3) REF1,OFFSET1,REF2
%                  (4) REF,OFFSET1,OFFSET2
%                  (5) OFFSET1,REF2,OFFSET2
%                  (6) OFFSET1,OFFSET2
%                  (7) REF1,REF2
%                  (8) REF1,OFFSET1
%                  (9) OFFSET1
%                 (10) REF1
%                 (11) 
%
%     that defines explicitly or implicitly the starting and stopping 
%     points of the window.  This syntax attempts to match that of SAC.
%
%     REF is a string that refers to a reference value for the independent
%     component (usually time).  It can be any valid numeric header field 
%     or 'z' or 'x' or 'n'.
%      'z' is the zero position for the record - 0.
%      'x' indicates that the following offset is a sample number (first 
%       sample is 1).  Note that this defaults to sample 0 without an 
%       offset.
%      'n' indicates that the following offset is the length of the window
%       in number of samples.  Note that this defaults to length 0 without
%       an offset.  Also note that this is only valid as a second REF.
%
%     OFFSET is a numeric value giving the offset to be added to REF to 
%     define the position of the window start/stop with respect to the 
%     reference.  The offset can be a vector of values, one for each record
%     to be read, to define different offsets for each record (note that 
%     REF can not be a vector of reference positions currently).
%
%     REF1,OFFSET1 define the starting position of the window and
%     REF2,OFFSET2 define the ending position.  If a REF is given without
%     an OFFSET (forms 2, 3, 7, 10), the OFFSET defaults to 0 (no offset).  
%     If an OFFSET is given without a REF (forms 5, 6, 9), the REF defaults
%     to 'z' (zero) unless the window is of form (4) in which case both 
%     offsets share the same reference position.  If no ending position 
%     information is given (forms 8, 9, 10), REF2,OFFSET2 default to 'e' 
%     and 0 (end of the record).  If no window parameters are given (form
%     11), the window defaults to the entire record ('b',0,'e',0) - ie no 
%     window.
%
%     CUTIM(...,'FILL',TRUE|FALSE) turns on/off the filling of data gaps to
%     allow records that don't extend for the entire window to do so by
%     padding them with zeros.  Does not work with unevenly sampled data.
%     This option is useful for operations that require records to have the
%     same length.  See option 'FILLER' to alter the fill value.  By 
%     default 'FILL' is false (no fill).
%
%     CUTIM(...,'FILLER',VALUE) changes the fill value to VALUE.  Adjusting
%     this value does NOT turn option 'FILL' to true.  By default 'FILLER'
%     is set to 0 (zero).
%
%     [DATA,FAILED]=CUTIM(...,'TRIM',TRUE|FALSE) turns on/off deleting
%     records that have no data (after cutting) as well as spectral and xyz
%     records (which are not supported by CUTIM).  Optional output FAILED 
%     gives the indices of these records (with respect to their position in
%     the input DATA).  By default 'TRIM' is set to true (deletes records).
%
%    Notes:
%     - Windowing of spectral and xyz files are not supported.  By default
%       they are deleted from DATA (see option 'TRIM' to change this 
%       behavior).
%     - Windows with a start position after an end position will return
%       empty records (with headers updated accordingly) rather than
%       returning an error.
%     - Multiple component data is supported.
%     - FILL only works with evenly sampled data.
%
%    System requirements: Matlab 7
%
%    Data requirements: Time Series and General X vs Y data only
%     
%    Header changes: B, E, NPTS, DELTA, ODELTA
%                    DEPMEN, DEPMIN, DEPMAX
%
%    Usage: data=cutim(data,pdw)
%           data=cutim(data,pdw,'fill',true|false)
%           data=cutim(data,pdw,'fill',true|false,'filler',value)
%           [data,failed]=cutim(data,pdw,...'trim',true|false)
%
%    Examples:
%     cut a 400 sample window starting from the 33rd sample
%       data=cutim(data,'x',33,'n',400);
%
%     cut out a 90s window around t1
%       data=cutim(data,'t1',-30,'t1',60);
%
%     cut one hour starting at the origin time, 
%     padding incomplete records with zeros
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

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 23, 2008 - works with GLGC, GENUMDESC
%        Feb. 25, 2008 - bugfix
%        Feb. 28, 2008 - minor improvements
%        Feb. 29, 2008 - better leven check
%        Mar.  2, 2008 - dataless support
%        Mar.  4, 2008 - documentation update
%        Apr. 17, 2008 - PDW now like SAC
%        Apr. 18, 2008 - bugfix
%        May  12, 2008 - uses new dep* formula
%        June 12, 2008 - documentation update
%        June 23, 2008 - major documentation update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 23, 2008 at 22:45 GMT

% todo:
% multi ref support
% use chkhdr

% input check
error(nargchk(1,11,nargin))

% check data structure
error(seischk(data,'x'))

% parse cut parameters
[ref1,ref2,offset1,offset2,fill,filler,trim]=cutparam(varargin{:});

% number of records
nrecs=length(data);

% expand scalars
if(numel(offset1)==1); offset1(1:nrecs)=offset1; end
if(numel(offset2)==1); offset2(1:nrecs)=offset2; end

% make column vector
offset1=offset1(:);
offset2=offset2(:);

% cut parameter checks
if(numel(offset1)~=nrecs || numel(offset2)~=nrecs)
    error('SAClab:cutim:badInputSize',...
        'Number of elements in OFFSET not correct')
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
failed=false(nrecs,1);
for i=1:nrecs
    % check for unsupported filetypes
    if(strcmpi(iftype(i),'General XYZ (3-D) file'))
        failed(i)=true;
        warning('SAClab:cutim:illegalFiletype',...
            'illegal operation on xyz file');
        continue;
    elseif(any(strcmpi(iftype(i),{'Spectral File-Real/Imag'...
            'Spectral File-Ampl/Phase'})))
        failed(i)=true;
        warning('SAClab:cutim:illegalFiletype',...
            'illegal operation on spectral file');
        continue;
    end
    
    % evenly spaced
    if(strcmpi(leven(i),'true'))
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
            
            % empty window - add to failed list
            if(isempty(data(i).x))
                data(i)=ch(data(i),'b',0,'e',0,'npts',0,'delta',0,...
                    'depmen',0,'depmin',0,'depmax',0);
                failed(i)=true;
                continue;
            end
            
            % fix header
            data(i)=ch(data(i),'b',b(i)+(bp(i)-1)*delta(i),...
                'e',b(i)+(ep2-1)*delta(i),'npts',size(data(i).x,1),...
                'depmen',mean(data(i).x(:)),...
                'depmin',min(data(i).x(:)),...
                'depmax',max(data(i).x(:)));
        else
            % empty window - add to failed list
            if(isempty(data(i).x))
                data(i)=ch(data(i),'b',0,'e',0,'npts',0,'delta',0,...
                    'depmen',0,'depmin',0,'depmax',0);
                failed(i)=true;
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
                % empty window - add to failed list
                data(i).x=data(i).x([],:); data(i).t=[];
                data(i)=ch(data(i),'b',0,'e',0,'npts',0,'delta',0,...
                    'depmen',0,'depmin',0,'depmax',0,'odelta',0);
                failed(i)=true;
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
                % empty window - add to failed list
                data(i).x=data(i).x([],:); data(i).t=[];
                data(i)=ch(data(i),'b',0,'e',0,'npts',0,'delta',0,...
                    'depmen',0,'depmin',0,'depmax',0,'odelta',0);
                failed(i)=true;
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
            failed(i)=true;
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

% removed failed/empty cut records
if(trim); data(failed)=[]; end

end
