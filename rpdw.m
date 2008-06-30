function [data,failed]=rpdw(data,varargin)
%RPDW    Read a partial data window from datafiles into SAClab
%
%    Description: RPDW(DATA,PDW) reads a partial data window PDW from 
%     datafiles on disk into the SAClab structure DATA and returns the 
%     updated structure.  This provides a mechanism similar to the SAC 
%     'cut' command to limit memory/cpu usage related to reading in large 
%     datafiles where only a segment is needed.  Note that header fields
%     will be updated to match the windowed data.
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
%     RPDW(...,'FILL',TRUE|FALSE) turns on/off the filling of data gaps to
%     allow records that don't extend for the entire window to do so by
%     padding them with zeros.  This is useful for operations that require
%     records to all have the same length.  See option 'FILLER' to alter 
%     the fill value.  By default 'FILL' is false (no fill).
%
%     RPDW(...,'FILLER',VALUE) changes the fill value to VALUE.  Adjusting
%     this value does NOT turn option 'FILL' to true.  By default 'FILLER'
%     is set to 0 (zero).
%
%     [DATA,FAILED]=RPDW(...,'TRIM',TRUE|FALSE) turns on/off deleting
%     records that failed to read or would have no data.  This includes
%     spectral and xyz files (which are not supported by RPDW).  Optional 
%     output FAILED gives the indices of these records (with respect to 
%     their position in the input DATA).  By default 'TRIM' is set to true
%     (deletes records).
%
%    Notes:
%     - Partial reads of spectral and xyz files are not supported.  They
%       will not be read in with RPDW.  By default they are deleted from
%       DATA (see option 'TRIM' to change this behavior).
%     - Windows with a start position after an end position will return
%       empty records (with headers updated accordingly) rather than
%       returning an error.
%     - Multiple component data is supported.
%     - Partial reads of unevenly sampled data is not supported.  They are
%       passed to RDATA and CUTIM instead.
%     - FILL only works with evenly sampled data.
%
%    System requirements: Matlab 7
%
%    Data requirements: DATA has header, endian, name and version fields.
%                       Time Series and General X vs Y data only.
%
%    Header changes: B, E, NPTS, DELTA, ODELTA
%                    DEPMEN, DEPMIN, DEPMAX
%
%    Usage: data=rpdw(data,pdw)
%           data=rpdw(data,pdw,'fill',true|false)
%           data=rpdw(data,pdw,'fill',true|false,'filler',value)
%           [data,failed]=rpdw(data,pdw,...'trim',true|false)
%
%    Examples:
%     ALL THE FOLLOWING EXAMPLES REQUIRE YOU TO 
%     REMEMBER TO READ IN THE HEADER INFO FIRST:
%      data=rh('*.SAC')
%
%     read in only the first 300 samples of records:
%      data=rpdw(data,'x',1,300)
%     
%     read in a 90 second window around the t1 arrival, 
%     padding data gaps with zeros as necessary:
%      data=rpdw(data,'t1',-30,60,'fill',true)
%
%     read in only the 123rd sample:
%      data=rpdw(data,'x',123,123)
%     or
%      data=rpdw(data,'x',123,'n',1)
%
%    See also: cutim, rh, rdata, rseis, wh, wseis

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 23, 2008 - works with GLGC, GENUMDESC
%        Mar.  2, 2008 - dataless support
%        Mar.  4, 2008 - documentation update
%        Apr. 17, 2008 - PDW now like SAC
%        Apr. 18, 2008 - bugfix
%        May  12, 2008 - uses new dep* formula
%        June 12, 2008 - documentation update
%        June 23, 2008 - major documentation update, bugfix for uneven
%                        files, some other cleanups for readability
%        June 30, 2008 - fixed dataless support, .dep & .ind rather than .x
%                        & .t, fix for single point data
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 23, 2008 at 18:45 GMT

% todo:
% multi ref support?
% use chkhdr?

% input check
error(nargchk(1,11,nargin))

% check data structure
error(seischk(data,'name','endian'))

% parse cut parameters
[ref1,ref2,offset1,offset2,fill,filler,trim]=cutparam(varargin{:});

% number of records
nrecs=numel(data);

% expand scalars
if(numel(offset1)==1); offset1(1:nrecs)=offset1; end
if(numel(offset2)==1); offset2(1:nrecs)=offset2; end

% make column vector
offset1=offset1(:);
offset2=offset2(:);

% cut parameter checks
if(numel(offset1)~=nrecs || numel(offset2)~=nrecs)
    error('SAClab:rpdw:badInputSize',...
        'Number of elements in OFFSET not correct')
end

% check leven
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
even=find(strcmpi(leven,'true')).';
uneven=strcmpi(leven,'false');

% get more header info
iftype=genumdesc(data,'iftype');
warning('off','SAClab:gh:fieldInvalid')
[b,npts,delta,ncmp]=gh(data,'b','npts','delta','ncmp');
warning('on','SAClab:gh:fieldInvalid')

% clean up and check ncmp
ncmp(isnan(ncmp))=1;
if(any(ncmp<1 | fix(ncmp)~=ncmp))
    error('SAClab:rpdw:badNumCmp',...
        'field ncmp must be a positive integer')
end

% window start point/value
bt=[]; bp=[];
if(strcmpi(ref1,'z'))
    bt=offset1;
elseif(strcmpi(ref1,'x'))
    bp=round(offset1);
else
    bt=gh(data,ref1)+offset1;
end

% window stop point/value
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
    error('SAClab:rpdw:nanHeaderField','Window position is NaN!')
end

% grab header setup (so we know how to read the file)
vers=unique([data.version]);
nver=length(vers);
h(nver)=seisdef(vers(nver));
for i=1:nver-1
    h(i)=seisdef(vers(i));
end

% allocate bad records matrix
failed=false(nrecs,1);

% let rdata/cutim handle unevenly sampled records (minus file deletion)
if(any(uneven))
    [data(uneven),failed(uneven)]=rdata(data(uneven),false);
    if(any(uneven & ~failed))
        uneven=(uneven & ~failed);
        [data(uneven),failed(uneven)]=cutim(data(uneven),...
            ref1,offset1(uneven),ref2,offset2(uneven),...
            'fill',fill,'filler',filler,'trim',false);
    end
end

% loop through each evenly spaced file
for i=even
    % header version index
    v=(data(i).version==vers);
    
    % check for unsupported filetypes
    if(strcmpi(iftype(i),'General XYZ (3-D) file'))
        failed(i)=true;
        warning('SAClab:rpdw:illegalFiletype',...
            'illegal operation on xyz file');
        continue;
    elseif(any(strcmpi(iftype(i),{'Spectral File-Real/Imag'...
            'Spectral File-Ampl/Phase'})))
        failed(i)=true;
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
        failed(i)=true;
        continue;
    end
    
    % preallocate data record with NaNs, deallocate timing
    data(i).dep=nan(nnp,ncmp(i),h(v).data.store);
    data(i).ind=[];
    
    % skip if new npts==0
    if(nnp<1)
        data(i)=ch(data(i),'b',nan,'e',nan,'npts',0,'delta',nan,...
            'depmen',nan,'depmin',nan,'depmax',nan);
        failed(i)=true;
        continue;
    end
    
    % loop through each component
    for j=1:ncmp(i)
        % move to first byte of window and read
        try
            fseek(fid,h(v).data.startbyte+...
                h(v).data.bytesize*((j-1)*npts(i)+nbp-1),'bof');
            data(i).dep(:,j)=fread(fid,nnp,['*' h(v).data.store]);
        catch
            warning('SAClab:rpdw:readFailed',...
                'Read in of data failed: %s',data(i).name);
            failed(i)=true;
            break;
        end
    end
    
    % close file
    fclose(fid);
    
    % skip to next record if read failed
    if(failed(i)); continue; end
    
    % to fill
    if(fill)
        % add filler
        data(i).dep=[ones(1-bp(i),ncmp(i))*filler; ...
            data(i).dep; ...
            ones(ep2-npts(i),ncmp(i))*filler];
        
        % empty window - add to failed list
        if(isempty(data(i).dep))
            data(i)=ch(data(i),'b',nan,'e',nan,'npts',0,'delta',nan,...
                'depmen',nan,'depmin',nan,'depmax',nan);
            failed(i)=true;
            continue;
        end
        
        % fix header
        data(i)=ch(data(i),'b',b(i)+(bp(i)-1)*delta(i),...
            'e',b(i)+(ep2-1)*delta(i),'npts',size(data(i).dep,1),...
            'depmen',mean(data(i).dep(:)),...
            'depmin',min(data(i).dep(:)),...
            'depmax',max(data(i).dep(:)));
    % not to fill
    else
        % empty window - add to failed list
        if(isempty(data(i).dep))
            data(i)=ch(data(i),'b',nan,'e',nan,'npts',0,'delta',nan,...
                'depmen',nan,'depmin',nan,'depmax',nan);
            failed(i)=true;
            continue;
        end
        
        % fix header
        data(i)=ch(data(i),'b',b(i)+(nbp-1)*delta(i),...
            'e',b(i)+(nep-1)*delta(i),'npts',size(data(i).dep,1),...
            'depmen',mean(data(i).dep(:)),...
            'depmin',min(data(i).dep(:)),...
            'depmax',max(data(i).dep(:)));
    end
    % single point fix
    if(size(data(i).dep,1)==1); data(i)=ch(data(i),'delta',nan); end
end

% remove failed records
if(trim); data(failed)=[]; end

end
