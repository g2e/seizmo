function [data,failed]=rpdw(data,varargin)
%RPDW    Reads partial data window of datafiles into SEIZMO data structure
%
%    Description: RPDW(DATA,PDW) reads a partial data window PDW from 
%     SEIZMO compatible datafiles into the SEIZMO structure DATA, returning
%     the updated structure.  This provides a mechanism similar to the SAC 
%     'cut' command to limit memory/cpu usage related to reading in large 
%     datafiles where only a segment is needed.  Note that the datafile 
%     headers must exist/be read into DATA before running RPDW (use RH).
%     Header fields will be updated to match the windowed data.
%     
%     PDW is a set of several arguments that reference header fields in 
%     DATA to define the data window.  PDW must be one of the following 
%     forms:
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
%     REF is a string that refers to a reference field/value for the 
%     independent component (usually time).  It can be any valid numeric 
%     header field or 'z', 'x', or 'n'.
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
%     REF can not be a vector of reference positions).
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
%     records that are unsupported (spectral and xyz files), that failed to
%     be read, or would have no data in the window.  Optional output FAILED
%     returns a logical array (equal in size to DATA) with elements set to 
%     TRUE for records that encountered problems.  By default 'TRIM' is set
%     to TRUE (deletes records).
%
%    Notes:
%     - Partial reads of spectral and xyz files are not supported.  They
%       will not be read in with RPDW.  By default they are deleted from
%       DATA (see option 'TRIM' to change this behavior).
%     - Windows with a start position after an end position will return
%       empty records (with headers updated accordingly) rather than
%       returning an error.  These records will be deleted unless 'trim' is
%       set to false.
%     - Multiple component data is supported.
%     - Partial reads of unevenly sampled data is not supported.  They are 
%       passed to RDATA and CUTIM instead.
%     - FILL only works with evenly sampled data.
%     - Records with a negative DELTA will return as empty.
%
%    System requirements: Matlab 7
%
%    Header changes: B, E, NPTS, DELTA, ODELTA, LEVEN
%                    DEPMEN, DEPMIN, DEPMAX
%
%    Usage:    data=rpdw(data,pdw)
%              data=rpdw(data,pdw,'fill',true|false)
%              data=rpdw(data,pdw,'fill',true|false,'filler',value)
%              [data,failed]=rpdw(data,pdw,...'trim',true|false)
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
%        Mar.  4, 2008 - doc update
%        Apr. 17, 2008 - PDW now like SAC
%        Apr. 18, 2008 - bugfix
%        May  12, 2008 - uses new dep* formula
%        June 12, 2008 - doc update
%        June 23, 2008 - major doc update, bugfix for uneven
%                        files, some other cleanups for readability
%        June 30, 2008 - fixed dataless support, .dep & .ind rather than .x
%                        & .t, fix for single point data
%        Sep. 22, 2008 - minor doc update, error msg fixes
%        Sep. 27, 2008 - doc update, LEVEN for dataless & 1pnt, updated for
%                        GET_N_CHECK & VINFO, checks for datafile size,
%                        updated RDATA call, single CH call
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 28, 2008 at 03:35 GMT

% todo:
% - dataless, 0pt & 1pt handling
%   - don't change leven or delta
%   - how to read/cut under such conditions
% - cmp option to subselect components to read
%   - how to specify individually

% input check
error(nargchk(1,11,nargin))

% check data structure
error(seischk(data))

% number of records
nrecs=numel(data);

% parse cut parameters
option=cutparameters(nrecs,varargin{:});

% header info
ncmp=gncmp(data);
[b,delta,e,depmen,depmin,depmax,npts]=...
    gh(data,'b','delta','e','depmen','depmin','depmax','npts');
iftype=genumdesc(data,'iftype');
leven=glgc(data,'leven');
error(lgcchk('leven',leven(npts>1)))

% index by spacing
even=strcmpi(leven,'true');
uneven=strcmpi(leven,'false');
eveni=find(even).';

% look out for undefined leven
undef=(~even & ~uneven); 

% estimated filesize from header
est_bytes=seissize(data);

% headers setup
[h,vi]=vinfo(data);

% window start point/value
bt=[]; bp=[];
if(strcmpi(option.REF1,'z'))
    bt=option.OFFSET1;
elseif(strcmpi(option.REF1,'x'))
    bp=round(option.OFFSET1);
else
    bt=gh(data,option.REF1)+option.OFFSET1;
end

% window stop point/value
et=[]; ep=[];
if(strcmpi(option.REF2,'z'))
    et=option.OFFSET2;
elseif(strcmpi(option.REF2,'x') || strcmpi(option.REF2,'n'))
    ep=round(option.OFFSET2);
else
    et=gh(data,option.REF2)+option.OFFSET2;
end

% check for nans
if(any(isnan(bt)) || any(isnan(bp)) || any(isnan(et)) || any(isnan(ep)))
    error('seizmo:rpdw:nanHeaderField','Window position is NaN!')
end

% allocate bad records matrix
failed=false(nrecs,1);

% let rdata/cutim handle unevenly sampled records (minus file deletion)
if(any(uneven | undef))
    % read in unevenly sampled datafiles
    uneven=(uneven | undef);
    [data(uneven),failed(uneven)]=rdata(data(uneven),'trim',false);
    if(any(uneven & ~failed))
        % cut unevenly sampled datafiles
        uneven=(uneven & ~failed);
        [data(uneven),failed(uneven)]=cutim(data(uneven),...
            option.REF1,option.OFFSET1(uneven),option.REF2,option.OFFSET2(uneven),...
            'fill',option.FILL,'filler',option.FILLER,'trim',false);
    end
end

% loop through each evenly spaced file
for i=eveni
    % skip dataless
    if(npts(i)==0)
        failed(i)=true; 
        warning('seizmo:rpdw:noData',...
            'Illegal operation on dataless file!');
        continue; 
    end
    
    % look out for negative delta
    if(delta(i)<=0)
        failed(i)=true;
        warning('seizmo:rpdw:negDelta',...
            'DELTA for record %d is zero/negative, deleting!',i);
        continue;
    end
    
    % check for unsupported filetypes
    if(strcmpi(iftype(i),'General XYZ (3-D) file'))
        failed(i)=true;
        warning('seizmo:rpdw:illegalFiletype',...
            'Illegal operation on xyz file!');
        continue;
    elseif(any(strcmpi(iftype(i),{'Spectral File-Real/Imag'...
            'Spectral File-Ampl/Phase'})))
        failed(i)=true;
        warning('seizmo:rpdw:illegalFiletype',...
            'Illegal operation on spectral file!');
        continue;
    end
    
    % closest points to begin and end
    if(~strcmpi(option.REF1,'x'))
        bp(i)=round((bt(i)-b(i))/delta(i))+1;
    end
    if(strcmpi(option.REF2,'n'))
        ep2=bp(i)+ep(i)-1;
    elseif(strcmpi(option.REF2,'x'))
        ep2=ep(i);
    else
        ep2=round((et(i)-b(i))/delta(i))+1;
    end
    
    % boundary conditions
    nbp=max([bp(i) 1]);
    nep=min([ep2 npts(i)]);
    nnp=max([nep-nbp+1 0]);
    
    % skip if new npts==0
    if(nnp==0)
        b(i)=nan; npts(i)=0; leven{i}='nan'; delta(i)=nan;
        e(i)=nan; depmen(i)=nan; depmin(i)=nan; depmax(i)=nan;
        failed(i)=true;
        continue;
    end
    
    % open file
    fid=fopen(data(i).name,'r',data(i).endian);
    
    % check that it opened
    if(fid<0)
        warning('seizmo:rpdw:badFID',...
            'File not openable: %s !',data(i).name);
        failed(i)=true;
        continue;
    end
    
    % file size
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    
    % byte size check
    if(bytes>est_bytes(i))
        % size big enough but inconsistent - read anyways (SAC bugfix)
        % SAC BUG: converting a spectral file to a time series file does
        % not deallocate the second component, thus the written file has
        % twice as much data.
        warning('seizmo:rpdw:badFileSize',...
            ['Filesize of file %s does not match header info!\n'...
            '%d (estimated) > %d (on disk) --> Reading Anyways!'...
            'This is usually caused by a SAC bug and can be ignored.'],...
            data(i).name,est_bytes(i),bytes);
    elseif(bytes<est_bytes(i))
        % size too small - skip
        fclose(fid);
        warning('seizmo:rpdw:badFileSize',...
            ['Filesize of file %s does not match header info!\n'...
            '%d (estimated) < %d (on disk) --> Skipping!'],...
            data(i).name,est_bytes(i),bytes);
        failed(i)=true;
        continue;
    end
    
    % preallocate data record with NaNs, deallocate timing
    data(i).dep=nan(nnp,ncmp(i),h(vi(i)).data.store);
    data(i).ind=[];
    
    % read in each component
    for j=1:ncmp(i)
        % move to first byte of window and read
        try
            fseek(fid,h(vi(i)).data.startbyte+...
                h(vi(i)).data.bytesize*((j-1)*npts(i)+nbp-1),'bof');
            temp=fread(fid,nnp,['*' h(vi(i)).data.store]);
            if(numel(temp)~=nnp)
                error('seizmo:rpdw:readFailed',...
                    'Read in of data failed: %s !',data(i).name);
            else
                data(i).dep(:,j)=temp;
                clear temp
            end
        catch
            warning('seizmo:rpdw:readFailed',...
                'Read in of data failed: %s !',data(i).name);
            failed(i)=true;
            clear temp
            break;
        end
    end
    
    % close file
    fclose(fid);
    
    % skip to next record if read failed
    if(failed(i)); continue; end
    
    % to fill
    if(option.FILL(i))
        % add filler
        data(i).dep=[ones(1-bp(i),ncmp(i))*option.FILLER(i); ...
            data(i).dep; ...
            ones(ep2-npts(i),ncmp(i))*option.FILLER(i)];
        
        % fix header
        b(i)=b(i)+(bp(i)-1)*delta(i);
        e(i)=b(i)+(ep2-1)*delta(i);
    % not to fill
    else
        % fix header
        b(i)=b(i)+(nbp-1)*delta(i);
        e(i)=b(i)+(nep-1)*delta(i);
    end
    
    % more header fix
    npts(i)=size(data(i).dep,1);
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
    
    % 1pnt fix
    if(size(data(i).dep,1)==1)
        delta(i)=nan; leven{i}='nan';
    end
end

% update headers
data(even)=ch(data(even),'b',b(even),'e',e(even),'delta',delta(even),...
    'leven',leven(even),'npts',npts(even),...
    'depmen',depmen(even),'depmin',depmin(even),'depmax',depmax(even));

% remove failed records
if(option.TRIM); data(failed)=[]; end

end
