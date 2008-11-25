function [data,failed]=readdatawindow(data,varargin)
%READDATAWINDOW    Read data window into SEIZMO data structure
%
%    Description: READDATAWINDOW(DATA,PDW) reads a data window PDW from 
%     SEIZMO compatible datafiles into the SEIZMO structure DATA, returning
%     the updated structure.  This provides a mechanism similar to the SAC 
%     'cut' command to limit memory/cpu usage related to reading in large 
%     datafiles where only a segment is needed.  Note that the datafile 
%     headers must exist/be read into DATA before running READDATAWINDOW
%     (use READHEADER). Header fields will be updated to match the windowed
%     data.
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
%     READDATAWINDOW(...,'CMPLIST',LIST) allows reading in only specific
%     components of records.  Should be either a row vector of indices or
%     ':'.  May also be an array of rows of indices or ':' (use a cell
%     array to get a mixture.  Default is ':' (all components).
%
%     READDATAWINDOW(...,'FILL',TRUE|FALSE) turns on/off the filling of
%     data gaps to allow records that don't extend for the entire window to
%     do so by padding them with zeros.  This is useful for operations that
%     require records to all have the same length.  See option 'FILLER' to
%     alter the fill value.  By default 'FILL' is false (no fill).
%
%     READDATAWINDOW(...,'FILLER',VALUE) changes the fill value to VALUE.
%     Adjusting this value does NOT turn option 'FILL' to true.  By default
%     'FILLER' is set to 0 (zero).
%
%     [DATA,FAILED]=READDATAWINDOW(...,'TRIM',TRUE|FALSE) turns on/off
%     deleting records that are unsupported for windowing (spectral and xyz
%     files), that failed to be read, or would have no data in the window.
%     Optional output FAILED returns a logical array (equal in size to
%     DATA) with elements set to TRUE for records that encountered
%     problems.  By default 'TRIM' is set to TRUE (deletes records).
%
%    Notes:
%     - Partial reads of spectral and xyz files are not supported.  They
%       will not be read in.  By default they are deleted from DATA (see
%       option 'TRIM' to change this behavior).
%     - Windows with a start position after an end position will return
%       empty records (with headers updated accordingly) rather than
%       returning an error.  These records will be deleted unless 'trim' is
%       set to false.
%     - Multiple component data is supported.
%     - Partial reads of unevenly sampled data is not supported.  They are 
%       passed to READDATA and CUT instead.
%     - FILL only works with evenly sampled data.
%     - Records with a negative DELTA will return as empty.
%
%    Tested on: Matlab r2007b
%
%    Header changes: B, E, NPTS, DELTA, NCMP, DEPMEN, DEPMIN, DEPMAX
%
%    Usage:    data=readdatawindow(data,pdw)
%              data=readdatawindow(...,'cmplist',list,...)
%              data=readdatawindow(...,'fill',logical,...)
%              data=readdatawindow(...,'filler',value,...)
%              [data,failed]=readdatawindow(...,'trim',logical,...)
%
%    Examples:
%     ALL THE FOLLOWING EXAMPLES REQUIRE YOU TO 
%     REMEMBER TO READ IN THE HEADER INFO FIRST:
%      data=readheader('*.SAC')
%
%     read in only the first 300 samples of records:
%      data=readdatawindow(data,'x',1,300)
%     
%     read in a 90 second window around the t1 arrival, 
%     padding data gaps with zeros as necessary:
%      data=readdatawindow(data,'t1',-30,60,'fill',true)
%
%     read in only the 123rd sample:
%      data=readdatawindow(data,'x',123,123)
%     or
%      data=readdatawindow(data,'x',123,'n',1)
%
%    See also: cut, readheader, readdata, readseizmo, writeheader,
%              writeseizmo, bseizmo, getheader, changeheader, listheader

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
%        Nov. 18, 2008 - updated to new name schema (now READDATAWINDOW),
%                        fixed fill bug
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 18, 2008 at 05:35 GMT

% todo:

% input check
error(nargchk(1,11,nargin));

% headers setup (also checks struct)
[h,vi]=versioninfo(data);

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% estimated filesize from header
est_bytes=seizmosize(data);

% number of records
nrecs=numel(data);

% parse cut parameters
option=cutparameters(nrecs,varargin{:});

% header info
ncmp=getncmp(data);
[b,delta,e,npts]=getheader(data,'b','delta','e','npts');
iftype=getenumdesc(data,'iftype');
leven=getlgc(data,'leven');

% index by spacing
even=strcmpi(leven,'true');
uneven=~even;

% allocate bad records matrix
failed=false(nrecs,1);

% let rdata/cutim handle unevenly sampled records (minus file deletion)
if(any(uneven))
    % read in unevenly sampled datafiles
    [data(uneven),failed(uneven)]=readdata(data(uneven),'trim',false);
    if(any(uneven & ~failed))
        % cut unevenly sampled datafiles
        uneven=(uneven & ~failed);
        [data(uneven),failed(uneven)]=cut(data(uneven),...
            option.REF1,option.OFFSET1(uneven),...
            option.REF2,option.OFFSET2(uneven),...
            'fill',option.FILL,'filler',option.FILLER,...
            'cmplist',option.CMPLIST,'trim',false);
    end
end

% window start point
if(strcmpi(option.REF1,'z'))
    bp=round((option.OFFSET1-b)./delta)+1;
elseif(strcmpi(option.REF1,'x'))
    bp=round(option.OFFSET1);
else
    bp=round((getheader(data,option.REF1)+option.OFFSET1-b)./delta)+1;
end

% window end point
if(strcmpi(option.REF2,'z'))
    ep=round((option.OFFSET2-b)./delta)+1;
elseif(strcmpi(option.REF2,'x'))
    ep=round(option.OFFSET2);
elseif(strcmpi(option.REF2,'n'))
    ep=bp+round(option.OFFSET2)-1;
else
    ep=round((getheader(data,option.REF2)+option.OFFSET2-b)./delta)+1;
end

% boundary conditions
nbp=max(bp,1);
nep=min(ep,npts);
nnp=max(nep-nbp+1,0);
if(any(~nnp))
    nbp(~nnp)=nan;
    nep(~nnp)=nan;
end

% get new b/e
fill=option.FILL;
nofill=~fill;
% to fill
if(any(fill))
    e(fill)=b(fill)+(ep(fill)-1).*delta(fill);
    b(fill)=b(fill)+(bp(fill)-1).*delta(fill);
    npts(fill)=max(ep(fill)-bp(fill)+1,0);
% not to fill
elseif(any(nofill))
    e(nofill)=b(nofill)+(nep(nofill)-1).*delta(nofill);
    b(nofill)=b(nofill)+(nbp(nofill)-1).*delta(nofill);
    npts(nofill)=nnp(nofill);
end

% fail evenly spaced records with 0 pts
failed(npts==0 & even)=true;

% loop through each evenly spaced file
[depmen,depmin,depmax]=swap(nan(nrecs,1));
for i=find(even).'
    % check for unsupported filetypes
    if(strcmpi(iftype(i),'General XYZ (3-D) file'))
        failed(i)=true;
        warning('seizmo:readdatawindow:illegalFiletype',...
            'Illegal operation on xyz file!');
        continue;
    elseif(any(strcmpi(iftype(i),{'Spectral File-Real/Imag'...
            'Spectral File-Ampl/Phase'})))
        failed(i)=true;
        warning('seizmo:readdatawindow:illegalFiletype',...
            'Illegal operation on spectral file!');
        continue;
    end
    
    % construct fullname
    name=fullfile(data(i).location,data(i).name);
    
    % open file
    fid=fopen(name,'r',data(i).byteorder);
    
    % check that it opened
    if(fid<0)
        warning('seizmo:readdatawindow:badFID',...
            'Record: %d, File not openable: %s !',i,name);
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
        warning('seizmo:readdatawindow:badFileSize',...
            ['Record: %d, File: %s'...
            'Filesize does not match header info!\n'...
            '%d (estimated) > %d (on disk) --> Reading Anyways!'...
            'This is usually caused by a SAC bug and can be ignored.'],...
            i,name,est_bytes(i),bytes);
    elseif(bytes<est_bytes(i))
        % size too small - skip
        fclose(fid);
        warning('seizmo:readdatawindow:badFileSize',...
            ['Record: %d, File: %s'...
            'Filesize of file does not match header info!\n'...
            '%d (estimated) < %d (on disk) --> Skipping!'],...
            i,name,est_bytes(i),bytes);
        failed(i)=true;
        continue;
    end
    
    % deal with components
    cmp=1:ncmp(i);
    cmp=cmp(option.CMPLIST{i});
    cmp=cmp(:).';
    ncmp(i)=numel(cmp);
    
    % preallocate data record with NaNs, deallocate timing
    data(i).dep=nan(nnp(i),ncmp(i),h(vi(i)).data.store);
    if(isfield(data,'ind')); data(i).ind=[]; end
    
    % read in each component
    if(nnp(i))
        for j=1:ncmp(i)
            % move to first byte of window and read
            fseek(fid,h(vi(i)).data.startbyte+h(vi(i)).data.bytesize...
                *((cmp(j)-1)*npts(i)+nbp(i)-1),'bof');
            data(i).dep(:,j)=fread(fid,nnp(i),['*' h(vi(i)).data.store]);
        end
    end
    
    % data read in
    data(i).hasdata=true;
    
    % close file
    fclose(fid);
    
    % add filler
    if(option.FILL(i))
        data(i).dep=...
            [ones(min(1,ep(i))-bp(i)+1,ncmp(i))*option.FILLER(i);...
            data(i).dep;...
            ones(ep(i)-max(npts(i),bp(i))+1,ncmp(i))*option.FILLER(i)];
    end
    
    % more header fix
    if(npts(i)>0)
        depmen(i)=mean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
    end
end

% update headers
warning('off','seizmo:changeheader:fieldInvalid')
data(even)=changeheader(data(even),'b',b(even),'e',e(even),...
    'npts',npts(even),'ncmp',ncmp(even),...
    'depmen',depmen(even),'depmin',depmin(even),'depmax',depmax(even));
warning('off','seizmo:changeheader:fieldInvalid')

% remove failed records
if(option.TRIM); data(failed)=[]; end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
