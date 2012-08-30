function [data,failed]=readdatawindow(data,varargin)
%READDATAWINDOW    Read data window into SEIZMO data structure
%
%    Usage:    data=readdatawindow(data,pdw)
%              data=readdatawindow(...,'cmplist',list,...)
%              data=readdatawindow(...,'fill',logical,...)
%              data=readdatawindow(...,'filler',value,...)
%              data=readdatawindow(...,'iztype',iztype,...)
%              [data,failed]=readdatawindow(...,'trim',logical,...)
%
%    Description:
%     READDATAWINDOW(DATA,PDW) reads a data window PDW from SEIZMO
%     compatible datafiles into the SEIZMO structure DATA, returning the
%     updated structure.  This provides a mechanism similar to the SAC
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
%     REF is a string that refers to a reference value for the independent
%     component (usually time).  It can be any valid numeric header field, 
%     'z', 'x','n', or an absolute time (allowed formats defined below). 
%     Where:
%      - 'z' is the zero position for the record == 0.
%      - 'x' indicates that the following offset is a sample number (first 
%        sample is 1).  Note that this defaults to sample 0 when given
%        without an offset.
%      - 'n' indicates that the following offset is the length of the
%        window in number of samples.  Note that this defaults to length 0
%        when given without an offset.  Also note that this is only valid
%        as a second REF.
%      - Absolute times can be a __ROW__ vector of [year dayofyear],
%        [year month dayofmonth], [year dayofyear hour minute seconds], or
%        [year month dayofmonth hour minute seconds]. Absolute times may
%        also be a string as given by the KZDTTM or KZDATE fields or a
%        string that can be interpreted by the Matlab function DATEVEC.
%        Note that DATEVEC interpretation does not properly handle times
%        within a leap second.
%
%     OFFSET is a numeric value giving the offset to be added to REF to 
%     define the position of the window start/stop with respect to the 
%     reference.  The offset can be a __COLUMN__ vector of values, one for
%     each record to be read, to define different offsets for each record
%     (note that REF can not be a vector of reference positions).
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
%     DATA=READDATAWINDOW(...,'CMPLIST',LIST,...) allows reading in only
%     specific components of records.  Should be either a row vector of
%     indices or ':'.  May also be an array of rows of indices or ':' (use
%     a cell array to get a mixture.  Default is ':' (all components).
%
%     DATA=READDATAWINDOW(...,'FILL',TRUE|FALSE,...) turns on/off the
%     filling of data gaps to allow records that don't extend for the
%     entire window to do so by padding them with zeros.  This is useful
%     for operations that require records to all have the same length.  See
%     option 'FILLER' to alter the fill value.  By default 'FILL' is false
%     (no fill).
%
%     DATA=READDATAWINDOW(...,'FILLER',VALUE,...) changes the fill value to
%     VALUE.  Adjusting this value does NOT turn option 'FILL' to true.  By
%     default 'FILLER' is set to 0 (zero).
%
%     DATA=READDATAWINDOW(...,'IZTYPE',IZTYPE,...) sets the header field
%     'iztype' of the output records to IZTYPE.    This value is passed
%     directly to CHANGEHEADER as 'changeheader(DATA,'iztype',IZTYPE)'.
%     The default IZTYPE is [] and changes nothing.
%
%     [DATA,FAILED]=READDATAWINDOW(...,'TRIM',TRUE|FALSE,...) turns on/off
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
%    Header changes: B, E, NPTS, DELTA, NCMP, DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % ALL THE FOLLOWING EXAMPLES REQUIRE YOU TO 
%     % REMEMBER TO READ IN THE HEADER INFO FIRST:
%     data=readheader('*.SAC')
%
%     % read in only the first 300 samples of records:
%     data=readdatawindow(data,'x',1,300)
%     
%     % read in a 90 second window around the t1 arrival, 
%     % padding data gaps with zeros as necessary:
%     data=readdatawindow(data,'t1',-30,60,'fill',true)
%
%     % read in only the 123rd sample:
%     data=readdatawindow(data,'x',123,123)
%     % or
%     data=readdatawindow(data,'x',123,'n',1)
%
%     % Read in data from only the last week:
%     data=readdatawindow(data,datevec(now-7),datevec(now));
%
%    See also: CUT, READHEADER, READDATA, READSEIZMO, WRITEHEADER,
%              WRITESEIZMO, BSEIZMO, GETHEADER, CHANGEHEADER, LISTHEADER

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
%        Mar. 12, 2009 - fixed 3 more fill bugs :(
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Feb.  3, 2010 - update for XDIR fixes, better SAC bug workaround,
%                        proper SEIZMO handling, seizmoverbose support
%        Feb. 11, 2011 - mass nargchk fix
%        Nov.  2, 2011 - doc update, allow absolute time input
%        Dec.  1, 2011 - IZTYPE option added
%        Jan. 30, 2012 - better getheader usage, seizmosize override for
%                        speed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 30, 2012 at 15:05 GMT

% todo:

% input check
error(nargchk(1,11,nargin));

% headers setup (also checks struct)
[h,vi]=versioninfo(data);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt reading in data window
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % parse cut parameters
    option=cutparameters(nrecs,varargin{:});

    % header info
    [b,delta,e,npts,ncmp,iftype,leven]=getheader(data,...
        'b','delta','e','npts','ncmp','iftype id','leven lgc');
    
    % estimated filesize from header
    [est_bytes,hbytes,dbytes]=seizmosize(h,vi,npts,ncmp,iftype,leven);

    % find spectral/xyz
    amph=strcmpi(iftype,'iamph');
    rlim=strcmpi(iftype,'irlim');
    xyz=strcmpi(iftype,'ixyz');
    
    % find (un)even
    uneven=strcmpi(leven,'false');
    even=~uneven;

    % allocate bad records matrix
    failed=false(nrecs,1);
    if(any(amph | rlim))
        failed(amph | rlim)=true;
        warning('seizmo:readdatawindow:illegalFiletype',...
            ['Record(s):\n' sprintf('%d ',find(amph | rlim)) ...
            '\nIllegal operation on spectral record(s)!']);
    elseif(any(xyz))
        failed(xyz)=true;
        warning('seizmo:readdatawindow:illegalFiletype',...
            ['Record(s):\n' sprintf('%d ',find(xyz)) ...
            '\nIllegal operation on xyz record(s)!']);
    end
    
    % debugging
    %option
    %[option.OFFSET1 option.OFFSET2]
    
    % check for abs time strings
    % - require it is 8+ characters
    % - allow kzdate, kzdttm, & datevec recongnised formats
    if(numel(option.REF1>7) && ischar(option.REF1))
        if(numel(option.REF1)==29 && all(option.REF1([12 16])=='()'))
            option.REF1=str2double({option.REF1(1:4) option.REF1(13:15) ...
                option.REF1(18:19) option.REF1(21:22) option.REF1(24:29)});
        elseif(numel(option.REF1)==16 && all(option.REF1([12 16])=='()'))
            option.REF1=str2double({option.REF1(1:4) option.REF1(13:15)});
        else
            try
                option.REF1=datevec(option.REF1);
            catch
                if(numel(option.REF1)>8)
                    error('seizmo:readdatawindow:malformedAbsTimeStr',...
                        'Malformed abs time string: ''%s''',option.REF1);
                end
            end
        end
    end
    if(numel(option.REF2>7) && ischar(option.REF2))
        if(numel(option.REF2)==29 && all(option.REF2([12 16])=='()'))
            option.REF2=str2double({option.REF2(1:4) option.REF2(13:15) ...
                option.REF2(18:19) option.REF2(21:22) option.REF2(24:29)});
        elseif(numel(option.REF2)==16 && all(option.REF2([12 16])=='()'))
            option.REF2=str2double({option.REF2(1:4) option.REF2(13:15)});
        else
            try
                option.REF2=datevec(option.REF2);
            catch
                if(numel(option.REF2)>8)
                    error('seizmo:readdatawindow:malformedAbsTimeStr',...
                        'Malformed abs time string: ''%s''',option.REF2);
                end
            end
        end
    end
    
    % debugging
    %option
    
    % convert abs time to relative time from each record's reference time
    if(isnumeric(option.REF1))
        option.OFFSET1=option.OFFSET1 ...
            +timediff(cell2mat(getheader(data,'z')),option.REF1,'UTC');
        option.REF1='z';
    end
    if(isnumeric(option.REF2))
        option.OFFSET2=option.OFFSET2 ...
            +timediff(cell2mat(getheader(data,'z')),option.REF2,'UTC');
        option.REF2='z';
    end
    
    % debugging
    %option
    %[option.OFFSET1 option.OFFSET2]

    % let readdata/cut handle unevenly sampled records (minus file deletion)
    if(any(uneven))
        % detail message
        if(verbose)
            disp('Reading Unevenly Sampled Record(s) & Cutting Window');
        end
        
        % read in unevenly sampled datafiles
        [data(uneven),failed(uneven)]=readdata(data(uneven),'trim',false);
        if(any(uneven & ~failed))
            % cut unevenly sampled datafiles
            uneven=(uneven & ~failed);
            if(isempty(option.IZTYPE)); iztype=option.IZTYPE;
            else iztype=option.IZTYPE(uneven);
            end
            [data(uneven),failed(uneven)]=cut(data(uneven),...
                option.REF1,option.OFFSET1(uneven),...
                option.REF2,option.OFFSET2(uneven),...
                'fill',option.FILL,'filler',option.FILLER,...
                'cmplist',option.CMPLIST,'trim',false,...
                'iztype',iztype);
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
    
    % detail message
    if(verbose)
        disp('Reading Data Window of Evenly Sampled Record(s)');
        print_time_left(0,nrecs);
    end

    % loop through each evenly spaced file
    [depmen,depmin,depmax]=deal(nan(nrecs,1));
    neven=sum(even & ~failed);
    for i=find(even & ~failed).'
        % construct fullname
        name=[data(i).path data(i).name];

        % open file
        fid=fopen(name,'r',data(i).byteorder);

        % check that it opened
        if(fid<0)
            warning('seizmo:readdatawindow:badFID',...
                'Record: %d, File not openable: %s !',i,name);
            failed(i)=true;
            % detail message
            if(verbose); print_time_left(i,neven,true); end
            continue;
        end

        % file size
        fseek(fid,0,'eof');
        bytes=ftell(fid);
        
        % byte size check
        if(bytes>est_bytes(i))
            if(bytes==(hbytes(i)+2*dbytes(i)))
                % SAC BUG: converting a spectral record to a time series
                %      record and then writing it out, writes out a 2nd
                %      component that contains gibberish.
            else
                % size big enough but inconsistent - skip
                warning('seizmo:readdatawindow:badFileSize',...
                    ['Record: %d\nFile: %s\n'...
                    'Filesize does not match header info!\n'...
                    '%d (estimated) > %d (on disk) --> Skipping!'],...
                    i,name,est_bytes(i),bytes);
                failed(i)=true;
                % detail message
                if(verbose); print_time_left(i,neven,true); end
                continue;
            end
        elseif(bytes<est_bytes(i))
            % size too small - skip
            fclose(fid);
            warning('seizmo:readdatawindow:badFileSize',...
                ['Record: %d, File: %s'...
                'Filesize of file does not match header info!\n'...
                '%d (estimated) < %d (on disk) --> Skipping!'],...
                i,name,est_bytes(i),bytes);
            failed(i)=true;
            % detail message
            if(verbose); print_time_left(i,neven,true); end
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
                data(i).dep(:,j)=...
                    fread(fid,nnp(i),['*' h(vi(i)).data.store]);
            end
        end

        % data read in
        data(i).hasdata=true;

        % close file
        fclose(fid);

        % add filler
        if(option.FILL(i))
            data(i).dep=...
                [ones(min(0,ep(i))-bp(i)+1,ncmp(i))*option.FILLER(i);...
                data(i).dep;...
                ones(ep(i)-max(npts(i),bp(i)-1),ncmp(i))*option.FILLER(i)];
        end

        % update dep*
        if(numel(data(i).dep))
            depmen(i)=nanmean(data(i).dep(:));
            depmin(i)=min(data(i).dep(:));
            depmax(i)=max(data(i).dep(:));
        else
            failed(i)=true;
        end
        
        % detail message
        if(verbose); print_time_left(i,neven); end
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

    % update headers
    if(isempty(option.IZTYPE)); iztype=option.IZTYPE;
    else iztype=option.IZTYPE(even);
    end
    data(even)=changeheader(data(even),'iztype',iztype,...
        'b',b(even),'e',e(even),'npts',npts(even),'ncmp',ncmp(even),...
        'depmen',depmen(even),'depmin',depmin(even),'depmax',depmax(even));

    % remove failed records
    if(option.TRIM); data(failed)=[]; end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end


function [bytes,hbytes,dbytes]=seizmosize(h,vi,npts,ncmp,iftype,leven)
% override version (for speed)

% trim down header info
h=[h.data];

% filetype
count=strcmpi(iftype,'itime')...
    +strcmpi(iftype,'ixy')+strcmpi(iftype,'ixyz')...
    +2*strcmpi(iftype,'irlim')+2*strcmpi(iftype,'iamph');

% multi-component
count=ncmp.*count;

% uneven sampling
count=count+strcmpi(leven,'false');

% final tally
hbytes=[h(vi).startbyte].';
dbytes=count.*npts.*[h(vi).bytesize].';
bytes=hbytes+dbytes;

end
