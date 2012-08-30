function [data,failed]=cut(data,varargin)
%CUT    Cut a window out of SEIZMO records
%
%    Usage:    data=cut(data,pdw)
%              data=cut(...,'cmplist',list,...)
%              data=cut(...,'fill',logical,...)
%              data=cut(...,'filler',value,...)
%              data=cut(...,'iztype',iztype,...)
%              [data,failed]=cut(...,'trim',logical,...)
%
%    Description:
%     DATA=CUT(DATA,PDW) cuts the records in DATA to the window limits
%     defined by PDW and outputs the updated dataset.  Header fields are
%     updated to match the windowed data.  Works with unevenly sampled
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
%        Absolute times are assumed to be in UTC.  Note that DATEVEC
%        interpretation does not properly handle times within a leap
%        second.
%
%     OFFSET is a numeric value giving the offset to be added to REF to 
%     define the position of the window start/stop with respect to the 
%     reference.  The offset can be a __COLUMN__ vector of values (one
%     per record) to define different offsets for each record (Note that
%     REF cannot be a vector of reference positions!).
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
%     DATA=CUT(...,'CMPLIST',LIST,...) allows specifying only certain
%     components to be included in the output (for multicomponent files).
%     Basically any components not in the list are cut.  LIST should be
%     either a row vector of indices or ':'.  It may also be an array of
%     rows of indices or ':' (use a cell array to give a mixture.  Default
%     is ':' (all components).
%
%     DATA=CUT(...,'FILL',TRUE|FALSE,...) turns on/off the filling of data
%     gaps to allow records that don't extend for the entire window to do
%     so by padding them with zeros.  Does not work with unevenly sampled
%     data.  This option is useful for operations that require records to
%     have the same length.  See option 'FILLER' to alter the fill value. 
%     By default 'FILL' is false (no fill).
%
%     DATA=CUT(...,'FILLER',VALUE,...) changes the fill value to VALUE.
%     Adjusting this value does NOT turn option 'FILL' to true.  By default
%     'FILLER' is set to 0 (zero).
%
%     DATA=CUT(...,'IZTYPE',IZTYPE,...) sets the header field 'iztype' of
%     the output records to IZTYPE.    This value is passed directly to
%     CHANGEHEADER as 'changeheader(DATA,'iztype',IZTYPE)'.  The default
%     IZTYPE is [] and changes nothing.
%
%     [DATA,FAILED]=CUT(...,'TRIM',TRUE|FALSE) turns on/off deleting
%     records that have no data (after cutting) as well as spectral and xyz
%     records (which are not supported by CUT).  Optional output FAILED 
%     gives the indices of these records (with respect to their position in
%     the input DATA).  By default 'TRIM' is set to true (deletes records).
%
%    Notes:
%     - Windowing of spectral and xyz records is not supported.  By default
%       they are deleted from DATA (see option 'TRIM' to change this 
%       behavior).
%     - Windows with a start position after an end position will return
%       empty records (if 'TRIM' is set to false) with headers updated
%       accordingly rather than returning an error.
%     - Multiple component data is supported.
%     - FILL only works with evenly sampled data.
%     
%    Header changes: B, E, NPTS, DELTA, NCMP, DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Cut a 400 sample window starting with the 33rd sample:
%     data=cut(data,'x',33,'n',400);
%
%     % Cut out a 90s window around t1:
%     data=cut(data,'t1',-30,'t1',60);
%
%     % Cut one hour starting at the origin time, 
%     % padding incomplete records with zeros:
%     data=cut(data,'o',0,3600,'fill',true);
%
%     % Cut records to first 300 seconds:
%     data=cut(data,'b',0,300);
%
%     % Cut from 300 to 500 seconds relative to reference time:
%     data=cut(data,300,500);
%
%     % Cut out data from the last week:
%     data=cut(data,datevec(now-7),datevec(now));
%
%    See also: READDATAWINDOW, SYNCHRONIZE, TIMESHIFT

%     Version History:
%        Oct. 30, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Jan. 28, 2008 - new sachp support
%        Feb.  6, 2008 - renamed to cutim
%        Feb. 23, 2008 - works with GLGC, GENUMDESC
%        Feb. 25, 2008 - bugfix
%        Feb. 28, 2008 - minor improvements
%        Feb. 29, 2008 - better leven check
%        Mar.  2, 2008 - dataless support
%        Mar.  4, 2008 - doc update
%        Apr. 17, 2008 - PDW now like SAC
%        Apr. 18, 2008 - bugfix
%        May  12, 2008 - uses new dep* formula
%        June 12, 2008 - doc update
%        June 23, 2008 - major doc update
%        June 30, 2008 - fixed dataless support, .dep & .ind rather than .x
%                        & .t, improved checks
%        Sep. 15, 2008 - minor doc update
%        Sep. 22, 2008 - error msg fixes
%        Nov. 24, 2008 - fixed some bugs for uneven, fixed fill bug, one
%                        changeheader call, support new cutparameters,
%                        better checking, doc update
%        Mar. 12, 2009 - fixed 3 more fill bugs :(
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Jan. 28, 2010 - proper SEIZMO handling, seizmoverbose support,
%                        better warning messages
%        Feb.  2, 2010 - versioninfo caching
%        Mar.  8, 2010 - dropped versioninfo caching
%        Feb. 11, 2011 - mass nargchk fix
%        Nov.  2, 2011 - doc update, allow absolute time input
%        Dec.  1, 2011 - IZTYPE option added
%        Feb.  1, 2012 - better getheader usage
%        Mar.  8, 2012 - minor doc update (abs times are UTC)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  8, 2012 at 15:05 GMT

% todo:

% input check
error(nargchk(1,inf,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% make sure global settings are kept
try
    % check headers
    data=checkheader(data);

    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % parse cut parameters
    option=cutparameters(nrecs,varargin{:});

    % header info
    [b,delta,e,npts,ncmp,iftype,leven]=getheader(data,...
        'b','delta','e','npts','ncmp','iftype id','leven lgc');
    
    % find spectral/xyz
    amph=strcmpi(iftype,'iamph');
    rlim=strcmpi(iftype,'irlim');
    xyz=strcmpi(iftype,'ixyz');
    uneven=strcmpi(leven,'false');
    even=~uneven;
    
    % check iftype
    failed=false(nrecs,1);
    if(any(amph | rlim))
        failed(amph | rlim)=true;
        warning('seizmo:cut:illegalFiletype',...
            ['Record(s):\n' sprintf('%d ',find(amph | rlim)) ...
            '\nIllegal operation on spectral record(s)!']);
    elseif(any(xyz))
        failed(xyz)=true;
        warning('seizmo:cut:illegalFiletype',...
            ['Record(s):\n' sprintf('%d ',find(xyz)) ...
            '\nIllegal operation on xyz record(s)!']);
    end
    
    % check uneven & fill
    if(any(uneven & option.FILL))
        warning('seizmo:cut:noFillUneven',...
            ['Record(s):\n' sprintf('%d ',find(uneven & option.FILL)) ...
            '\nCannot fill unevenly sampled record(s)!']);
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
                    error('seizmo:cut:malformedAbsTimeStr',...
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
                    error('seizmo:cut:malformedAbsTimeStr',...
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

    % window start point/time
    if(strcmpi(option.REF1,'z'))
        bt=option.OFFSET1;
        bp=round((bt-b)./delta)+1;
    elseif(strcmpi(option.REF1,'x'))
        bp=round(option.OFFSET1);
    else
        bt=getheader(data,option.REF1)+option.OFFSET1;
        bp=round((bt-b)./delta)+1;
    end

    % window end point/time
    if(strcmpi(option.REF2,'z'))
        et=option.OFFSET2;
        ep=round((et-b)./delta)+1;
    elseif(strcmpi(option.REF2,'x'))
        ep=round(option.OFFSET2);
    elseif(strcmpi(option.REF2,'n'))
        ep=bp+round(option.OFFSET2)-1;
    else
        et=getheader(data,option.REF2)+option.OFFSET2;
        ep=round((et-b)./delta)+1;
    end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
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
    disp('Cutting Record(s)');
    print_time_left(0,nrecs);
end

% loop through each file
[depmen,depmin,depmax]=deal(nan(nrecs,1));
for i=find(~failed')
    % evenly spaced
    if(even(i))
        % cut
        if(nnp(i))
            data(i).dep=data(i).dep(nbp(i):nep(i),option.CMPLIST{i});
        else
            data(i).dep=data(i).dep([],option.CMPLIST{i});
        end
        
        % get ncmp
        ncmp(i)=size(data(i).dep,2);
        
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
    else % unevenly spaced
        % get begin point for window
        if(~strcmpi(ref1,'x'))
            % find first point after beginning of window
            temp=find(data(i).ind>=bt(i),1);
            
            % handle no points after begin
            if(isempty(temp))
                % empty window
                bp(i)=npts(i)+1;
            % 'round' to closer point
            elseif(temp>1)
                % index of point on other side of begin time
                temp2=temp-1;
                
                % figure out which is closer
                if(abs(bt(i)-data(i).ind(temp))...
                        <=abs(bt(i)-data(i).ind(temp2)))
                    bp(i)=temp;
                else
                    bp(i)=temp2;
                end
            % 1st point is closest
            else
                bp(i)=temp;
            end
        end
        
        % get end point for window
        if(~strcmpi(ref2,'n') && ~strcmpi(ref2,'x'))
            % find last point before end of window
            temp=find(data(i).ind<=et(i),1,'last');
            
            % handle no points before end
            if(isempty(temp))
                % empty window
                ep(i)=0;
            % 'round' to closest point
            elseif(temp<npts(i))
                % index of point on other side of end time
                temp2=temp+1;
                
                % figure out which is closer
                if(abs(et(i)-data(i).ind(temp))...
                        <=abs(et(i)-data(i).ind(temp2)))
                    ep(i)=temp;
                else
                    ep(i)=temp2;
                end
            % last point is closest
            else
                ep(i)=temp;
            end
        end
        
        % boundary conditions
        nbp(i)=max([bp(i) 1]);
        nep(i)=min([ep(i) npts(i)]);
        npts(i)=max(nep(i)-nbp(i)+1,0);
        
        % cut
        data(i).dep=data(i).dep(nbp(i):nep(i),option.CMPLIST{i});
        data(i).ind=data(i).ind(nbp(i):nep(i),1);
        
        % get ncmp
        ncmp(i)=size(data(i).dep,2);
        
        % handle empty
        if(~npts(i))
            b(i)=nan; e(i)=nan; 
            failed(i)=true;
            
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end
        
        % update timing/dep*
        b(i)=data(i).ind(1);
        e(i)=data(i).ind(end);
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        if(npts(i)>1)
            delta(i)=(e(i)-b(i))/(npts(i)-1);
        end
    end
    
    % detail message
    if(verbose); print_time_left(i,nrecs); end
end

% detail message
if(verbose && (isempty(i) || i~=nrecs))
    print_time_left(nrecs,nrecs);
end

% get new b/e/npts for evenly spaced
fe=option.FILL & even;
nfe=~option.FILL & even;
% to fill, even spaced
if(any(fe))
    e(fe)=b(fe)+(ep(fe)-1).*delta(fe);
    b(fe)=b(fe)+(bp(fe)-1).*delta(fe);
    npts(fe)=max(ep(fe)-bp(fe)+1,0);
% not to fill, even spaced
elseif(any(nfe))
    e(nfe)=b(nfe)+(nep(nfe)-1).*delta(nfe);
    b(nfe)=b(nfe)+(nbp(nfe)-1).*delta(nfe);
    npts(nfe)=nnp(nfe);
end

% try update headers
try
    % toggle checking off
    seizmocheck_state(false);
    
    % update headers
    data=changeheader(data,'b',b,'e',e,'delta',delta,'npts',npts,...
        'ncmp',ncmp,'depmen',depmen,'depmin',depmin,'depmax',depmax,...
        'iztype',option.IZTYPE);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% removed failed/empty cut records
if(option.TRIM); data(failed)=[]; end

end
