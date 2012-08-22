function [data,ssidx]=syncedsets(data,idx1,idx2,idx3,varargin)
%SYNCEDSETS    Output synchronized sets of SEIZMO records
%
%    Usage:    [data,ssidx]=syncedsets(data,idx1,idx2,idx3)
%              [...]=syncedsets(...,'minoverlap',value,...)
%              [...]=syncedsets(...,'overlapunits',units,...)
%              [...]=syncedsets(...,'adjust',method,...)
%              [...]=syncedsets(...,'shiftmax',value,...)
%              [...]=syncedsets(...,'shiftunits',units,...)
%              [...]=syncedsets(...,'interpolate',method,...)
%              [...]=syncedsets(...,'useabsolutetiming',logical,...)
%              [...]=syncedsets(...,'timing',standard,...)
%              [...]=syncedsets(...,'allocate',size,...)
%              [...]=syncedsets(...,'verbose',logical,...)
%              [...]=syncedsets(...,'debug',logical,...)
%
%    Description:
%     [DATA,SSIDX]=SYNCEDSETS(DATA,IDX1,IDX2,IDX3) takes records in SEIZMO
%     dataset DATA and indexing from FINDTRIPLETS or HORZPAIRS and creates
%     a new dataset of synchronized record sets where the input record sets
%     overlapped.  Sets are defined by the given indexing and usually
%     correspond to a distinct seismic station.  This is useful for
%     rotation, particle motion visualization, etc.  Sets must overlap by
%     at least 2 samples (see options MINOVERLAP and OVERLAPUNITS to alter
%     this) and overlap is checked based on their absolute timing (see
%     options USEABSOLUTETIMING and TIMING options to alter this).  Synced
%     sets can be easily indexed using the SSIDX array (equal size to the
%     output DATA).  For instance, SSIDX elements equal to 1 correspond to
%     records in the first synced set and so on.  Ordering of the
%     components in the synced sets follows that of the input indexing.
%
%     [...]=SYNCEDSETS(...,'MINOVERLAP',VALUE,...) sets the minimum time
%     overlap to consider.  By default MINOVERLAP is 2 intervals (the units
%     to may be changed to seconds using the OVERLAPUNITS option below).
%
%     [...]=SYNCEDSETS(...,'OVERLAPUNITS',UNITS,...) changes the units of
%     the MINOVERLAP option.  By default UNITS is 'INTERVALS'.  This can
%     be changed to 'SECONDS' if that is more useful.
%
%     [...]=SYNCEDSETS(...,'ADJUST',METHOD,...) allows changing which
%     evenly-sampled records of a set is shifted/interpolated to time-align
%     with the other.  There are six choices: 'FIRST' 'LAST' 'LONGER'
%     'SHORTER' 'ONE' & 'TWO'.  The default is 'SHORTER' (which adjusts the
%     shorter records to time-align with the longer).
%
%     [...]=SYNCEDSETS(...,'SHIFTMAX',VALUE,...) allows changing the cap on
%     when the record-to-be-adjusted (see ADJUST option) is interpolated or
%     shifted to time-align with the other record.  By default the
%     MAXVALUE is 0.01 intervals -- meaning the adjusted record is only
%     shifted if the time change is very minor, otherwise the data will be
%     interpolated at the aligned sample times.  A MAXVALUE of 0.5
%     intervals will always shift the data to the new times without
%     interpolating new values.  Really the choice depends on how sensitive
%     you think your data is to time shifts and/or how much you trust the
%     timing of the adjusted record.  If you're trying to get relative
%     arrival times of P recordings then you probably are worried about
%     minor shifts in timing (which begs the question of why you are
%     dealing with such crappy data in the first place).  The default is
%     basically saying the data timing is pretty accurate and any new data
%     from new time points should be interpolated unless the time shift is
%     damn small and really would not change anything.
%
%     [...]=SYNCEDSETS(...,'SHIFTUNITS',UNITS,...) changes the units of the
%     SHIFTMAX option.  By default UNITS is 'INTERVALS'.  This can be
%     changed to 'SECONDS' if that is more useful.
%
%     [...]=SYNCEDSETS(...,'INTERPOLATE',METHOD,...) allows changing the
%     interpolation method.  The choices are basically those allowed in
%     Matlab's INTERP1 command: 'spline' 'pchip' 'linear' and 'nearest'.
%     The default is 'spline', which is continuous in the 1st and 2nd
%     derivatives.  Look out for artifacting if you use one of the other
%     options and are going to differentiate the data later.
%
%     [...]=SYNCEDSETS(...,'USEABSOLUTETIMING',LOGICAL,...) allows turning
%     on/off the usage of the reference time fields to figure out the
%     timing of data.  This can be safely turned off if all the data share
%     the same reference time.  Leave it on if your reference times vary
%     with each record.
%
%     [...]=SYNCEDSETS(...,'TIMING',STANDARD,...) changes the timing
%     standard assumed for the reference time.  The choices are: 'UTC' and
%     'TAI'.  The default is 'UTC', which has leap second awareness.  This
%     is useful for dealing with data that have had UTC leap seconds
%     properly inserted.  Proper handling of leap seconds requires that the
%     records' have their reference time at the actual UTC time.  If the
%     recording equipment doesn't actually handle leap seconds then some
%     time adjustment is/was needed for the data.  See LEAPSECONDS for more
%     info.  The 'TAI' option is useful for data without any leap second
%     concerns.
%
%     [...]=SYNCEDSETS(...,'ALLOCATE',SIZE,...) sets the temporary space
%     initially allocated for rotating records.  This is just a guess of
%     the maximum number of synced sets records created.  The default value
%     is equal to the number of input records.  Not really worth changing.
%
%     [...]=SYNCEDSETS(...,'VERBOSE',LOGICAL,...) turns on/off the progress
%     bar.  Useful for seeing how far along we are & how long there is to
%     finish.  Default is TRUE (on).
%
%     [...]=SYNCEDSETS(...,'DEBUG',LOGICAL,...) turns on/off detailed
%     debugging messages.  Default is FALSE (off).
%
%    Notes:
%     - Run FIXDELTA first to take care of small differences in sample
%       rates caused by floating point inaccuracies!
%
%    Header changes: B, E, NPTS, DEPMEN, DEPMIN, DEPMAX,
%                    (see CHECKHEADER for more)
%
%    Examples:
%     % 
%
%    See also: FINDTRIPLETS, MAKETRIPLETS, UNMAKETRIPLETS, ROTATE3,
%              PLOTPM3, PMMOVIE3, HORZPAIRS

%     Version History:
%        Aug.  1, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  1, 2012 at 23:18 GMT

% todo:

% check nargin
if(nargin<4 || mod(nargin,2))
    error('seizmo:syncedsets:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data,'NONTIME_IFTYPE','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% attempt to make synced sets
try
    % number of records
    nrecs=numel(data);
    
    % check index inputs
    if(~isnumeric(idx1) || ~isnumeric(idx2) || ~isnumeric(idx3))
        error('seizmo:syncedsets:badInput',...
            'Index matrices must be numeric!');
    elseif(~isequal(numel(idx1),numel(idx2),numel(idx3)))
        error('seizmo:syncedsets:badInput',...
            'Index matrices are not the same size!');
    elseif(numel(idx1)==0)
        error('seizmo:syncedsets:badInput',...
            'Index matrices empty! No sets to sync!');
    elseif(~all(idx1==fix(idx1) & idx2==fix(idx2) & idx3==fix(idx3)))
        error('seizmo:syncedsets:badInput',...
            'Index matrices must be positive indices!');
    elseif(any([max(idx1) max(idx2) max(idx3)]>nrecs))
        error('seizmo:syncedsets:badInput',...
            'Index matrices have too high of values!');
    elseif(any([min(idx1) min(idx2) min(idx3)]<1))
        error('seizmo:syncedsets:badInput',...
            'Index matrices must be positive indices!');
    end
    
    % parse parameters
    option=parse_syncedsets_parameters(varargin{:});
    
    % get header fields
    if(option.USEABSOLUTETIMING)
        [b,e,delta,npts,nz,leven]=getheader(data,...
            'b','e','delta','npts','nz','leven lgc');
    else
        [b,e,delta,npts,leven]=getheader(data,...
            'b','e','delta','npts','leven lgc');
        nz=nan(nrecs,6);
    end
    
    % get filenames (for debug output)
    if(option.DEBUG)
        name=strcat({data.path}.',{data.name}.');
    end
    
    % get start and end of records in absolute time
    if(option.USEABSOLUTETIMING)
        if(strcmp(option.TIMING,'utc'))
            ab=gregorian2modserial(utc2tai(...
                [nz(:,1:4) nz(:,5)+nz(:,6)/1000+b]));
            ae=gregorian2modserial(utc2tai(...
                [nz(:,1:4) nz(:,5)+nz(:,6)/1000+e]));
        else
            ab=gregorian2modserial([nz(:,1:4) nz(:,5)+nz(:,6)/1000+b]);
            ae=gregorian2modserial([nz(:,1:4) nz(:,5)+nz(:,6)/1000+e]);
        end
    else
        ab=[zeros(nrecs,1) b];
        ae=[zeros(nrecs,1) e];
    end
    
    % number of sets
    nsets=max(idx2);
    
    % allocate output dataset (setting as empty)
    ndata=data(nan(0,0));
    
    % allocate header field arrays
    nnz=nan(0,6); nab=nan(0,2); nae=nab;
    nnpts=nan(0,1); ndepmin=nnpts; ndepmen=nnpts; ndepmax=nnpts;
    
    % detail message
    if(option.VERBOSE)
        disp('Making Synchronized Record Set(s)');
        print_time_left(0,nsets);
    end
    
    % loop over pairs
    for i=1:nsets
        % separate indices based on component
        sidx=idx1(idx2==i); ns=numel(sidx);
        cmpmax=max(idx3(idx2==i));
        cidx=cell(cmpmax,1);
        for j=1:cmpmax
            cidx{j}=idx1(idx2==i & idx3==j);
        end
        
        % detail message
        if(option.DEBUG)
            fprintf('\n\nProcessing Set: %d\n',i);
            disp(['Members: ' sprintf('%d ',sidx)]);
            fprintf('Number in Set: %d\n',ns);
        end
        
        % adjust absolute times of set to be based on same day
        day1=ab(sidx(1),1); % may not be the best choice but it is easy
        cb=cell(cmpmax,1); ce=cb;
        for j=1:cmpmax
            cb{j}=ab(cidx{j},2)+86400*(ab(cidx{j},1)-day1);
            ce{j}=ae(cidx{j},2)+86400*(ae(cidx{j},1)-day1);
        end
        
        % split action based on leven/delta
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % unevenly spaced or delta unmatched
        if(any(strcmpi(leven(sidx),'false')) ...
                || numel(unique(delta(sidx)))~=1)
            % first make sure independent component exists
            data(sidx)=makeuneven(data(sidx));
            
            % get independent component
            t=cell(cmpmax,1);
            for j=1:cmpmax
                t{j}=cell(nc(j),1);
                for k=1:nc(j)
                    t{j}{k}=data(cidx{j}(k)).ind(:)...
                    -data(cidx{j}(k)).ind(1)+cb{j}(k);
                end
            end
            
            % loop over all possible pairings
            [ndata,ninfo,ssidx]=findsets(...
                data,name,delta,cidx,ab,ae,true,option);
            
            for j=1:nc(1)
                [oi,ob,oe,ol,ot]=findset(cb{1}(j),ce{1}(j),t{1}{j},...
                    cb,ce,t,true,2,option);
                if(ol(1)>0)
                    % detail message
                    if(option.DEBUG)
                        % lots of info
                        % form strings for each record
                        tmpstr=cell(2,cmpmax);
                        for k=1:cmpmax
                            tmpi=cidx{k}(oi(k));
                            tmpstr{2,k}=sprintf(...
                                [' %d - %s\n'...
                                ' begin: Day: %d Second: %f\n'...
                                ' end:   Day: %d Second: %f\n'...
                                ' avg. delta: %f seconds\n'...
                                ' npts:  %d\n'],tmpi,name{tmpi},...
                                ab(tmpi,1),ab(tmpi,2),...
                                ae(tmpi,1),ae(tmpi,2),...
                                delta(tmpi),npts(tmpi));
                        end
                        tmpstr(1,:)=[{'\nSyncing Records:\n'} ...
                            repmat({'and\n'},1,maxcmp-1)];
                        fprintf([tmpstr{:} ...
                            'Overlap: %f seconds (%d samples)\n\n'],...
                            ol(1),ol(2));
                    end
                    die
                    % increment counters
                    c1=c1+2; c2=c2+2;
                    
                    % copy records
                    ndata([c1 c2],1)=data([cidx1(j) cidx2(k)]);
                    
                    % interpolate values
                    ndata(c1).dep=interp1(...
                        t1{j},data(cidx1(j)).dep,t,option.INTERPOLATE);
                    ndata(c2).dep=interp1(...
                        t2{k},data(cidx2(k)).dep,t,option.INTERPOLATE);
                    
                    % get updated header values
                    nnpts([c1 c2])=numel(ndata(c1).dep);
                    ndelta([c1 c2])=(t(end)-t(1))/nnpts(c1);
                    nnz([c1 c2],:)=nz([cidx1(j) cidx2(k)],:);
                    nleven([c1 c2])={'false'};
                    
                    % get updated times
                    nab([c1 c2],1)=day1+floor(t(1)/86400);
                    nae([c1 c2],1)=day1+floor(t(end)/86400);
                    nab([c1 c2],2)=mod(t(1),86400);
                    nae([c1 c2],2)=mod(t(end),86400);
                    
                    % dep*
                    if(nnpts(c1))
                        ndepmin(c1)=min(ndata(c1).dep(:));
                        ndepmen(c1)=mean(ndata(c1).dep(:));
                        ndepmax(c1)=max(ndata(c1).dep(:));
                        ndepmin(c2)=min(ndata(c2).dep(:));
                        ndepmen(c2)=mean(ndata(c2).dep(:));
                        ndepmax(c2)=max(ndata(c2).dep(:));
                    end
                    
                    % detail message
                    if(option.DEBUG)
                        % output reversal flag
                        if(rf2); tmp='REVERSE';
                        else tmp='NORMAL';
                        end
                        
                        % lots of info
                        fprintf(...
                            ['\nOutput Records:\n'...
                            ' %d - %s\n'...
                            ' azimuth: %f degrees\n'...
                            ' begin: Day: %d Second: %f\n'...
                            ' end:   Day: %d Second: %f\n'...
                            ' avg. delta: %f seconds\n'...
                            ' npts:  %d\n'...
                            'and\n'...
                            ' %d - %s\n'...
                            ' azimuth: %f degrees (%s)\n'...
                            ' begin: Day: %d Second: %f\n'...
                            ' end:   Day: %d Second: %f\n'...
                            ' avg. delta: %f seconds\n'...
                            ' npts:  %d\n\n'],...
                            c1,name{cidx1(j)},pto1,...
                            nab(c1,1),nab(c1,2),...
                            nae(c1,1),nae(c1,2),...
                            ndelta(c1),nnpts(c1),...
                            c2,name{cidx2(k)},pto2,tmp,...
                            nab(c2,1),nab(c2,2),...
                            nae(c2,1),nae(c2,2),...
                            ndelta(c2),nnpts(c2));
                    end
                end
            end
        else % evenly spaced, delta matched
            % get independent component
            t1=cell(n1,1); t2=cell(n2,1);
            for j=1:n1
                t1{j}=cb1(j)+(0:npts(cidx1(j))-1).'*delta(pidx(1));
            end
            for j=1:n2
                t2{j}=cb2(j)+(0:npts(cidx2(j))-1).'*delta(pidx(1));
            end
            
            % loop over all possible pairings
            for j=1:n1
                for k=1:n2
                    % get overlap info
                    [ob,oe,ol]=get_overlap(cb1(j),ce1(j),cb2(k),ce2(k),...
                        delta(pidx(1)),false);
                    
                    % skip if unacceptable overlap
                    switch lower(option.OVERLAPUNITS)
                        case 'intervals'
                            if(ol(2)<option.MINOVERLAP); continue; end
                        case 'seconds'
                            if(ol(1)<option.MINOVERLAP); continue; end
                    end
                    
                    % detail message
                    if(option.DEBUG)
                        % input reversal flag
                        if(rf1); tmp='REVERSE';
                        else tmp='NORMAL';
                        end
                        
                        % lots of info
                        fprintf(...
                            ['\nRotating Record:\n'...
                            ' %d - %s\n'...
                            ' azimuth: %f degrees\n'...
                            ' begin: Day: %d Second: %f\n'...
                            ' end:   Day: %d Second: %f\n'...
                            ' delta: %f seconds\n'...
                            ' npts:  %d\n'...
                            'with\n'...
                            ' %d - %s\n'...
                            ' azimuth: %f degrees (%s)\n'...
                            ' begin: Day: %d Second: %f\n'...
                            ' end:   Day: %d Second: %f\n'...
                            ' delta: %f seconds\n'...
                            ' npts:  %d\n'...
                            'Overlap: %f seconds (%f intervals)\n\n'],...
                            cidx1(j),name{cidx1(j)},cmpaz(cidx1(j)),...
                            ab(cidx1(j),1),ab(cidx1(j),2),...
                            ae(cidx1(j),1),ae(cidx1(j),2),...
                            delta(cidx1(j)),npts(cidx1(j)),...
                            cidx2(k),name{cidx2(k)},cmpaz(cidx2(k)),tmp,...
                            ab(cidx2(k),1),ab(cidx2(k),2),...
                            ae(cidx2(k),1),ae(cidx2(k),2),...
                            delta(cidx2(k)),npts(cidx2(k)),ol(1),ol(2));
                    end
                    
                    % increment counters
                    c1=c1+2; c2=c2+2;
                    
                    % copy records
                    ndata([c1 c2],1)=data([cidx1(j) cidx2(k)]);
                    
                    % syncedsets
                    [t,ndata(c1).dep,ndata(c2).dep]=syncedsets_even(...
                        t1{j},data(cidx1(j)).dep,cmpaz(cidx1(j)),pto1,...
                        t2{k},data(cidx2(k)).dep,rf1,rf2,ob,oe,...
                        delta(pidx(1)),option);
                    
                    % get updated header values
                    nnz([c1 c2],:)=nz([cidx1(j) cidx2(k)],:);
                    nab([c1 c2],1)=day1; nae([c1 c2],1)=day1;
                    nab([c1 c2],2)=t(1); nae([c1 c2],2)=t(end);
                    ndelta([c1 c2])=delta(pidx(1));
                    nleven([c1 c2])=leven(pidx(1));
                    
                    % dep*
                    nnpts([c1 c2])=numel(ndata(c1).dep);
                    if(nnpts(c1))
                        ndepmin(c1)=min(ndata(c1).dep(:));
                        ndepmen(c1)=mean(ndata(c1).dep(:));
                        ndepmax(c1)=max(ndata(c1).dep(:));
                        ndepmin(c2)=min(ndata(c2).dep(:));
                        ndepmen(c2)=mean(ndata(c2).dep(:));
                        ndepmax(c2)=max(ndata(c2).dep(:));
                    end
                    
                    % detail message
                    if(option.DEBUG)
                        % output reversal flag
                        if(rf2); tmp='REVERSE';
                        else tmp='NORMAL';
                        end
                        
                        % lots of info
                        fsprintf(...
                            ['\nOutput Records:\n'...
                            ' %d - %s\n'...
                            ' azimuth: %f degrees\n'...
                            ' begin: Day: %d Second: %f\n'...
                            ' end:   Day: %d Second: %f\n'...
                            ' delta: %f seconds\n'...
                            ' npts:  %d\n'...
                            'and\n'...
                            ' %d - %s\n'...
                            ' azimuth: %f degrees (%s)\n'...
                            ' begin: Day: %d Second: %f\n'...
                            ' end:   Day: %d Second: %f\n'...
                            ' delta: %f seconds\n'...
                            ' npts:  %d\n\n'],...
                            c1,name{cidx1(j)},pto1,...
                            nab(c1,1),nab(c1,2),...
                            nae(c1,1),nae(c1,2),...
                            ndelta(c1),nnpts(c1),...
                            c2,name{cidx2(k)},pto2,tmp,...
                            nab(c2,1),nab(c2,2),...
                            nae(c2,1),nae(c2,2),...
                            ndelta(c2),nnpts(c2));
                    end
                end
            end
        end
        
        % detail message
        if(option.DEBUG)
            fprintf('\nFinished Rotating Group: %d',i);
            disp(['Members: ' sprintf('%d ',oldc2+1:c2)]);
            fprintf('Number in Group: %d\n',c2-oldc2);
        end
        
        % detail message
        if(option.VERBOSE); print_time_left(i,npairs); end
    end
    
    % get relative times from absolute
    if(option.USEABSOLUTETIMING)
        if(strcmp(option.TIMING,'utc'))
            naz=gregorian2modserial(utc2tai(...
                [nnz(:,1:4) nnz(:,5)+nnz(:,6)/1000]));
            nb=(nab(:,1)-naz(:,1))*86400+(nab(:,2)-naz(:,2));
            ne=(nae(:,1)-naz(:,1))*86400+(nae(:,2)-naz(:,2));
        else
            naz=gregorian2modserial(...
                [nnz(:,1:4) nnz(:,5)+nnz(:,6)/1000]);
            nb=(nab(:,1)-naz(:,1))*86400+(nab(:,2)-naz(:,2));
            ne=(nae(:,1)-naz(:,1))*86400+(nae(:,2)-naz(:,2));
        end
    else
        nb=nab(:,2);
        ne=nae(:,2);
    end
    
    % error if no output records
    if(~numel(ndata))
        error('seizmo:syncedsets:noRotatedRecords',...
            'No rotatable horizontal pairs found!');
    end
    
    % update header
    data=changeheader(ndata,'b',nb,'e',ne,'npts',nnpts,...
        'kcmpnm',nkcmpnm,'cmpaz',ncmpaz,'leven',nleven,'delta',ndelta,...
        'depmin',ndepmin,'depmax',ndepmax,'depmen',ndepmen);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end


function [option]=parse_syncedsets_parameters(varargin)
% parses/checks syncedsets p/v pairs

% valid values for string options
valid.OVERLAPUNITS={'seconds' 'intervals'};
valid.INTERPOLATE={'spline' 'pchip' 'linear' 'nearest'};
valid.ADJUST={'longer' 'shorter' 'first' 'last' 'one' 'two'};
valid.SHIFTUNITS={'seconds' 'intervals'};
valid.TIMING={'utc' 'tai'};

% defaults
option.MINOVERLAP=2; % minimum overlap of pair to syncedsets (in samples)
option.OVERLAPUNITS='intervals'; % seconds/intervals
option.INTERPOLATE='spline'; % spline/pchip/linear/nearest
option.ADJUST='shorter'; % longer/shorter/first/last
option.SHIFTMAX=0.01; % interval: 0-0.5 , seconds: 0+
option.SHIFTUNITS='intervals'; % seconds/intervals
option.TIMING='utc'; % utc/tai
option.USEABSOLUTETIMING=true; % true/false
option.ALLOCATE=0; % size of tmp space
option.VERBOSE=seizmoverbose; % default to seizmoverbose state
option.DEBUG=seizmodebug; % default to seizmodebug state

% require options specified by strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:syncedsets:badInput',...
        'Not all options are specified with a string!');
end

% parse options
for i=1:2:nargin
    switch lower(varargin{i})
        case {'minoverlap' 'minover' 'minlap' 'minov' 'mino' 'mo'}
            option.MINOVERLAP=varargin{i+1};
        case {'overlapunits' 'ounits' 'ounit' 'ou'}
            option.OVERLAPUNITS=varargin{i+1};
        case {'interpolate' 'interp' 'int' 'i'}
            option.INTERPOLATE=varargin{i+1};
        case {'adjust' 'adj' 'a'}
            option.ADJUST=varargin{i+1};
        case {'shiftunits' 'shiftu' 'sunit' 'su'}
            option.SHIFTUNITS=varargin{i+1};
        case {'shiftmax' 'shift' 'smax' 'max' 'sh'}
            option.SHIFTMAX=varargin{i+1};
        case {'timing' 'time'}
            option.TIMING=varargin{i+1};
        case {'useabsolutetiming' 'useabs' 'abs' 'abstiming'}
            option.USEABSOLUTETIMING=varargin{i+1};
        case {'allocate' 'all'}
            option.ALLOCATE=varargin{i+1};
        case {'verbose' 'v'}
            option.VERBOSE=varargin{i+1};
        case {'debug' 'd'}
            option.DEBUG=varargin{i+1};
        otherwise
            error('seizmo:syncedsets:badInput',...
                'Unknown option: %s',varargin{i});
    end
end

% check options
fields=fieldnames(option);
for i=1:numel(fields)
    % get value of field and do a basic check
    value=option.(fields{i});
    
    % specific checks
    switch lower(fields{i})
        case {'shiftmax' 'minoverlap'}
            if(~isnumeric(value) || ~isscalar(value))
                error('seizmo:syncedsets:badInput',...
                    '%s must be a scalar real number!',fields{i});
            end
        case {'overlapunits' 'interpolate' 'adjust' 'shiftunits' 'timing'}
            if(~ischar(value) || size(value,1)~=1 ...
                    || ~any(strcmpi(value,valid.(fields{i}))))
                error('seizmo:syncedsets:badInput',...
                    ['%s option must be one of the following:\n'...
                    sprintf('%s ',valid.(fields{i}){:})],fields{i});
            end
        case 'allocate'
            if(~isnumeric(value) || fix(value)~=value)
                error('seizmo:syncedsets:badInput',...
                    'ALLOCATE must be a scalar integer!');
            end
        case {'useabsolutetiming' 'verbose' 'debug'}
            if(~islogical(value) || ~isscalar(value))
                error('seizmo:syncedsets:badInput',...
                    '%s option must be a logical!',fields{i});
            end
        otherwise
            error('seizmo:syncedsets:badInput',...
                'Unknown option: %s !',fields{i});
    end
end

% turn off verbose if debugging
if(option.DEBUG); option.VERBOSE=false; end

end


function [ok,oi,ob,oe,ol,ot]=findset(cb1,ce1,t1,cb,ce,t,uneven,k,option)
% loop over this level looking for overlap
% if found then ...
%   check if at final level
%     if so then return values
%   otherwise recurse at next level
% otherwise return 0 overlap

% k=level
if(k>numel(t)); error('whoops'); end
for i=1:numel(t{k})
    if(uneven)
        % get overlapping samples in uneven
        [ob,oe,ol,ot]=get_overlap(cb1,ce1,...
            cb{k}(i),ce{k}(i),[t1; t{k}{i}],uneven);
        % skip if unacceptable overlap
        switch lower(option.OVERLAPUNITS)
            case 'intervals' % actually samples here
                if(ol(2)<option.MINOVERLAP); continue; end
            case 'seconds'
                if(ol(1)<option.MINOVERLAP); continue; end
        end
        % so overlap is good...final level or not?
        if(k==numel(t))
            % make set
        else % recurse
            [ok,oi,ob,oe,ol,ot]=findset(ob,oe,ot,...
                cb,ce,t,uneven,k+1,option);
        end
    else
        
    end
end

end


function [ob,oe,ol,varargout]=get_overlap(b1,e1,b2,e2,delta,uneven)
% get overlap range/length

if(uneven) % unevenly sampled
    % what samples are within range
    t=unique(delta(delta>=b1 & delta<=e1 & delta>=b2 & delta<=e2));
    ob=t(1);
    oe=t(end);
    ol=[t(end)-t(1) numel(t)];
    varargout{1}=t;
else % evenly sampled
    ob=max(b1,b2);
    oe=min(e1,e2);
    ol=[oe-ob (oe-ob)/delta];
end

end


function [data,ab,ae]=make_syncedsets_uneven(...
    data,idx,name,ab,ae,day1,delta,ot,t,option)

% debug message
% records in
% index, name, start/end, delta, overlap

% interpolate values
for i=1:numel(data)
    data(i).dep=interp1(t{i},data(i).dep,ot,option.INTERPOLATE);
    data(i).ind=ot; % FIXME!?!?!
end

% get updated header values
nnpts([c1 c2])=numel(ndata(c1).dep);
ndelta([c1 c2])=(t(end)-t(1))/nnpts(c1);
nnz([c1 c2],:)=nz([cidx1(j) cidx2(k)],:);
nleven([c1 c2])={'false'};

% get updated times
nab([c1 c2],1)=day1+floor(t(1)/86400);
nae([c1 c2],1)=day1+floor(t(end)/86400);
nab([c1 c2],2)=mod(t(1),86400);
nae([c1 c2],2)=mod(t(end),86400);

% dep*
if(nnpts(c1))
    ndepmin(c1)=min(ndata(c1).dep(:));
    ndepmen(c1)=mean(ndata(c1).dep(:));
    ndepmax(c1)=max(ndata(c1).dep(:));
    ndepmin(c2)=min(ndata(c2).dep(:));
    ndepmen(c2)=mean(ndata(c2).dep(:));
    ndepmax(c2)=max(ndata(c2).dep(:));
end
end


function [t,x1,x2]=make_syncedsets_even(t1,x1,t2,x2,ob,oe,delta,option)

% points basically in overlap
bidx1=find(abs(t1-ob)<=delta/2,1,'first');
bidx2=find(abs(t2-ob)<=delta/2,1,'first');
npts=round((oe-ob)/delta)+1;

% what is the necessary shift
shift=mod(abs(t1(1)-t2(1)),delta);
if(shift>delta/2); shift=delta-shift; end
if(option.DEBUG)
    fprintf('Misalignment: %f seconds (%f intervals)\n',...
        shift,shift/delta);
end

% shiftmax
switch lower(option.SHIFTUNITS)
    case 'seconds'
        shiftmax=option.SHIFTMAX;
    case 'intervals'
        shiftmax=option.SHIFTMAX*delta;
end

% which do we adjust
switch lower(option.ADJUST)
    case 'first'
        if(t1(1)==t2(1))
            move1=t1(end)<t2(end);
        else
            move1=t1(1)<t2(1);
        end
    case 'last'
        if(t1(1)==t2(1))
            move1=t1(end)>t2(end);
        else
            move1=t1(1)>t2(1);
        end
    case 'longer'
        move1=numel(t1)>numel(t2);
    case 'shorter'
        move1=numel(t1)<numel(t2);
    case 'one'
        move1=true;
    otherwise % two
        move1=false;
end

% get new times
if(move1)
    t=t2(bidx2+(0:npts-1));
else % move 2
    t=t1(bidx1+(0:npts-1));
end

% get x1,x2
if(shift<=shiftmax)
    % ok we are just shifting one record to time align
    if(option.DEBUG); disp('Adjust Method: SHIFT'); end
    x1=x1(bidx1+(0:npts-1));
    x2=x2(bidx2+(0:npts-1));
else % interpolate
    % interpolating one record to time align
    if(option.DEBUG); disp('Adjust Method: INTERPOLATE'); end
    if(move1)
        x1=interp1(t1,x1,t,option.INTERPOLATE,'extrap');
        x2=x2(bidx2+(0:npts-1));
    else % interpolate 2
        x1=x1(bidx1+(0:npts-1));
        x2=interp1(t2,x2,t,option.INTERPOLATE,'extrap');
    end
end

end

