function [data]=rotate(data,varargin)
%ROTATE    Rotates SEIZMO records that are horizontal pairs
%
%    Usage:    data=rotate(data)
%              data=rotate(...,'to',azimuth,...)
%              data=rotate(...,'reverse',logical,...)
%              data=rotate(...,'minoverlap',value,...)
%              data=rotate(...,'overlapunits',units,...)
%              data=rotate(...,'kcmpnm1',string,...)
%              data=rotate(...,'kcmpnm2',string,...)
%              data=rotate(...,'adjust',method,...)
%              data=rotate(...,'shiftmax',value,...)
%              data=rotate(...,'shiftunits',units,...)
%              data=rotate(...,'interpolate',method,...)
%              data=rotate(...,'useabsolutetiming',logical,...)
%              data=rotate(...,'timing',standard,...)
%              data=rotate(...,'requiredcharfields',fields,...)
%              data=rotate(...,'requiredrealfields',fields,...)
%              data=rotate(...,'allocate',size,...)
%              data=rotate(...,'verbose',logical,...)
%              data=rotate(...,'debug',logical,...)
%
%    Description:
%     DATA=ROTATE(DATA) rotates the orthogonal horizontal pairs of records
%     in SEIZMO struct DATA to the great circle path (ie to the radial and
%     transverse directions).  See HORZPAIRS for details on how records are
%     determined to be horizontal pairs.  Pairs are returned such that
%     records 1 & 2 are a pair, 3 & 4 are a pair, etc.  The KCMPNM field is
%     altered to ??R for odd indexed records (indicating that these are the
%     radial component) and ??T for even indexed records (transverse
%     component) in output struct DATA.  Rotated pairs must overlap by at
%     least 2 samples (see options MINOVERLAP and OVERLAPUNITS to alter
%     this) and overlap is checked based on their absolute timing (see
%     options USEABSOLUTETIMING and TIMING options to alter this).
%     
%     DATA=ROTATE(...,'TO',AZIMUTH,...) sets the azimuth to rotate the
%     first component of a horizontal pair to.  AZIMUTH may be a header
%     field (only one field allowed!).  For instance, the default is
%     AZIMUTH set to 'GCP' (virtual field rotated 180deg from 'BAZ').  The
%     secondary component is output to AZIMUTH+90deg (see option REVERSE to
%     alter this).  Make sure AZIMUTH is consistent for a pair or ROTATE
%     will generate an error accordingly.
%
%     DATA=ROTATE(...,'REVERSE',LOGICAL,...) indicates if the secondary
%     component leads or follows the first component by 90deg.  The default
%     is false (leads by 90deg).  For example, if the output azimuth of the
%     1st component is 0 then the second will be 270 if REVERSE is TRUE.
%
%     DATA=ROTATE(...,'MINOVERLAP',VALUE,...) minimum time overlap to
%     consider two records for rotation.  By default MINOVERLAP is 2
%     intervals (the units to may be changed to seconds using the
%     OVERLAPUNITS option below).
%
%     DATA=ROTATE(...,'OVERLAPUNITS',UNITS,...) changes the units
%     of the MINOVERLAP option.  By default UNITS is 'INTERVALS'.  This can
%     be changed to 'SECONDS' if that is more useful.
%
%     DATA=ROTATE(...,'KCMPNM1',STRING,...)
%     DATA=ROTATE(...,'KCMPNM2',STRING,...)
%     changes the string of the output records' 'kcmpnm' header field.
%     Note that KCMPNM1 alters the kcmpnm of the 1st component (the one
%     oriented to the TO option's azimuth) while KCMPNM2 alters the 2nd
%     component.  KCMPNM gives the character(s) occupying the 3rd character
%     position and on.  So setting KCMPNM1 to 'N' would output odd-indexed
%     records with their 'kcmpnm' field set to ??N.  The first two
%     characters are preserved.  The default for KCMPNM1 is 'R' and KCMPNM2
%     is 'T'.
%
%     DATA=ROTATE(...,'ADJUST',METHOD,...) allows changing which record
%     out of a rotatible pair is shifted/interpolated to time-align with
%     the other.  There are six choices: 'FIRST' 'LAST' 'LONGER' 'SHORTER'
%     'ONE' & 'TWO'.  The default is 'SHORTER' (which adjusts the shorter
%     record to time-align with the longer).
%
%     DATA=ROTATE(...,'SHIFTMAX',VALUE,...) allows changing the cap on when
%     the record-to-be-adjusted (see ADJUST option) is interpolated or
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
%     DATA=ROTATE(...,'SHIFTUNITS',UNITS,...) allows changing the units
%     of the SHIFTMAX option.  By default UNITS is 'INTERVALS'.  This can
%     be changed to 'SECONDS' if that is more useful.
%
%     DATA=ROTATE(...,'INTERPOLATE',METHOD,...) allows changing the
%     interpolation method.  The choices are basically those allowed in
%     Matlab's INTERP1 command: 'spline' 'pchip' 'linear' and 'nearest'.
%     The default is 'spline', which is continuous in the 1st and 2nd
%     derivatives.  Look out for artifacting if you use one of the other
%     options and are going to differentiate the data later.
%
%     DATA=ROTATE(...,'USEABSOLUTETIMING',LOGICAL,...) allows turning
%     on/off the usage of the reference time fields to figure out the
%     timing of data.  This can be safely turned off if all the data share
%     the same reference time.  Leave it on if your reference times vary
%     with each record.
%
%     DATA=ROTATE(...,'TIMING',STANDARD,...) allows changing the timing
%     standard assumed for the reference time.  The choices are: 'UTC' and
%     'TAI'.  The default is 'UTC', which has leap second awareness.  This
%     is useful for dealing with data that have had UTC leap seconds
%     properly inserted (basically ROTATE won't even see the data overlap
%     because UTC times are converted to a leapless standard).  Proper
%     handling of leap seconds requires that the records' have their
%     reference time at the actual UTC time.  If the recording equipment
%     doesn't actually handle leap seconds then some time adjustment is/was
%     needed for the data.  See LEAPSECONDS for more info.  The 'TAI'
%     option is useful for data without any leap second concerns.
%
%     DATA=ROTATE(...,'REQUIREDCHARFIELDS',FIELDS,...) changes the
%     character fields required to be equal between records before they are
%     paired.  The list is a cellstring array.  The default is: {}.  Note
%     that pairing based on KNETWK, KSTNM, KHOLE, KCMPNM header fields is
%     hard-coded.
%
%     DATA=ROTATE(...,'REQUIREDREALFIELDS',FIELDS,...) changes the
%     numerical fields required to be equal between records before they are
%     paired.  The list must be a cellstring array.  The default is:
%     {'delta'}.  Note that CMPINC, LEVEN and NCMP are also required but
%     cannot be removed from the list.
%
%     DATA=ROTATE(...,'ALLOCATE',SIZE,...) sets the temporary space
%     initially allocated for rotating records.  This is just a guess of
%     the maximum number of rotated records created.  The default value is
%     equal to the number of input records.  Not really worth changing.
%
%     DATA=ROTATE(...,'VERBOSE',LOGICAL,...) turns on/off the rotate
%     progress bar.  Useful for seeing how far along we are & how long
%     there is to finish.  Default is TRUE (on).
%
%     DATA=ROTATE(...,'DEBUG',LOGICAL,...) turns on/off detailed debugging
%     messages.  Default is FALSE (off).
%
%    Notes:
%     - ROTATE requires several header fields to match in order to decide
%       which records should be rotated against one another.  At the bare
%       minimum the fields KNETWK, KSTNM, KHOLE, LEVEN, and NCMP must
%       match.  KCMPNM should match for the first 2 characters while CMPINC
%       must be 90.  See HORZPAIRS and its output to diagnose if your
%       records have sufficient header info.
%     - Run FIXDELTA first to take care of small differences in sample
%       rates caused by floating point inaccuracies!
%
%    Header changes: B, E, NPTS, DEPMEN, DEPMIN, DEPMAX, CMPAZ,
%                    DELTA, LEVEN (for unevenly sampled)
%                    (see CHECKHEADER for more)
%
%    Examples:
%     % Do the default (rotates to the great circle path):
%     rdata=rotate(data);
%
%     % Rotate to north (and east) + require rotated records are at least
%     % 100 seconds in length + correctly change kcmpnm of output records:
%     rdata=rotate(data,'to',0,'kcmpnm1','N','kcmpnm2','E',...
%                       'minoverlap',100,'overlapunits','seconds');
%
%     % Adjust some timing options to make rotation significantly faster
%     % (records must share the same reference time and we allow shifting
%     % the timing of records by up to delta/2):
%     rdata=rotate(data,'useabsolutetiming',false,'shiftmax',0.5);
%
%    See also: HORZPAIRS, ROTATE_CORRELATIONS

%     Version History:
%        Feb. 24, 2010 - initial version
%        Feb. 25, 2010 - significantly more debugging messages, some code
%                        cleaning
%        Aug. 21, 2010 - error usage fix, better checkheader usage, updated
%                        undef checks, minor fixes of debug messages
%        Sep. 29, 2010 - warn & skip rather than error on bad TO info
%        Oct.  6, 2010 - error if no output records
%        Dec. 13, 2011 - doc update
%        Jan. 28, 2012 - pass char to strnlen
%        Jan. 30, 2012 - drop SEIZMO global
%        Feb.  7, 2012 - fix no horizontals bug
%        Mar. 13, 2012 - use getheader improvements
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 23:18 GMT

% todo:
% - a 'norotate' would be nice
%   - ACTUALLY THIS IS JUST A 'THRU' OPTION SET TO 0
% - rotate 3cmp?

% check nargin
if(mod(nargin-1,2))
    error('seizmo:rotate:badNumInputs',...
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

% attempt rest
try
    % number of records
    nrecs=numel(data);
    
    % valid values for string options
    valid.OVERLAPUNITS={'seconds' 'intervals'};
    valid.INTERPOLATE={'spline' 'pchip' 'linear' 'nearest'};
    valid.ADJUST={'longer' 'shorter' 'first' 'last' 'one' 'two'};
    valid.SHIFTUNITS={'seconds' 'intervals'};
    valid.TIMING={'utc' 'tai'};
    
    % defaults
    option.TO='gcp'; % rotate to backazimuth
    option.MINOVERLAP=2; % minimum overlap of pair to rotate (in samples)
    option.OVERLAPUNITS='intervals'; % seconds/intervals
    option.KCMPNM1='R'; % string to change 1st kcmpnm(3:end) to
    option.KCMPNM2='T'; % string to change 2nd kcmpnm(3:end) to
    option.REQUIREDCHARFIELDS={}; % char fields to define pair
    option.REQUIREDREALFIELDS={'delta'}; % real fields to define pair
    option.REVERSE=false; % do not reverse orientation of 2nd cmp
    option.INTERPOLATE='spline'; % spline/pchip/linear/nearest
    option.ADJUST='shorter'; % longer/shorter/first/last
    option.SHIFTMAX=0.01; % interval: 0-0.5 , seconds: 0+
    option.SHIFTUNITS='intervals'; % seconds/intervals
    option.TIMING='utc'; % utc/tai
    option.USEABSOLUTETIMING=true; % true/false
    option.ALLOCATE=nrecs; % size of tmp space
    option.VERBOSE=seizmoverbose; % default to seizmoverbose state
    option.DEBUG=seizmodebug; % default to seizmodebug state
    
    % get options from command line
    for i=1:2:nargin-1
        if(~ischar(varargin{i}))
            error('seizmo:rotate:badInput',...
                'Options must be specified as a string!');
        end
        if(~isempty(varargin{i+1}))
            option.(upper(varargin{i}))=varargin{i+1};
        end
    end
    
    % check options
    fields=fieldnames(option);
    for i=1:numel(fields)
        % get value of field
        value=option.(fields{i});
        
        % specific checks
        switch lower(fields{i})
            case 'to'
                % handle header field case
                if(ischar(value))
                    value=getheader(data,value);
                end
                
                % no undefined allowed
                if(~isreal(value) || any(isnan(value) | isinf(value)))
                    error('seizmo:rotate:badInput',...
                        ['%s option must be a real or a ' ...
                        'real-valued header field!'],fields{i});
                end
                
                % expand scalar values
                if(isscalar(value))
                    value(1:nrecs,1)=value;
                end
                
                % force accuracy only to the thousandth
                % - this prevents fixed precision issues
                %to1=round(1000*value)/1000;
                to1=value; % more instructive of metadata problems
            case {'shiftmax' 'minoverlap'}
                if(~isreal(value) || ~isscalar(value))
                    error('seizmo:rotate:badInput',...
                        '%s must be a scalar real number!',fields{i});
                end
            case {'overlapunits' 'interpolate' ...
                    'adjust' 'shiftunits' 'timing'}
                if(~ischar(value) || size(value,1)~=1 ...
                        || ~any(strcmpi(value,valid.(fields{i}))))
                    error('seizmo:rotate:badInput',...
                        ['%s option must be one of the following:\n'...
                        sprintf('%s ',valid.(fields{i}){:})],fields{i});
                end
            case {'reverse' 'useabsolutetiming' 'verbose' 'debug'}
                if(~islogical(value) || ~isscalar(value))
                    error('seizmo:rotate:badInput',...
                        '%s option must be a logical!',fields{i});
                end
            case {'requiredcharfields' 'requiredrealfields'}
                % fix char arrays
                if(ischar(value))
                    value=cellstr(value);
                    option.(fields{i})=value;
                end
                if(~iscellstr(value))
                    error('seizmo:rotate:badInput',...
                        '%s option must be a cellstr of header fields!',...
                        fields{i});
                end
            case {'kcmpnm1' 'kcmpnm2'}
                if(~ischar(value) || size(value,1)~=1)
                    error('seizmo:rotate:badInput',...
                        '%s option must be a string!',fields{i});
                end
            case 'allocate'
                if(~isreal(value) || fix(value)~=value)
                    error('seizmo:rotate:badInput',...
                        'ALLOCATE must be a scalar integer!');
                end
            otherwise
                error('seizmo:rotate:badInput',...
                    'Unknown option: %s !',fields{i});
        end
    end
    
    % turn off verbose if debugging
    if(option.DEBUG); option.VERBOSE=false; end
    
    % get azimuth of secondary output component
    if(option.REVERSE)
        % second component trails by 90deg
        rf2=true;
        to2=mod(to1-90,360);
    else
        % second component leads by 90deg
        rf2=false;
        to2=mod(to1+90,360);
    end
    
    % get header fields
    if(option.USEABSOLUTETIMING)
        [b,e,delta,npts,nz,kcmpnm,cmpaz,leven]=getheader(data,...
            'b','e','delta','npts','nz','kcmpnm','cmpaz','leven lgc');
    else
        [b,e,delta,npts,kcmpnm,cmpaz,leven]=getheader(data,...
            'b','e','delta','npts','kcmpnm','cmpaz','leven lgc');
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
    
    % get horizontal pairs
    [idx1,idx2,idx3]=horzpairs(data,...
        'requiredcharfields',option.REQUIREDCHARFIELDS,...
        'requiredrealfields',option.REQUIREDREALFIELDS);
    npairs=max(idx2);
    
    % error if no output records
    if(isempty(npairs))
        error('seizmo:rotate:noRotatedRecords',...
            'No rotatable horizontal pairs found!');
    end
    
    % allocate output dataset (setting as empty)
    ndata=data(nan(option.ALLOCATE,0)); c1=-1; c2=0;
    
    % allocate header field arrays
    nnz=nan(option.ALLOCATE,6);
    nab=nan(option.ALLOCATE,2); nae=nab;
    nleven=cell(option.ALLOCATE,1); nkcmpnm=nleven;
    ndelta=nan(option.ALLOCATE,1); nnpts=ndelta;
    ndepmin=ndelta; ndepmen=ndelta; ndepmax=ndelta; ncmpaz=ndelta;
    
    % detail message
    if(option.VERBOSE)
        disp('Rotating Horizontal Record(s)');
        print_time_left(0,npairs);
    end
    
    % loop over pairs
    for i=1:npairs
        % separate indices based on component
        pidx=idx1(idx2==i); np=numel(pidx);
        cidx1=idx1(idx2==i & idx3==1); n1=numel(cidx1);
        cidx2=idx1(idx2==i & idx3==2); n2=numel(cidx2);
        
        % get reversal flag for input azimuths
        if(abs(mod(cmpaz(cidx1(1))-90,360)-cmpaz(cidx2(1)))<1)
            % second component trails by 90deg
            rf1=true;
        else
            % second component leads by 90deg
            rf1=false;
        end
        
        % check "to" azimuths
        if(numel(unique(to1(pidx)))>1)
            warning('seizmo:rotate:badInput',...
                ['Record(s):\n' sprintf('%d ',pidx) ...
                '\nFilename(s):\n' sprintf('%s\n',data(pidx).name) ...
                '\nAzimuth(s): ' sprintf('%f ',unique(to1(pidx))) ...
                '\nOnly 1 azimuth allowed to rotate to per pair!' ...
                '\n(For great circle path rotations this means the ' ...
                '\n station/earthquake location info is inconsistent!)']);
            if(option.VERBOSE); print_time_left(i,npairs,true); end
            continue;
        end
        pto1=unique(to1(pidx));
        pto2=unique(to2(pidx));
        
        % last new record (for detail messages)
        oldc2=c2;
        
        % detail message
        if(option.DEBUG)
            fprintf('\n\nProcessing Group: %d\n',i);
            disp(['Members: ' sprintf('%d ',pidx)]);
            fprintf('Number in Group: %d\n',np);
        end
        
        % adjust absolute times of pair to be based on same day
        day1=ab(pidx(1),1); % may not be the best choice but it is easy
        cb1=ab(cidx1,2)+86400*(ab(cidx1,1)-day1);
        ce1=ae(cidx1,2)+86400*(ae(cidx1,1)-day1);
        cb2=ab(cidx2,2)+86400*(ab(cidx2,1)-day1);
        ce2=ae(cidx2,2)+86400*(ae(cidx2,1)-day1);
        
        % rotate subfunction based on leven/delta
        % unevenly spaced or delta unmatched
        if(any(strcmpi(leven(pidx),'false')) ...
                || numel(unique(delta(pidx)))~=1)
            % get independent component
            t1=cell(n1,1); t2=cell(n2,1);
            for j=1:n1
                t1{j}=data(cidx1(j)).ind(:)...
                    -data(cidx1(j)).ind(1)+cb1(j);
            end
            for j=1:n2
                t2{j}=data(cidx2(j)).ind(:)...
                    -data(cidx2(j)).ind(1)+cb2(j);
            end
            
            % loop over all possible pairings
            for j=1:n1
                for k=1:n2
                    % get overlap info
                    [ob,oe,ol,t]=...
                        get_overlap(cb1(j),ce1(j),cb2(k),ce2(k),...
                        [t1{j}; t2{k}],true);
                    
                    % skip if unacceptable overlap
                    switch lower(option.OVERLAPUNITS)
                        case 'intervals' % actually samples
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
                            ' avg. delta: %f seconds\n'...
                            ' npts:  %d\n'...
                            'with\n'...
                            ' %d - %s\n'...
                            ' azimuth: %f degrees (%s)\n'...
                            ' begin: Day: %d Second: %f\n'...
                            ' end:   Day: %d Second: %f\n'...
                            ' avg. delta: %f seconds\n'...
                            ' npts:  %d\n'...
                            'Overlap: %f seconds (%d samples)\n\n'],...
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
                    
                    % interpolate values
                    ndata(c1).dep=interp1(...
                        t1{j},data(cidx1(j)).dep,t,option.INTERPOLATE);
                    ndata(c2).dep=interp1(...
                        t2{k},data(cidx2(k)).dep,t,option.INTERPOLATE);

                    % now rotate!!!
                    [ndata(c1).dep,ndata(c2).dep]=go_rotate(...
                        ndata(c1).dep,ndata(c2).dep,...
                        cmpaz(cidx1(j)),rf1,pto1,rf2);
                    
                    % update independent cmp
                    ndata(c1).ind=t+data(cidx1(j)).ind(1)-cb1(j);
                    ndata(c2).ind=t+data(cidx2(k)).ind(1)-cb2(k);
                    
                    % get updated header values
                    nnpts([c1 c2])=numel(ndata(c1).dep);
                    ndelta([c1 c2])=(t(end)-t(1))/nnpts(c1);
                    nnz([c1 c2],:)=nz([cidx1(j) cidx2(k)],:);
                    nkcmpnm([c1 c2])=kcmpnm([cidx1(j) cidx2(k)]);
                    nleven([c1 c2])={'false'};
                    ncmpaz([c1 c2])=[pto1 pto2];
                    
                    % get updated times
                    nab([c1 c2],1)=day1+floor(t(1)/86400);
                    nae([c1 c2],1)=day1+floor(t(end)/86400);
                    nab([c1 c2],2)=mod(t(1),86400);
                    nae([c1 c2],2)=mod(t(end),86400);
                    
                    % dep*
                    if(nnpts(c1))
                        ndepmin(c1)=min(ndata(c1).dep(:));
                        ndepmen(c1)=nanmean(ndata(c1).dep(:));
                        ndepmax(c1)=max(ndata(c1).dep(:));
                        ndepmin(c2)=min(ndata(c2).dep(:));
                        ndepmen(c2)=nanmean(ndata(c2).dep(:));
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
                t1{j}=cb1(j)+...
                    (0:delta(pidx(1)):delta(pidx(1))*(npts(cidx1(j))-1)).';
            end
            for j=1:n2
                t2{j}=cb2(j)+...
                    (0:delta(pidx(1)):delta(pidx(1))*(npts(cidx2(j))-1)).';
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
                    
                    % rotate
                    [t,ndata(c1).dep,ndata(c2).dep]=rotate_even(...
                        t1{j},data(cidx1(j)).dep,cmpaz(cidx1(j)),pto1,...
                        t2{k},data(cidx2(k)).dep,rf1,rf2,ob,oe,...
                        delta(pidx(1)),option);
                    
                    % get updated header values
                    nnz([c1 c2],:)=nz([cidx1(j) cidx2(k)],:);
                    nab([c1 c2],1)=day1; nae([c1 c2],1)=day1;
                    nab([c1 c2],2)=t(1); nae([c1 c2],2)=t(end);
                    nkcmpnm([c1 c2])=kcmpnm([cidx1(j) cidx2(k)]);
                    ndelta([c1 c2])=delta(pidx(1));
                    nleven([c1 c2])=leven(pidx(1));
                    ncmpaz([c1 c2])=[pto1 pto2];
                    
                    % dep*
                    nnpts([c1 c2])=numel(ndata(c1).dep);
                    if(nnpts(c1))
                        ndepmin(c1)=min(ndata(c1).dep(:));
                        ndepmen(c1)=nanmean(ndata(c1).dep(:));
                        ndepmax(c1)=max(ndata(c1).dep(:));
                        ndepmin(c2)=min(ndata(c2).dep(:));
                        ndepmen(c2)=nanmean(ndata(c2).dep(:));
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
    
    % trim off any remaining tmp space
    ndata(c2+1:end)=[]; nnz(c2+1:end,:)=[];
    nab(c2+1:end,:)=[]; nae(c2+1:end,:)=[];
    nleven(c2+1:end)=[]; nkcmpnm(c2+1:end)=[];
    ndelta(c2+1:end)=[]; nnpts(c2+1:end)=[]; ncmpaz(c2+1:end)=[];
    ndepmin(c2+1:end)=[]; ndepmen(c2+1:end)=[]; ndepmax(c2+1:end)=[];
    
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
    
    % fix kcmpnm
    nkcmpnm=cellstr(strnlen(char(nkcmpnm),2));
    nkcmpnm(1:2:end)=strcat(nkcmpnm(1:2:end),option.KCMPNM1);
    nkcmpnm(2:2:end)=strcat(nkcmpnm(2:2:end),option.KCMPNM2);
    
    % error if no output records
    if(~numel(ndata))
        error('seizmo:rotate:noRotatedRecords',...
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

function [t,x1,x2]=rotate_even(...
    t1,x1,az1,naz1,t2,x2,rf1,rf2,ob,oe,delta,option)

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
    t=t2(bidx2+(0:(npts-1)));
else % move 2
    t=t1(bidx1+(0:(npts-1)));
end

% get x1,x2
if(shift<=shiftmax)
    % ok we are just shifting one record to time align
    if(option.DEBUG); disp('Adjust Method: SHIFT'); end
    x1=x1(bidx1+(0:(npts-1)));
    x2=x2(bidx2+(0:(npts-1)));
else % interpolate
    % interpolating one record to time align
    if(option.DEBUG); disp('Adjust Method: INTERPOLATE'); end
    if(move1)
        x1=interp1(t1,x1,t,option.INTERPOLATE,'extrap');
        x2=x2(bidx2+(0:(npts-1)));
    else % interpolate 2
        x1=x1(bidx1+(0:(npts-1)));
        x2=interp1(t2,x2,t,option.INTERPOLATE,'extrap');
    end
end

% now rotate!!!
[x1,x2]=go_rotate(x1,x2,az1,rf1,naz1,rf2);

end

function [x1,x2]=go_rotate(x1,x2,az1,rf1,naz1,rf2)
% rotate x1,x2 where x1 points to az1 and rf1 indicates if
% x2 follows by 90deg.  naz1 gives the azimuth to rotate x1
% to and rf2 indicates if the output x2 follows the output
% x1 by 90deg (rf = true = follows)

% get rotation transformation matrix
m=rtm(az1,rf1,naz1,rf2);

% rotate by angle accounting for orientation reversals
x1=m*[x1 x2].';
x2=x1(2,:).';
x1=x1(1,:).';

end

function [m]=rtm(azin,rf1,azout,rf2)
%RTM    returns 2d rotation matrix
%
% azin  = cmp1 input azimuth in degrees
% rf1   = flag indicating if input cmp2 azimuth trails by 90 deg
% azout = cmp1 output azimuth in degrees
% rf2   = flag indicating if output cmp2 azimuth trails by 90 deg

% azimuth difference in radians
theta=pi/180*(azin-azout);

% rotation matrix
m=[           cos(theta)     (-1)^(1+rf1)*sin(theta);
   (-1)^(rf2)*sin(theta)   (-1)^(rf1+rf2)*cos(theta)];

end

function [ob,oe,ol,varargout]=get_overlap(b1,e1,b2,e2,delta,flag)
% get overlap range/length

if(flag) % unevenly sampled
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
