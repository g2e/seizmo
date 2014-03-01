function [data]=meld(data,varargin)
%MELD    Merge SEIZMO records
%
%    Usage:    data=meld(data)
%              data=meld(...,'tolerance',tolerance,...)
%              data=meld(...,'toleranceunits',units,...)
%              data=meld(...,'adjust',method,...)
%              data=meld(...,'overlap',method,...)
%              data=meld(...,'unevenoverlap',method,...)
%              data=meld(...,'gap',method,...)
%              data=meld(...,'shiftmax',value,...)
%              data=meld(...,'shiftunits',units,...)
%              data=meld(...,'interpolate',method,...)
%              data=meld(...,'filler',filler,...)
%              data=meld(...,'only',type,...)
%              data=meld(...,'skip',type,...)
%              data=meld(...,'useabsolutetiming',logical,...)
%              data=meld(...,'timing',standard,...)
%              data=meld(...,'requiredcharfields',fields,...)
%              data=meld(...,'requiredrealfields',fields,...)
%              data=meld(...,'allocate',size,...)
%              data=meld(...,'verbose',logical,...)
%              data=meld(...,'debug',logical,...)
%
%    Description:
%     DATA=MELD(DATA) will take all records in DATA and merge any pairs
%     that are within +/-0.02 seconds of being continuous.  The output
%     dataset will contain the merged records with all duplicate or partial
%     records removed.  A 'pair' must have a number of fields that are
%     identical to one another (see options REQUIREDCHARFIELDS and
%     REQUIREDREALFIELDS for a list of those).  By default the records are
%     merged end-to-end so no data is added or deleted -- just shifted to
%     be continuous in time.  Basically this is for merging data that is
%     already continuous as well as for eliminating small gaps & overlaps
%     (aka "time tears") associated with clock drift and digitization which
%     do not reflect the amount of drift over the data being merged.
%
%     DATA=MELD(...,'TOLERANCEUNITS',UNITS,...) allows changing the
%     units of the TOLERANCE option.  By default UNITS is 'SECONDS'.  This
%     may be changed to 'INTERVALS' (an interval being the time between 2
%     samples) if that is more useful.
%
%     DATA=MELD(...,'TOLERANCE',TOLERANCE,...) allows changing the
%     magnitude of the time tears that can be merged.  The default value
%     for TOLERANCE is 0.02 seconds (you can change the units of TOLERANCE
%     using the TOLERANCEUNITS option above).  This will limit the merges
%     to records that are within 0.02 seconds of being continuous (so gaps
%     and overlaps of up to 0.02 seconds are worked on).  Setting TOLERANCE
%     to 0.1 will therefore allow merging records with a gap or overlap
%     of up to 0.1 seconds (assuming option TOLERANCEUNITS is 'SECONDS').
%     TOLERANCE can also be a two-element vector, so that gaps or overlaps
%     around a certain magnitude can be targeted (use only positive numbers
%     as all time tears have a positive value).  This is particularly
%     suited for removing leap seconds that have been inserted into what
%     should be continuous data (see the Examples section below).
%
%     DATA=MELD(...,'ADJUST',METHOD,...) allows changing which record
%     out of a mergible pair is shifted/interpolated to time-align with the
%     other.  There are six choices: 'FIRST' 'LAST' 'LONGER' 'SHORTER'
%     'ONE' & 'TWO'.  The default is 'SHORTER' (which adjusts the shorter
%     record to time-align with the longer).  Method 'ONE' adjusts the
%     record with a lower index, while 'TWO' adjusts the higher.  The best
%     method probably varies with the situation.
%
%     DATA=MELD(...,'OVERLAP',METHOD,...) allows changing how
%     overlaps are merged.  There are 5 choices: 'SEQUENTIAL', 'TRUNCATE',
%     'ADD', 'AVERAGE' & 'BLEND'.  The default is 'SEQUENTIAL', which
%     shifts the timing of one of the records (as chosen by the ADJUST
%     option) so they no longer overlap and then combines the two records.
%     This is useful for cases where the data is actually continuous but a
%     time tear was inserted to deal with accrued time error over a longer
%     time section than the records that are being merged (like if the
%     digitizer has a clock drift).  The 'TRUNCATE' option allows for
%     deleting overlapping data from one of the records (as chosen by the
%     ADJUST option).  This is useful for merging records that do have some
%     redundant data (this does happen for data from IRIS).  The 'ADD'
%     option will merge the records much like the 'TRUNCATE' method, but
%     rather than dropping the overlapping data from one of the records,
%     the data is numerically added to the other records' data in the
%     overlapping segment.  This method is useful for merging what once
%     were sequential records that have been filtered with the final
%     conditions appended to those records (so that this will add the final
%     conditions to adjacent records).  The 'AVERAGE' option does a simple
%     average of the data in the overlapping segment so the records mend
%     together better.  This might be useful for some situations but seems
%     like a bad idea to me.  The 'BLEND' option is like the 'AVERAGE'
%     option but with a linear blend over the overlap to improve the
%     transition.  Like 'AVERAGE' this method seems like a bad idea to me
%     but you are the expert!
%
%     DATA=MELD(...,'UNEVENOVERLAP',METHOD,...) allows changing how
%     repeated points of uneven records are merged.  There are 4 choices:
%     'TRUNCATE', 'ADD', 'AVERAGE' & 'IGNORE'.  The default is 'TRUNCATE'.
%     The 'IGNORE' method keeps all the repeated points.  See the OVERLAP
%     option above for more details on the other methods.
%
%     DATA=MELD(...,'GAP',METHOD,...) allows changing how gaps are
%     merged.  There are three choices: 'SEQUENTIAL' 'INTERPOLATE' and
%     'FILL'.  The default is 'SEQUENTIAL', which basically just shifts
%     the timing of one of the records (as chosen by the ADJUST option) so
%     there no longer is a gap and then just combines the two records.  The
%     'INTERPOLATE' option allows for interpolation of the data gap.  The
%     style of interpolation can be set with the INTERPOLATE option (the
%     default is 'spline').  This option is useful for small data gaps that
%     are real (maybe from a small loss of data in a buffer) or artificial
%     (useful for interpolating over a removed glitch).  The 'FILL' option
%     will insert a filler value into the data gap (default is zeros).
%     This is useful for combining data with large gaps.  The filler can be
%     changed with the FILLER option.
%
%     DATA=MELD(...,'SHIFTMAX',MAXVALUE,...) allows changing the threshold
%     on when the remaining portion (the part which is not overlapping) of
%     a record-to-be-adjusted (see ADJUST option) is interpolated or
%     shifted to align (interval-wise) with the other record.  This option
%     only applies to the overlaps and gaps that ARE NOT to be made
%     sequential.  This means that if you are doing a truncation, gap
%     interpolation or gap filling this IS in effect.  A record in these
%     cases can be shifted at most one half the sample interval.  By
%     default MAXVALUE is 0.01 intervals -- meaning the adjusted record is
%     only shifted in time (without changing the dependent data) if the
%     time change necessary to align is less than a hundredth of the sample
%     interval.  Otherwise the data will be interpolated at the aligned
%     sample times (which is obviously slower due to the computation).  A
%     MAXVALUE of 0.5 intervals will always shift the data to the new times
%     without interpolating new values.  Really the choice depends on how
%     sensitive you think your data is to time shifts and/or how much you
%     trust the timing of the adjusted record.  If you're trying to get
%     relative arrival times of P recordings then you probably are worried
%     about minor shifts in timing (which begs the question of why you are
%     dealing with such crappy data in the first place).  The default is
%     basically saying the data timing is pretty accurate and any new data
%     from new time points should be interpolated unless the time shift is
%     damn small and really would not change anything.
%
%     DATA=MELD(...,'SHIFTUNITS',UNITS,...) allows changing the units
%     of the SHIFTMAX option.  By default UNITS is 'INTERVALS'.  This can
%     be changed to 'SECONDS' if that is more useful.
%
%     DATA=MELD(...,'INTERPOLATE',METHOD,...) allows changing the
%     interpolation method.  The choices are those allowed by the INTERP1
%     command: 'spline' 'pchip' 'linear' and 'nearest'.  The default is
%     'spline', which is continuous in the 1st and 2nd derivatives.  Look
%     out for artifacting if you use one of the other options and are going
%     to differentiate the data later.
%
%     DATA=MELD(...,'FILLER',FILLER,...) allows changing the filler
%     when the GAP option is set to 'FILL'.  The default is zero.  Can be
%     any real number or NaN.
%
%     DATA=MELD(...,'ONLY',TYPE,...) allows specifying which type of time
%     tear is merged: 'sequential' 'gaps' or 'overlaps'.  The default is []
%     (empty) and allows all types.  Do not use the ONLY option with the
%     SKIP option.  Note that 'sequential' applies to records that are
%     already sequential (have a time tear of 0), not gaps/overlaps to be
%     made sequential!
%
%     DATA=MELD(...,'SKIP',TYPE,...) allows specifying which type of time
%     tear is skipped: 'sequential' 'gaps' or 'overlaps'.  The default is
%     [] (empty) and skips none.  Do not use the SKIP option with the ONLY
%     option.  Note that 'sequential' only applies to records that are
%     already sequential (have a time tear of 0), not gaps/overlaps to be
%     made sequential!
%
%     DATA=MELD(...,'USEABSOLUTETIMING',LOGICAL,...) allows turning
%     on/off the usage of the reference time fields to figure out the
%     timing of data.  This can be safely turned off if all the data share
%     the same reference time (can be done for you using SYNCHRONIZE).
%     Leave it on if your reference times vary with each record.
%
%     DATA=MELD(...,'TIMING',STANDARD,...) allows changing the timing
%     standard assumed for the reference time.  The choices are: 'UTC' and
%     'TAI'.  The default is 'UTC' (oddly the standard used by digitizers,
%     databases and data formats) and does have UTC leap second awareness.
%     This is useful for merging data that have had UTC leap seconds
%     properly inserted (MELD actually won't even see the data "overlap"
%     often shown in programs that don't have leap second awareness because
%     the UTC times are converted to a leapless standard to carry out the
%     merge).  Proper handling of leap seconds requires that the records'
%     have their reference time at the actual UTC time.  If the recording
%     equipment doesn't actually handle leap seconds (and is thus a second
%     off until the clock is updated) then records after a leap second will
%     have a reference time 1 second early (because all leap seconds thus
%     far have been an additional second).  Some manual time adjustment may
%     be needed for the data if this is the case or see the Examples
%     section below for a simple usage form to deal exclusively with these
%     problematic time tears.  See LEAPSECONDS for even more info.  The
%     'TAI' option is useful for data without leap second concerns.
%
%     DATA=MELD(...,'REQUIREDCHARFIELDS',FIELDS,...) allows changing
%     the character fields required to be equal between records before
%     checking if they can be merged.  The list is a cellstring array.  The
%     default is: {'knetwk' 'kstnm' 'khole' 'kcmpnm'}.
%
%     DATA=MELD(...,'REQUIREDREALFIELDS',FIELDS,...) allows changing
%     the numerical fields required to be equal between records before
%     checking if they can be merged.  The list must be a cellstring array.
%     The default is: {'cmpinc' 'cmpaz'}.  Note that LEVEN and NCMP are
%     also required but cannot be removed from the list.  Adding DELTA to
%     the list will prevent creation of unevenly sampled records.
%
%     DATA=MELD(...,'ALLOCATE',SIZE,...) sets the temporary space
%     initially allocated for merged records.  This is just a guess of the
%     maximum number of merged records created for a group.  The default
%     value is 10.  Not really worth changing.
%
%     DATA=MELD(...,'VERBOSE',LOGICAL,...) turns on/off showing the meld
%     progress bar.  Useful for seeing how far along we are & how long
%     there is to finish.  Default is TRUE (on).
%
%     DATA=MELD(...,'DEBUG',LOGICAL,...) turns on/off detailed debugging
%     messages.  Default is FALSE (off).
%
%    Notes:
%     - MELD is a rather complicated function 'under the hood' as it tries
%       to allow for most of the sane to pathological ways I could come up
%       with to combine records.  The defaults should be good for the most
%       common case that the records should be continuous but if you have
%       you are faced with a more complicated situation, MELD likely can
%       save you some effort. Good luck!
%     - Biggest caveat -- merging multiple records together can return with
%       different solutions depending on the order of the records.  This is
%       because each merge operation is separate from the next (it works on
%       one pair of records at a time).  So which records gets truncated,
%       shifted, interpolated, or left alone depends on the order in which
%       the records are processed (starts at the lower indices) and the
%       ADJUST option.  So you may want to consider preparing the data a
%       bit before using meld if you are merging 3+ into 1+ records.
%       Sorting by start time and setting ADJUST to 'FIRST' would be enough
%       to give you a consistent result.
%     - Running FIXDELTA first to take care of small differences in sample
%       rates caused by floating point inaccuracies & clock drift allows
%       you to avoid the headache of handling merged data with variable
%       samplerates.  Using FIXDELTA can introduce timing issues though if
%       the clock drift is significant so user beware!
%     - If you receive an error while using MELD that you don't understand
%       or find a bug let me know.
%     - Want to speed MELD up?
%       - Are your reference times all the same? If yes, set
%         USEABSOLUTETIMING to FALSE.
%       - Do you care about leap seconds? If no, set TIMING to 'TAI'.
%       - Do you care about timing accuracy? If not too much, then consider
%         setting SHIFTMAX to 0.5 (with SHIFTUNITS set to INTERVALS).  This
%         allows nudging the timing of records by half an interval so that
%         they time-align without interpolating the data.  BIG speed jump
%         for the working with data that doesn't fit the default options.
%
%    Header changes: B, E, NPTS, DEPMEN, DEPMIN, DEPMAX
%                    (see CHECKHEADER for more)
%
%    Examples:
%     % Merge roughly 1 second gaps/overlaps:
%     data=meld(data,'tol',[0.99 1.01]);
%
%     % Merge just gaps:
%     data=meld(data,'only','gaps');
%
%     % Merge gaps by inserting NaNs:
%     data=meld(data,'gap','fill','filler',nan);
%
%     % Merge details (there is a lot!):
%     meld(data,'debug',true);
%
%     % Make & merge some data with different sample rates and overlap:
%     data=bseizmo(0:.1:1,rand(1,11),...
%                  1.1:.101:2.2,rand(1,11),...
%                  .5:.1:1.5,rand(1,11));
%     data=ch(data,'kstnm','R1');
%     mdata=merge(data,'debug',true);
%     plot0([data;mdata]);
%
%     % Compare defaults settings to fast (& cruder) options:
%     tic; meld(data); toc;
%     tic; meld(data,'shift',.5,'abs',false,'timing','tai'); toc;
%
%    See also: REMOVEDUPLICATES, CUT

%     Version History:
%        Dec.  6, 2008 - initial version
%        Dec.  8, 2008 - more options
%        Mar. 30, 2009 - major update: description added, tons of options, 
%                        handles the 'QSPA Case'
%        Apr.  1, 2009 - major update again: speedy multi-merge
%        Apr.  7, 2009 - slightly more detail in verbose output
%        Apr. 23, 2009 - fix seizmocheck for octave, move usage up
%        May  15, 2009 - minor doc update
%        May  28, 2009 - minor doc update
%        Sep. 30, 2009 - update for CMOD ==> LONMOD
%        Oct. 30, 2009 - significant update: overlap add/average method,
%                        improved time sequence code, handle partial
%                        pieces and dataless
%        Nov.  2, 2009 - toleranceunits, try/catch, seizmoverbose as
%                        default for verbose
%        Nov.  3, 2009 - minor doc update
%        Feb. 23, 2010 - moved verbose to debugging messages
%        Feb. 24, 2010 - fixed switch statements
%        Apr.  1, 2010 - progress bar shows nrecs
%        Feb. 11, 2011 - mass seizmocheck fix
%        Jan. 28, 2012 - drop SEIZMO global, better checkheader usage, made
%                        disp+sprintf use fprintf, parse rewrite allows
%                        more flexibility in option strings, only/skip
%                        options replace merge* options
%        Mar. 13, 2012 - use getheader improvements
%        Feb. 14, 2013 - use strcmpi for consistency (fixes a big bug in
%                        TIMING option (would treat 'UTC' as 'TAI' b/c it
%                        was not matching 'utc')
%        Feb. 20, 2014 - doc update, drop delta from requiredrealfields,
%                        progress bar now reprints to avoid clobbering
%                        debugging output, unevenoverlap option, uneven &
%                        multi-rate merging
%        Feb. 22, 2014 - overlap blend method, unevenoverlap ignore method
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:05 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:meld:badNumInputs',...
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

% attempt meld
try
    % parse p/v pairs
    option=parse_meld_parameters(varargin{:});

    % get full filenames (for debugging output)
    nrecs=numel(data);
    if(option.DEBUG)
        fullname=strcat({data.path}.',{data.name}.');
    end

    % get header fields
    if(option.USEABSOLUTETIMING)
        [b,e,ncmp,delta,npts,depmin,depmax,depmen,leven,...
            nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec]=getheader(data,...
            'b','e','ncmp','delta','npts','depmin','depmax','depmen',...
            'leven lgc','nzyear','nzjday','nzhour','nzmin','nzsec',...
            'nzmsec');
        dt=[nzyear nzjday nzhour nzmin nzsec nzmsec];
    else
        [b,e,ncmp,delta,npts,depmin,depmax,depmen,leven]=getheader(data,...
            'b','e','ncmp','delta','npts','depmin','depmax','depmen',...
            'leven lgc');
        dt=nan(nrecs,6);
    end
    szreal=size(option.REQUIREDREALFIELDS); reqreal=cell(szreal);
    szchar=size(option.REQUIREDCHARFIELDS); reqchar=cell(szchar);
    if(prod(szreal)~=0)
        [reqreal{:}]=getheader(data,option.REQUIREDREALFIELDS{:});
    end
    if(prod(szchar)~=0)
        [reqchar{:}]=getheader(data,option.REQUIREDCHARFIELDS{:});
    end
    leven=~strcmpi(leven,'false');

    % get start and end of records in absolute time
    if(option.USEABSOLUTETIMING)
        if(strcmpi(option.TIMING,'utc'))
            ab=gregorian2modserial(utc2tai(...
                [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000+b]));
            ae=gregorian2modserial(utc2tai(...
                [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000+e]));
        else
            ab=gregorian2modserial(...
                [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000+b]);
            ae=gregorian2modserial(...
                [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000+e]);
        end
    else
        ab=[zeros(nrecs,1) b];
        ae=[zeros(nrecs,1) e];
    end

    % change real to char
    for i=1:prod(szreal)
        reqreal{i}=num2str(reqreal{i},'%16.16e');
    end

    % make groups (require at least leven and ncmp to be the same)
    [f,h,h]=unique(char(strcat(strcat('',reqchar{:}),'_',...
        strcat('',reqreal{:}),'_',num2str(leven),'_',num2str(ncmp))),...
        'rows');
    ngrp=size(f,1);
    grppop=histc(h,1:ngrp);
    cumpop=cumsum(grppop);

    % add temp space to arrays
    alloc=(nrecs+1):(nrecs+option.ALLOCATE);
    data(nrecs+option.ALLOCATE).dep=[];
    ab(alloc,:)=nan; ae(alloc,:)=nan; dt(alloc,:)=nan;
    delta(alloc,1)=nan; npts(alloc,1)=nan; fullname(alloc,1)={''};
    depmen(alloc,1)=nan; depmin(alloc,1)=nan; depmax(alloc,1)=nan;

    % debug
    if(option.DEBUG)
        fprintf('Group IDs:\n')
        tmpids=strcat(num2str((1:ngrp)'),{'. '},f);
        fprintf('%s\n',tmpids{:});
        fprintf('\nLookup Table:\n')
        fprintf('(Record # - Group #)\n');
        fprintf('%d - %d\n',[1:nrecs; h.']);
    end
    
    % detail message
    if(option.VERBOSE)
        disp('Merging Record(s)');
        print_time_left(0,nrecs);
    end

    % loop through each group
    destroy=false(nrecs+option.ALLOCATE,1);
    for i=1:ngrp
        % get group member indices
        gidx=find(h==i);
        ng=numel(gidx);

        % backup for later
        origng=ng;

        % detail message
        if(option.DEBUG)
            fprintf('\nProcessing Group: %d\n',i);
            fprintf('Members: ');
            fprintf('%d ',gidx);
            fprintf('\nNumber in Group: %d\n',ng);
        end

        % find any exact timespan duplicates
        % - Assumes the underlying data is the same.
        %   - I could use dep* to check for this.
        bad=flagexactdupes(ab(gidx,:),ae(gidx,:));

        % detail message
        if(option.DEBUG)
            if(any(bad))
                fprintf('\nDeleting Duplicate(s):\n');
                fprintf(' %d',gidx(bad));
                fprintf('\nNumber Still in Group: %d\n',ng-sum(bad));
            end
        end

        % get rid of any exact timespan duplicates
        destroy(gidx(bad))=true;
        ng=ng-sum(bad);
        gidx(bad)=[];

        % no records to merge with
        if(ng==1)
            % detail message
            if(option.VERBOSE)
                if(option.DEBUG)
                    print_time_left(cumpop(i),nrecs,true);
                else
                    print_time_left(cumpop(i),nrecs);
                end
            end
            continue;
        end
        
        % separate arrays for adding new info
        history=num2cell(gidx).';
        newgidx=gidx;
        newidx=nrecs+1;
        newng=ng;

        % handle uneven / differing sample rates (NEEDS TO BE WRITTEN)
        if(~leven(gidx(1)) || numel(unique(delta(gidx)))~=1)
            % even or uneven?
            if(leven(gidx(1)))
                % detail message
                if(option.DEBUG)
                    fprintf('Group has multiple samplerates!\n');
                end
                
                % make even records uneven
                for j=1:ng
                    data(gidx(j)).ind=b(gidx(j))+(0:delta(gidx(j)):...
                        delta(gidx(j))*(npts(gidx(j))-1))';
                end
            else
                % detail message
                if(option.DEBUG)
                    fprintf('Group has uneven sampling!\n');
                end
            end
            
            % combine into 1 record with special handling for repeated pts
            [data(newidx),ab(newidx,:),ae(newidx,:),npts(newidx),...
                dt(newidx,:),fullname(newidx),history(newidx)]=...
                gomeld_uneven(data(gidx),ab(gidx,:),ae(gidx,:),...
                npts(gidx),dt(gidx,:),fullname(gidx),gidx,ng,option);
            
            % update arrays
            newng=newng+1;
            newgidx(newng)=newidx;
            delta(newidx)=diff(data(newidx).ind([end 1]))/(npts(newidx)-1);
            depmin(newidx)=min(data(newidx).dep(:));
            depmax(newidx)=max(data(newidx).dep(:));
            depmen(newidx)=nanmean(data(newidx).dep(:));
            leven(newidx)=false;
        else
            % all the same delta so share
            gdelta=delta(gidx(1));
            
            % loop until no files left unattempted
            attempted=npts(gidx)==0; % skip dataless
            while(any(~attempted))
                % get an unattempted file
                j=find(~attempted,1,'first');
                
                % go meld with other records
                [data(newidx),ab(newidx,:),ae(newidx,:),npts(newidx),...
                    dt(newidx,:),fullname(newidx),newhistory]=gomeld(...
                    data(gidx(j)),ab(gidx(j),:),ae(gidx(j),:),...
                    npts(gidx(j)),dt(gidx(j),:),fullname(gidx(j)),...
                    history(j),newidx,data(gidx),ab(gidx,:),ae(gidx,:),...
                    npts(gidx),dt(gidx,:),fullname(gidx),ng,...
                    history(1:ng),gidx,gdelta,option);
                
                % check meld history
                if(isequal(newhistory,history(j)))
                    attempted(j)=true;
                else
                    newng=newng+1;
                    attempted(ismember(gidx,[newhistory{:}]))=true;
                    newgidx(newng)=newidx;
                    delta(newidx)=gdelta;
                    history(newng)=newhistory;
                    depmin(newidx)=min(data(newidx).dep(:));
                    depmax(newidx)=max(data(newidx).dep(:));
                    depmen(newidx)=nanmean(data(newidx).dep(:));
                    leven(newidx)=true;
                    newidx=newidx+1;
                end
            end
        end

        % detail message
        if(option.DEBUG)
            fprintf('Finished Merging Group: %d\n',i);
            fprintf('Members: ');
            fprintf('%d ',newgidx);
            fprintf('\nNumber in Group: %d\n',newng);
        end

        % get longest records with unique time coverage
        good=~(flagdupes(ab(newgidx,:),ae(newgidx,:))...
            | flaghistorydupes(history));
        ngood=sum(good);
        goodidx=1:ngood;

        % detail message
        if(option.DEBUG)
            fprintf('Deleting Duplicate(s) and/or Partial Piece(s):\n');
            fprintf(' %d',newgidx(~good));
            fprintf('\nChanging Indices Of Good Record(s):\n');
            fprintf('%d ==> %d\n',[newgidx(good).'; ...
                newgidx(goodidx).']);
            fprintf('-------------------------------\n');
            fprintf('%d kept / %d made / %d original\n',...
                ngood,newng-ng,origng);
        end

        % get rid of any duplicates/partial pieces
        data(newgidx(goodidx))=data(newgidx(good));
        ab(newgidx(goodidx),:)=ab(newgidx(good),:);
        ae(newgidx(goodidx),:)=ae(newgidx(good),:);
        delta(newgidx(goodidx))=delta(newgidx(good));
        npts(newgidx(goodidx))=npts(newgidx(good));
        depmin(newgidx(goodidx))=depmin(newgidx(good));
        depmax(newgidx(goodidx))=depmax(newgidx(good));
        depmen(newgidx(goodidx))=depmen(newgidx(good));
        leven(newgidx(goodidx))=leven(newgidx(good));
        dt(newgidx(goodidx),:)=dt(newgidx(good),:);
        destroy(newgidx(ngood+1:end))=true;
        
        % detail message
        if(option.VERBOSE)
            if(option.DEBUG)
                print_time_left(cumpop(i),nrecs,true);
            else
                print_time_left(cumpop(i),nrecs);
            end
        end
    end

    % trim off temp space
    data(nrecs+1:end)=[];
    destroy(nrecs+1:end)=[];
    ab(nrecs+1:end,:)=[];
    ae(nrecs+1:end,:)=[];
    dt(nrecs+1:end,:)=[];
    delta(nrecs+1:end)=[];
    npts(nrecs+1:end)=[];
    depmin(nrecs+1:end)=[];
    depmax(nrecs+1:end)=[];
    depmen(nrecs+1:end)=[];
    leven(nrecs+1:end)=[];

    % get relative times from absolute
    if(option.USEABSOLUTETIMING)
        if(strcmpi(option.TIMING,'utc'))
            az=gregorian2modserial(utc2tai(...
                [dt(:,1:4) dt(:,5)+dt(:,6)/1000]));
            b=(ab(:,1)-az(:,1))*86400+(ab(:,2)-az(:,2));
            e=(ae(:,1)-az(:,1))*86400+(ae(:,2)-az(:,2));
        else
            az=gregorian2modserial(...
                [dt(:,1:4) dt(:,5)+dt(:,6)/1000]);
            b=(ab(:,1)-az(:,1))*86400+(ab(:,2)-az(:,2));
            e=(ae(:,1)-az(:,1))*86400+(ae(:,2)-az(:,2));
        end
    else
        b=ab(:,2);
        e=ae(:,2);
    end

    % update header
    data=changeheader(data,'b',b,'e',e,'delta',delta,'npts',npts,...
        'depmin',depmin,'depmax',depmax,'depmen',depmen,'leven',leven);

    % remove unwanted data
    data(destroy)=[];

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


function [mdata,mab,mae,mnpts,mdt,mname,mhistory]=gomeld_uneven(...
    data,ab,ae,npts,dt,name,gidx,ng,option)
%GOMELD_UNEVEN    Merge unevenly spaced records

% detail message
if(option.DEBUG)
    fprintf('\nMerging Uneven/Multi-Rate Records:\n');
    for i=1:ng
        fprintf(' %d - %s\n',gidx(i),name{i});
    end
end

% history is sorted by time order
[mhistory,mhistory]=sortrows(ab);

% Define succession based on ADJUST
switch lower(option.ADJUST)
    case 'first'
        king=flipud(mhistory);
    case 'last'
        king=mhistory;
    case 'longer'
        [king,king]=sort(npts,[],'descend');
    case 'shorter'
        [king,king]=sort(npts);
    case 'one'
        [king,king]=sort(gidx,[],'descend');
    otherwise % two
        [king,king]=sort(gidx);
end

% get reference times in mod serial
if(option.USEABSOLUTETIMING)
    if(strcmpi(option.TIMING,'utc'))
        az=gregorian2modserial(utc2tai([dt(:,1:4) dt(:,5)+dt(:,6)/1000]));
    else
        az=gregorian2modserial([dt(:,1:4) dt(:,5)+dt(:,6)/1000]);
    end
else
    az=zeros(nrecs,2);
end

% adjust .ind of subjects to be relative to the king
tdiff=(az(:,1)-az(king(1),1))*86400+(az(:,2)-az(king(1),2));
for i=1:ng
    data(i).ind=data(i).ind+tdiff(i);
end
    
% the king has a prince
mdata=data(king(1));
mab=ab(mhistory(1),:);
mae=sortrows(ae);
mae=mae(end,:);
mdt=dt(king(1),:);
mname=name(king(1));
mnpts0=sum(npts); % assuming no overlap
mdata.dep=cat(1,data(king).dep);
mdata.ind=cat(1,data(king).ind);

% how does the prince deal with disputes
switch lower(option.UNEVENOVERLAP)
    case 'truncate'
        % the prince solves the dispute based on nobility
        [mdata.ind,idx]=unique(mdata.ind,'first');
        mdata.dep=mdata.dep(idx);
    case 'add'
        % the prince puts the disputed positions together
        [times,idx]=sort(mdata.ind);
        uidx=[find(diff(times)); numel(times)]; % last occurance
        mdata.ind=times(uidx);
        deps=mdata.dep(idx);
        deps=deps(uidx);
        n=histc(times,mdata.ind);
        i=mat2cell(idx,n);
        for a=find(n>1)'
            deps(a)=sum(mdata.dep(i{a}));
        end
        mdata.dep=deps;
    case 'average'
        % the prince comprimises between the disputes
        [times,idx]=sort(mdata.ind);
        uidx=[find(diff(times)); numel(times)]; % last occurance
        mdata.ind=times(uidx);
        deps=mdata.dep(idx);
        deps=deps(uidx);
        n=histc(times,mdata.ind);
        i=mat2cell(idx,n);
        for a=find(n>1)'
            deps(a)=sum(mdata.dep(i{a}))/n(a);
        end
        mdata.dep=deps;
    case 'ignore'
        % the prince lets the dispute continue
        [mdata.ind,idx]=sort(mdata.ind);
        mdata.dep=mdata.dep(idx);
end

% update npts
mnpts=numel(mdata.dep);

% detail message
if(option.DEBUG)
    fprintf('%d Redundant Time Points Removed Using\n',mnpts0-mnpts);
    fprintf('UNEVENOVERLAP Method: %s\n',upper(option.UNEVENOVERLAP));
end

% the prince's geneology
mhistory={gidx(mhistory)'};

end


function [mdata,mab,mae,mnpts,mdt,mname,mhistory]=gomeld(...
    mdata,mab,mae,mnpts,mdt,mname,mhistory,mabsidx,...
    data,ab,ae,npts,dt,name,nrecs,history,absidx,delta,option)
%GOMELD    Merge evenly spaced records

% if record is unmerged then identify as the original
if(isscalar([mhistory{:}])); idx=[mhistory{:}];
else idx=mabsidx;
end

for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % skip records already merged in
    if(any(absidx(i)==[mhistory{:}])); continue; end
    
    % time tear (the shift required to make sequential)
    tdiff(1)=-delta+(ab(i,1)-mae(1))*86400+(ab(i,2)-mae(2));
    tdiff(2)=-delta+(mab(1)-ae(i,1))*86400+(mab(2)-ae(i,2));
    
    % use minimum shift (so we don't try to meld in strange ways)
    % lead = 1, record is after
    % lead = 2, record is before
    [lead,lead]=min(abs(tdiff));
    tdiff=tdiff(lead);
    
    % skip based on switches
    % SGO = 1x3 logical of SEQ, GAP, OV
    if((~option.SGO(1) && tdiff==0) ...
            || (~option.SGO(2) && tdiff>0) ...
            || (~option.SGO(3) && tdiff<0))
        continue;
    end
    
    % deal with units
    tol=option.TOLERANCE;
    if(strcmpi(option.TOLERANCEUNITS,'intervals')); tol=tol.*delta; end
    
    % check if within tolerance
    if(abs(tdiff)>=tol(1) && abs(tdiff)<=tol(2))
        % meld type
        if(~tdiff)
            meldfunc=@mseq;
            type='SEQUENTIAL';
        elseif(tdiff>0)
            meldfunc=@mgap;
            type='GAP';
        else
            meldfunc=@mlap;
            type='OVERLAP';
        end
    else
        % data not mergible
        continue;
    end
    
    % detail message
    if(option.DEBUG)
        fprintf(...
            ['\nMerging Record:\n'...
             ' %d - %s\n'...
             ' begin: Day: %d Second: %f\n'...
             ' end:   Day: %d Second: %f\n'...
             ' npts:  %d\n'...
             ' meld history: '...
             sprintf('%d ',[mhistory{:}]) '\n'...
             'with\n'...
             ' %d - %s\n'...
             ' begin: Day: %d Second: %f\n'...
             ' end:   Day: %d Second: %f\n'...
             ' npts:  %d\n'...
             ' meld history: '...
             sprintf('%d ',absidx(i)) '\n'...
             'Sample Interval: %f seconds\n'...
             'Time Tear:  %f seconds (%f samples)\n'...
             'Merge Type: %s\n'],...
            idx,mname{:},mab(1),mab(2),mae(1),mae(2),mnpts,...
            absidx(i),name{i},ab(i,1),ab(i,2),ae(i,1),ae(i,2),npts(i),...
            delta,tdiff,tdiff/delta,type);
    end
    
    % meld the records
    [mdata,mab,mae,mnpts,mdt,mname]=meldfunc(...
        [mdata; data(i)],[mab; ab(i,:)],[mae; ae(i,:)],delta,...
        [mnpts npts(i)],[mdt; dt(i,:)],option,lead,[idx absidx(i)],...
        [mname name(i)],tdiff);
    
    % update meld history
    if(lead==1); mhistory={[mhistory{:} history{i}]};
    else mhistory={[history{i} mhistory{:}]};
    end
    
    % detail message
    if(option.DEBUG)
        fprintf(...
            ['Output Record:\n'...
             ' %d - %s\n'...
             ' begin: Day: %d Second: %f\n'...
             ' end:   Day: %d Second: %f\n'...
             ' npts: %d\n'...
             ' meld history: ' sprintf('%d ',[mhistory{:}]) '\n'],...
            mabsidx,mname{:},mab(1),mab(2),mae(1),mae(2),mnpts);
    end
    
    % now recurse (try merging the new record with others)
    [mdata,mab,mae,mnpts,mdt,mname,mhistory]=gomeld(...
        mdata,mab,mae,mnpts,mdt,mname,mhistory,mabsidx,...
        data,ab,ae,npts,dt,name,nrecs,history,absidx,delta,option);
end

end


function [data,ab,ae,npts,dt,name]=mseq(...
    data,ab,ae,delta,npts,dt,option,first,idx,name,varargin)
%MSEQ    Merge sequential records

% who do we keep
last=3-first;
switch lower(option.ADJUST)
    case 'first'
        keep=last;
    case 'last'
        keep=first;
    case 'longer'
        [keep,keep]=min(npts);
    case 'shorter'
        [keep,keep]=max(npts);
    case 'one'
        keep=2;
    otherwise % two
        keep=1;
end
kill=3-keep;

% detail message
if(option.DEBUG)
    fprintf('Adjusting: %s (%d - %s)\n',...
        option.ADJUST,idx(kill),name{kill});
end

% adjust
if(keep==first)
    ab=ab(first,:);
    ae=fixmodserial(ae(first,:)+[0 npts(last)*delta]);
else
    ab=fixmodserial(ab(last,:)-[0 npts(first)*delta]);
    ae=ae(last,:);
end

% meld
data(keep).dep=[data(first).dep; data(last).dep];
data(kill)=[];

% update name, npts, dt
name=name(keep);
npts=sum(npts);
dt=dt(keep,:);
    
end


function [data,ab,ae,npts,dt,name]=mgap(...
    data,ab,ae,delta,npts,dt,option,first,idx,name,diff)
%MGAP    Merge gaps

% how much do we need to shift the samples
shift=lonmod(diff,delta);

% get max shift in seconds
if(isequal(option.SHIFTUNITS,'intervals'))
    maxshift=option.SHIFTMAX*delta;
else
    maxshift=option.SHIFTMAX;
end

% who do we keep
last=3-first;
switch lower(option.ADJUST)
    case 'first'
        keep=last;
    case 'last'
        keep=first;
    case 'longer'
        [keep,keep]=min(npts);
    case 'shorter'
        [keep,keep]=max(npts);
    case 'one'
        keep=2;
    otherwise % two
        keep=1;
end
kill=3-keep;

% detail message
if(option.DEBUG)
    fprintf('Gap Method: %s\n',option.GAP);
end

% how to handle gap?
switch lower(option.GAP)
    case 'sequential'
        % just shift the option.ADJUST record to make sequential
        [data,ab,ae,npts,dt,name]=mseq(...
            data,ab,ae,delta,npts,dt,option,first,idx,name);
    case 'interpolate'
        % detail message
        if(option.DEBUG)
            fprintf(...
                ['Interpolate Method: %s\n',...
                 'Adjusting: %s (%d - %s)\n'...
                 'Adjustment: %f seconds (%f samples)\n'],...
                option.INTERPOLATE,option.ADJUST,...
                idx(kill),name{kill},shift,shift/delta);
        end
        
        % update
        dt=dt(keep,:);
        name=name(keep);
        
        % shift or interpolate option.ADJUST record to align?
        if(abs(shift)<=maxshift)
            % detail message
            if(option.DEBUG)
                disp('Adjust Method: shift');
            end
            
            % shift which?
            if(kill==first)
                % shift timing of first record and interpolate gap
                [data,ab,ae,npts]=shiftfirstandinterpolategap(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            else
                % shift timing of last record and interpolate gap
                [data,ab,ae,npts]=shiftlastandinterpolategap(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            end
        else
            % detail message
            if(option.DEBUG)
                fprintf(...
                    ['Adjust Method: interpolate\n'...
                     'Interpolate Method: %s\n'],option.INTERPOLATE);
            end
            
            % interpolate which?
            if(kill==first)
                % interpolate first record and interpolate gap
                [data,ab,ae,npts]=interpolatefirstandinterpolategap(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            else
                % interpolate last record and interpolate gap
                [data,ab,ae,npts]=interpolatelastandinterpolategap(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            end
        end
    case 'fill'
        % detail message
        if(option.DEBUG)
            fprintf(...
                ['Filler: %d\n',...
                 'Adjusting: %s (%d - %s)\n'...
                 'Adjustment: %f seconds (%f samples)\n'],...
                option.FILLER,option.ADJUST,idx(kill),name{kill},...
                shift,shift/delta);
        end
        
        % update
        dt=dt(keep,:);
        name=name(keep);
        
        % shift or interpolate option.ADJUST record to align?
        if(abs(shift)<=maxshift)
            % detail message
            if(option.DEBUG)
                disp('Adjust Method: shift');
            end
            
            % shift which?
            if(kill==first)
                % shift timing of first record and fill gap
                [data,ab,ae,npts]=shiftfirstandfillgap(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            else
                % shift timing of last record and fill gap
                [data,ab,ae,npts]=shiftlastandfillgap(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            end
        else
            % detail message
            if(option.DEBUG)
                fprintf(...
                    ['Adjust Method: interpolate\n'...
                     'Interpolate Method: %s\n'],option.INTERPOLATE);
            end
            
            % interpolate which?
            if(kill==first)
                % interpolate first record and fill gap
                [data,ab,ae,npts]=interpolatefirstandfillgap(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            else
                % interpolate last record and fill gap
                [data,ab,ae,npts]=interpolatelastandfillgap(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            end
        end
end

end


function [data,ab,ae,npts,dt,name]=mlap(...
    data,ab,ae,delta,npts,dt,option,first,idx,name,diff)
%MLAP    Merge overlaps

% minimum time shift to align (always within +/-delta/2)
shift=lonmod(diff,delta);

% get max shift allowed in seconds (otherwise interpolate)
if(isequal(option.SHIFTUNITS,'intervals'))
    maxshift=option.SHIFTMAX*delta;
else
    maxshift=option.SHIFTMAX;
end

% who do we keep
last=3-first;
switch lower(option.ADJUST)
    case 'first'
        keep=last;
    case 'last'
        keep=first;
    case 'longer'
        [keep,keep]=min(npts);
    case 'shorter'
        [keep,keep]=max(npts);
    case 'one'
        keep=2;
    otherwise % two
        keep=1;
end
kill=3-keep;

% detail message
if(option.DEBUG)
    fprintf('Overlap Method: %s\n',option.OVERLAP);
end

% how to handle overlap?
switch lower(option.OVERLAP)
    case 'sequential'
        % just shift the option.ADJUST record to make sequential
        [data,ab,ae,npts,dt,name]=mseq(...
            data,ab,ae,delta,npts,dt,option,first,idx,name);
    case 'truncate'
        % truncate overlap from option.ADJUST record
        % detail message
        if(option.DEBUG)
            fprintf(...
                ['Adjusting:  %s (%d - %s)\n'...
                 'Adjustment: %f seconds (%f samples)\n'],...
                option.ADJUST,idx(kill),name{kill},shift,shift/delta);
        end
        
        % update
        dt=dt(keep,:);
        name=name(keep);
        
        % shift or interpolate option.ADJUST record to align?
        if(abs(shift)<=maxshift)
            % detail message
            if(option.DEBUG)
                disp('Adjust Method: shift');
            end
            
            % shift which?
            if(kill==first)
                % shift timing and truncate first record
                [data,ab,ae,npts]=shiftandtruncatefirst(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            else
                % shift timing and truncate last record
                [data,ab,ae,npts]=shiftandtruncatelast(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            end
        else
            % detail message
            if(option.DEBUG)
                fprintf(...
                    ['Adjust Method: interpolate\n'...
                     'Interpolate Method: %s\n'],option.INTERPOLATE);
            end
            
            % interpolate which?
            if(kill==first)
                % interpolate and truncate first record
                [data,ab,ae,npts]=interpolateandtruncatefirst(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            else
                % interpolate and truncate last record
                [data,ab,ae,npts]=interpolateandtruncatelast(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            end
        end
    case {'add' 'average'}
        % add/average overlap
        
        % divisor for adding (1) or averaging (2)
        d=1; if(strcmpi(option.OVERLAP,'average')); d=2; end
        
        % detail message
        if(option.DEBUG)
            fprintf(...
                ['Adjusting:  %s (%d - %s)\n'...
                 'Adjustment: %f seconds (%f samples)\n'],...
                option.ADJUST,idx(kill),name{kill},shift,shift/delta);
        end
        
        % update
        dt=dt(keep,:);
        name=name(keep);
        
        % shift or interpolate option.ADJUST record to align?
        if(abs(shift)<=maxshift)
            % detail message
            if(option.DEBUG)
                disp('Adjust Method: shift');
            end
            
            % shift which?
            if(kill==first)
                % shift timing and add first record
                [data,ab,ae,npts]=shiftandavgfirst(...
                    data,ab,ae,delta,npts,shift,first,last,option,d);
            else
                % shift timing and add last record
                [data,ab,ae,npts]=shiftandavglast(...
                    data,ab,ae,delta,npts,shift,first,last,option,d);
            end
        else
            % detail message
            if(option.DEBUG)
                fprintf(...
                    ['Adjust Method: interpolate\n'...
                     'Interpolate Method: %s\n'],option.INTERPOLATE);
            end
            
            % interpolate which?
            if(kill==first)
                % interpolate and add first record
                [data,ab,ae,npts]=interpolateandavgfirst(...
                    data,ab,ae,delta,npts,shift,first,last,option,d);
            else
                % interpolate and add last record
                [data,ab,ae,npts]=interpolateandavglast(...
                    data,ab,ae,delta,npts,shift,first,last,option,d);
            end
        end
    case 'blend'
        % blended average of overlap
        
        % detail message
        if(option.DEBUG)
            fprintf(...
                ['Adjusting:  %s (%d - %s)\n'...
                 'Adjustment: %f seconds (%f samples)\n'],...
                option.ADJUST,idx(kill),name{kill},shift,shift/delta);
        end
        
        % update
        dt=dt(keep,:);
        name=name(keep);
        
        % shift or interpolate option.ADJUST record to align?
        if(abs(shift)<=maxshift)
            % detail message
            if(option.DEBUG)
                disp('Adjust Method: shift');
            end
            
            % shift which?
            if(kill==first)
                % shift timing and add first record
                [data,ab,ae,npts]=shiftandblendfirst(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            else
                % shift timing and add last record
                [data,ab,ae,npts]=shiftandblendlast(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            end
        else
            % detail message
            if(option.DEBUG)
                fprintf(...
                    ['Adjust Method: interpolate\n'...
                     'Interpolate Method: %s\n'],option.INTERPOLATE);
            end
            
            % interpolate which?
            if(kill==first)
                % interpolate and add first record
                [data,ab,ae,npts]=interpolateandblendfirst(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            else
                % interpolate and add last record
                [data,ab,ae,npts]=interpolateandblendlast(...
                    data,ab,ae,delta,npts,shift,first,last,option);
            end
        end
end

end


function [time]=fixmodserial(time)
% fixes modserial times to have seconds between 0 and 86400
time=[time(:,1)+floor(time(:,2)./86400) mod(time(:,2),86400)];
end


function [data,ab,ae,npts]=shiftfirstandinterpolategap(...
    data,ab,ae,delta,npts,shift,first,last,option)
% shift first record and interpolate the data gap

% we need to shift the samples of the first record to align with the last
% then we interpolate samples in the gap between

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
%ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% figure out number of samples in gap
nsamples=round((ab(last,2)-sae(2)-delta)/delta);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% interpolate gap
lasttimes=ab(last,2)+(0:delta:delta*(npts(last)-1));
gaptimes=sae(2)+(delta:delta:delta*nsamples);
newtimes=sab(2)+(0:delta:delta*(npts(first)-1));
gapdata=interp1([newtimes lasttimes],...
    [data(first).dep; data(last).dep],...
    gaptimes,option.INTERPOLATE,'extrap');
gapdata=gapdata.';

% combine data
data(last).dep=[data(first).dep; gapdata; data(last).dep];
data(first)=[];

% set new ab, ae, npts
ab=fixmodserial(sab);
ae=ae(last,:);
npts=npts(1)+npts(2)+nsamples;

end


function [data,ab,ae,npts]=shiftlastandinterpolategap(...
    data,ab,ae,delta,npts,shift,first,last,option)
% shift last record and interpolate the data gap

% we need to shift the samples of the last record to align with the first
% then we interpolate samples in the gap between

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of last
sab=ab(last,:)-[0 shift];
sae=ae(last,:)-[0 shift];

% figure out number of samples in gap
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% interpolate gap
firsttimes=ab(first,2)+(0:delta:delta*(npts(first)-1));
gaptimes=ae(first,2)+(delta:delta:delta*nsamples);
newtimes=sab(2)+(0:delta:delta*(npts(last)-1));
gapdata=interp1([firsttimes newtimes],...
    [data(first).dep; data(last).dep],...
    gaptimes,option.INTERPOLATE,'extrap');
gapdata=gapdata.';

% combine data
data(first).dep=[data(first).dep; gapdata; data(last).dep];
data(last)=[];

% set new ab, ae, npts
ab=ab(first,:);
ae=fixmodserial(sae);
npts=npts(1)+npts(2)+nsamples;

end


function [data,ab,ae,npts]=interpolatefirstandinterpolategap(...
    data,ab,ae,delta,npts,shift,first,last,option)
% interpolate first record and gap

% we need to interpolate the samples of the first record to align with the
% last, then we interpolate the gap between

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
%ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% figure out number of samples in gap
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% interpolate first record and gap
firsttimes=ab(first,2)+(0:delta:delta*(npts(first)-1));
lasttimes=ab(last,2)+(0:delta:delta*(npts(last)-1));
gaptimes=sae(2)+(delta:delta:delta*nsamples);
newtimes=sab(2)+(0:delta:delta*(npts(first)-1));
newdata=interp1([firsttimes lasttimes],...
    [data(first).dep; data(last).dep],...
    [newtimes gaptimes],option.INTERPOLATE,'extrap');
newdata=newdata.';

% combine data
data(last).dep=[newdata; data(last).dep];
data(first)=[];

% set new ab, ae, npts
ab=fixmodserial(sab);
ae=ae(last,:);
npts=npts(1)+npts(2)+nsamples;

end


function [data,ab,ae,npts]=interpolatelastandinterpolategap(...
    data,ab,ae,delta,npts,shift,first,last,option)
% interpolate last record and gap

% we need to interpolate the samples of the last record to align with the
% first, then we interpolate the gap between

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of last
sab=ab(last,:)-[0 shift];
sae=ae(last,:)-[0 shift];

% figure out number of samples in gap
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% interpolate last record and gap
firsttimes=ab(first,2)+(0:delta:delta*(npts(first)-1));
lasttimes=ab(last,2)+(0:delta:delta*(npts(last)-1));
gaptimes=ae(first,2)+(delta:delta:delta*nsamples);
newtimes=sab(2)+(0:delta:delta*(npts(last)-1));
newdata=interp1([firsttimes lasttimes],...
    [data(first).dep; data(last).dep],...
    [gaptimes newtimes],option.INTERPOLATE,'extrap');
newdata=newdata.';

% combine data
data(first).dep=[data(first).dep; newdata];
data(last)=[];

% set new ab, ae, npts
ab=ab(first,:);
ae=fixmodserial(sae);
npts=npts(1)+npts(2)+nsamples;

end


function [data,ab,ae,npts]=...
    shiftfirstandfillgap(data,ab,ae,delta,npts,shift,first,last,option)
% shift first record and fill gap

% we need to shift the samples of the first record to align with the last
% then we fill the gap between

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
%ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% figure out number of samples in gap
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% combine data and add filler
ncmp=size(data(first).dep,2);
data(last).dep=[data(first).dep; ...
    option.FILLER.*ones(nsamples,ncmp); data(last).dep];
data(first)=[];

% set new ab, ae, npts
ab=fixmodserial(sab);
ae=ae(last,:);
npts=npts(1)+npts(2)+nsamples;

end


function [data,ab,ae,npts]=...
    shiftlastandfillgap(data,ab,ae,delta,npts,shift,first,last,option)
% shift last record and fill gap

% we need to shift the samples of the last record to align with the first
% then we fill the gap between

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of last
sab=ab(last,:)-[0 shift];
sae=ae(last,:)-[0 shift];

% figure out number of samples in gap
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% combine data and add filler
ncmp=size(data(first).dep,2);
data(first).dep=[data(first).dep;...
    option.FILLER.*ones(nsamples,ncmp); data(last).dep];
data(last)=[];

% set new ab, ae, npts
ab=ab(first,:);
ae=fixmodserial(sae);
npts=npts(1)+npts(2)+nsamples;

end


function [data,ab,ae,npts]=interpolatefirstandfillgap(...
    data,ab,ae,delta,npts,shift,first,last,option)
% interpolate first record and fill gap

% we need to interpolate the samples of the first record to align with the
% last, then we fill the gap

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
%ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% interpolate first record
oldtimes=ab(first,2)+(0:delta:delta*(npts(first)-1));
newtimes=sab(2)+(0:delta:delta*(npts(first)-1));
data(first).dep=...
    interp1(oldtimes,data(first).dep,newtimes,option.INTERPOLATE,'extrap');
data(first).dep=data(first).dep.';

% figure out number of samples in gap
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% combine data and fill gap
ncmp=size(data(first).dep,2);
data(last).dep=[data(first).dep;...
    option.FILLER.*ones(nsamples,ncmp); data(last).dep];
data(first)=[];

% set new ab, ae, npts
ab=fixmodserial(sab);
ae=ae(last,:);
npts=npts(1)+npts(2)+nsamples;

end


function [data,ab,ae,npts]=interpolatelastandfillgap(...
    data,ab,ae,delta,npts,shift,first,last,option)
% interpolate last record and fill gap

% we need to interpolate the samples of the last record to align with the
% first, then we fill the gap

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of last
sab=ab(last,:)-[0 shift];
sae=ae(last,:)-[0 shift];

% interpolate last record
oldtimes=ab(last,2)+(0:delta:delta*(npts(last)-1));
newtimes=sab(2)+(0:delta:delta*(npts(last)-1));
data(last).dep=...
    interp1(oldtimes,data(last).dep,newtimes,option.INTERPOLATE,'extrap');
data(last).dep=data(last).dep.';

% figure out number of samples in gap
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% combine data and fill gap
ncmp=size(data(first).dep,2);
data(first).dep=[data(first).dep;...
    option.FILLER.*ones(nsamples,ncmp); data(last).dep];
data(last)=[];

% set new ab, ae, npts
ab=ab(first,:);
ae=fixmodserial(sae);
npts=npts(1)+npts(2)+nsamples;

end


function [data,ab,ae,npts]=...
    shiftandtruncatefirst(data,ab,ae,delta,npts,shift,first,last,option)
% shift and truncate data of first record

% we need to shift the samples of the first record to align with the last
% then we drop the overlapping samples from the first

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% figure out number of samples dropped
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);
ns=min(nsamples,npts(last));

% detail message
if(option.DEBUG)
    fprintf('Samples Truncated: %d\n',nsamples);
end

% combine data and drop overlap from first
data(last).dep=[data(first).dep(1:end-nsamples,:); data(last).dep;
    data(first).dep(end-nsamples+ns+1:end,:)];
data(first)=[];

% set new ab, ae, npts
ab=fixmodserial(sab);
ae=fixmodserial([ab(first,1) max(ae(last,2),sae(2))]);
npts=npts(1)+npts(2)-ns;

end


function [data,ab,ae,npts]=...
    shiftandtruncatelast(data,ab,ae,delta,npts,shift,first,last,option)
% shift and truncate data of last record

% we need to shift the samples of the last record to align with the first
% then we drop the overlapping samples from the last 

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of last
sab=ab(last,:)-[0 shift];
sae=ae(last,:)-[0 shift];

% figure out number of samples dropped
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);
nsamples=min(nsamples,npts(last));

% detail message
if(option.DEBUG)
    fprintf('Samples Truncated: %d\n',nsamples);
end

% combine data and drop overlap from last
data(first).dep=[data(first).dep; data(last).dep(nsamples+1:end,:)];
data(last)=[];

% set new ab, ae, npts
ab=ab(first,:);
ae=fixmodserial([ab(first,1) max(ae(first,2),sae(2))]);
npts=npts(1)+npts(2)-nsamples;

end


function [data,ab,ae,npts]=interpolateandtruncatefirst(...
    data,ab,ae,delta,npts,shift,first,last,option)
% interpolate and truncate data of first record

% we need to interpolate the samples of the first record to align with the
% last, then we drop the overlapping samples from the first

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% interpolate first record
oldtimes=ab(first,2)+(0:delta:delta*(npts(first)-1));
newtimes=sab(2)+(0:delta:delta*(npts(first)-1));
data(first).dep=...
    interp1(oldtimes,data(first).dep,newtimes,option.INTERPOLATE,'extrap');
data(first).dep=data(first).dep.';

% figure out number of samples dropped
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);
ns=min(nsamples,npts(last));

% detail message
if(option.DEBUG)
    fprintf('Samples Truncated: %d\n',ns);
end

% combine data and drop overlap from first
data(last).dep=[data(first).dep(1:end-nsamples,:); data(last).dep;
    data(first).dep(end-nsamples+ns+1:end,:)];
data(first)=[];

% set new ab, ae, npts
ab=fixmodserial(sab);
ae=fixmodserial([ab(first,1) max(ae(last,2),sae(2))]);
npts=npts(1)+npts(2)-ns;

end


function [data,ab,ae,npts]=interpolateandtruncatelast(...
    data,ab,ae,delta,npts,shift,first,last,option)
% interpolate and truncate data of last record

% we need to interpolate the samples of the last record to align with the
% first, then we drop the overlapping samples from the last

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of last
sab=ab(last,:)-[0 shift];
sae=ae(last,:)-[0 shift];

% interpolate last record
oldtimes=ab(last,2)+(0:delta:delta*(npts(last)-1));
newtimes=sab(2)+(0:delta:delta*(npts(last)-1));
data(last).dep=...
    interp1(oldtimes,data(last).dep,newtimes,option.INTERPOLATE,'extrap');
data(last).dep=data(last).dep.';

% figure out number of samples dropped
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);
nsamples=min(nsamples,npts(last));

% detail message
if(option.DEBUG)
    fprintf('Samples Truncated: %d\n',nsamples);
end

% combine data and drop overlap from last
data(first).dep=[data(first).dep; data(last).dep(nsamples+1:end,:)];
data(last)=[];

% set new ab, ae, npts
ab=ab(first,:);
ae=fixmodserial([ab(first,1) max(ae(first,2),sae(2))]);
npts=npts(1)+npts(2)-nsamples;

end


function [data,ab,ae,npts]=...
    shiftandavgfirst(data,ab,ae,delta,npts,shift,first,last,option,dvd)
% shift and add data of first record

% we need to shift the samples of the first record to align with the last
% then we add the overlapping samples from the first to the last

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% figure out number of samples dropped
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);
e1=max(0,npts(first)-nsamples);
b2=min(nsamples,npts(last));
ns=b2+min(0,npts(first)-nsamples);

% detail message
if(option.DEBUG)
    fprintf('Samples Truncated: %d\n',nsamples);
end

% combine data
% - note that there may be partial pieces here
data(last).dep=[data(first).dep(1:e1,:); data(last).dep(1:b2-ns,:); ...
    (data(first).dep(e1+1:e1+ns,:)+data(last).dep(b2-ns+1:b2,:))/dvd; ...
    data(first).dep(e1+ns+1:end,:); data(last).dep(b2+1:end,:)];
data(first)=[];

% set new ab, ae, npts
ab=fixmodserial([ab(first,1) min(ab(last,2),sab(2))]);
ae=fixmodserial([ab(first,1) max(ae(last,2),sae(2))]);
npts=npts(1)+npts(2)-ns;

end


function [data,ab,ae,npts]=...
    shiftandavglast(data,ab,ae,delta,npts,shift,first,last,option,dvd)
% shift and add data of last record

% we need to shift the samples of the last record to align with the first
% then we add the overlapping samples from the last to the first

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of last
sab=ab(last,:)-[0 shift];
sae=ae(last,:)-[0 shift];

% figure out number of samples
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);
e1=max(0,npts(first)-nsamples);
b2=min(nsamples,npts(last));
ns=b2+min(0,npts(first)-nsamples);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% combine data
% - note that there may be partial pieces here
data(first).dep=[data(first).dep(1:e1,:); data(last).dep(1:b2-ns,:); ...
    (data(first).dep(e1+1:e1+ns,:)+data(last).dep(b2-ns+1:b2,:))/dvd; ...
    data(first).dep(e1+ns+1:end,:); data(last).dep(b2+1:end,:)];
data(last)=[];

% set new ab, ae, npts
ab=fixmodserial([ab(first,1) min(ab(first,2),sab(2))]);
ae=fixmodserial([ab(first,1) max(ae(first,2),sae(2))]);
npts=npts(1)+npts(2)-ns;

end


function [data,ab,ae,npts]=interpolateandavgfirst(...
    data,ab,ae,delta,npts,shift,first,last,option,dvd)
% interpolate and add data of first record

% we need to interpolate the samples of the first record to align with the
% last, then we add the overlapping samples from the first to the last

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% interpolate first record
oldtimes=ab(first,2)+(0:delta:delta*(npts(first)-1));
newtimes=sab(2)+(0:delta:delta*(npts(first)-1));
data(first).dep=...
    interp1(oldtimes,data(first).dep,newtimes,option.INTERPOLATE,'extrap');
data(first).dep=data(first).dep.';

% figure out number of samples
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);
e1=max(0,npts(first)-nsamples);
b2=min(nsamples,npts(last));
ns=b2+min(0,npts(first)-nsamples);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% combine data
% - note that there may be partial pieces here
data(last).dep=[data(first).dep(1:e1,:); data(last).dep(1:b2-ns,:); ...
    (data(first).dep(e1+1:e1+ns,:)+data(last).dep(b2-ns+1:b2,:))/dvd; ...
    data(first).dep(e1+ns+1:end,:); data(last).dep(b2+1:end,:)];
data(first)=[];

% set new ab, ae, npts
ab=fixmodserial([ab(first,1) min(ab(last,2),sab(2))]);
ae=fixmodserial([ab(first,1) max(ae(last,2),sae(2))]);
npts=npts(1)+npts(2)-ns;

end


function [data,ab,ae,npts]=interpolateandavglast(...
    data,ab,ae,delta,npts,shift,first,last,option,dvd)
% interpolate and add data of last record

% we need to interpolate the samples of the last record to align with the
% first, then we add the overlapping samples from the last to the first

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of last
sab=ab(last,:)-[0 shift];
sae=ae(last,:)-[0 shift];

% interpolate last record
oldtimes=ab(last,2)+(0:delta:delta*(npts(last)-1));
newtimes=sab(2)+(0:delta:delta*(npts(last)-1));
data(last).dep=...
    interp1(oldtimes,data(last).dep,newtimes,option.INTERPOLATE,'extrap');
data(last).dep=data(last).dep.';

% figure out number of samples
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);
e1=max(0,npts(first)-nsamples);
b2=min(nsamples,npts(last));
ns=b2+min(0,npts(first)-nsamples);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% combine data
% - note that there may be partial pieces here
data(first).dep=[data(first).dep(1:e1,:); data(last).dep(1:b2-ns,:); ...
    (data(first).dep(e1+1:e1+ns,:)+data(last).dep(b2-ns+1:b2,:))/dvd; ...
    data(first).dep(e1+ns+1:end,:); data(last).dep(b2+1:end,:)];
data(last)=[];

% set new ab, ae, npts
ab=fixmodserial([ab(first,1) min(ab(first,2),sab(2))]);
ae=fixmodserial([ab(first,1) max(ae(first,2),sae(2))]);
npts=npts(1)+npts(2)-ns;

end


function [data,ab,ae,npts]=...
    shiftandblendfirst(data,ab,ae,delta,npts,shift,first,last,option)
% shift and blend data of first record

% we need to shift the samples of the first record to align with the last
% then we blend the overlapping samples from the first to the last

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% figure out number of samples dropped
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);
e1=max(0,npts(first)-nsamples);
b2=min(nsamples,npts(last));
ns=b2+min(0,npts(first)-nsamples);

% detail message
if(option.DEBUG)
    fprintf('Samples Truncated: %d\n',nsamples);
end

% combine data
% - note that there may be partial pieces here
data(last).dep=[data(first).dep(1:e1,:); data(last).dep(1:b2-ns,:); ...
    (data(first).dep(e1+1:e1+ns,:).*(ns:-1:1).'...
    +data(last).dep(b2-ns+1:b2,:).*(1:ns).')/(ns+1); ...
    data(first).dep(e1+ns+1:end,:); data(last).dep(b2+1:end,:)];
data(first)=[];

% set new ab, ae, npts
ab=fixmodserial([ab(first,1) min(ab(last,2),sab(2))]);
ae=fixmodserial([ab(first,1) max(ae(last,2),sae(2))]);
npts=npts(1)+npts(2)-ns;

end


function [data,ab,ae,npts]=...
    shiftandblendlast(data,ab,ae,delta,npts,shift,first,last,option)
% shift and blend data of last record

% we need to shift the samples of the last record to align with the first
% then we blend the overlapping samples from the last to the first

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of last
sab=ab(last,:)-[0 shift];
sae=ae(last,:)-[0 shift];

% figure out number of samples
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);
e1=max(0,npts(first)-nsamples);
b2=min(nsamples,npts(last));
ns=b2+min(0,npts(first)-nsamples);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% combine data
% - note that there may be partial pieces here
data(first).dep=[data(first).dep(1:e1,:); data(last).dep(1:b2-ns,:); ...
    (data(first).dep(e1+1:e1+ns,:).*(ns:-1:1).'...
    +data(last).dep(b2-ns+1:b2,:).*(1:ns).')/(ns+1); ...
    data(first).dep(e1+ns+1:end,:); data(last).dep(b2+1:end,:)];
data(last)=[];

% set new ab, ae, npts
ab=fixmodserial([ab(first,1) min(ab(first,2),sab(2))]);
ae=fixmodserial([ab(first,1) max(ae(first,2),sae(2))]);
npts=npts(1)+npts(2)-ns;

end


function [data,ab,ae,npts]=interpolateandblendfirst(...
    data,ab,ae,delta,npts,shift,first,last,option)
% interpolate and blend data of first record

% we need to interpolate the samples of the first record to align with the
% last, then we blend the overlapping samples from the first to the last

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% interpolate first record
oldtimes=ab(first,2)+(0:delta:delta*(npts(first)-1));
newtimes=sab(2)+(0:delta:delta*(npts(first)-1));
data(first).dep=...
    interp1(oldtimes,data(first).dep,newtimes,option.INTERPOLATE,'extrap');
data(first).dep=data(first).dep.';

% figure out number of samples
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);
e1=max(0,npts(first)-nsamples);
b2=min(nsamples,npts(last));
ns=b2+min(0,npts(first)-nsamples);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% combine data
% - note that there may be partial pieces here
data(last).dep=[data(first).dep(1:e1,:); data(last).dep(1:b2-ns,:); ...
    (data(first).dep(e1+1:e1+ns,:).*(ns:-1:1).'...
    +data(last).dep(b2-ns+1:b2,:).*(1:ns).')/(ns+1); ...
    data(first).dep(e1+ns+1:end,:); data(last).dep(b2+1:end,:)];
data(first)=[];

% set new ab, ae, npts
ab=fixmodserial([ab(first,1) min(ab(last,2),sab(2))]);
ae=fixmodserial([ab(first,1) max(ae(last,2),sae(2))]);
npts=npts(1)+npts(2)-ns;

end


function [data,ab,ae,npts]=interpolateandblendlast(...
    data,ab,ae,delta,npts,shift,first,last,option)
% interpolate and blend data of last record

% we need to interpolate the samples of the last record to align with the
% first, then we blend the overlapping samples from the last to the first

% make sure all times share the same day
ae(first,:)=[ab(first,1) ae(first,2)+86400*(ae(first,1)-ab(first,1))];
ab(last,:)=[ab(first,1) ab(last,2)+86400*(ab(last,1)-ab(first,1))];
ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of last
sab=ab(last,:)-[0 shift];
sae=ae(last,:)-[0 shift];

% interpolate last record
oldtimes=ab(last,2)+(0:delta:delta*(npts(last)-1));
newtimes=sab(2)+(0:delta:delta*(npts(last)-1));
data(last).dep=...
    interp1(oldtimes,data(last).dep,newtimes,option.INTERPOLATE,'extrap');
data(last).dep=data(last).dep.';

% figure out number of samples
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);
e1=max(0,npts(first)-nsamples);
b2=min(nsamples,npts(last));
ns=b2+min(0,npts(first)-nsamples);

% detail message
if(option.DEBUG)
    fprintf('Samples Added: %d\n',nsamples);
end

% combine data
% - note that there may be partial pieces here
data(first).dep=[data(first).dep(1:e1,:); data(last).dep(1:b2-ns,:); ...
    (data(first).dep(e1+1:e1+ns,:).*(ns:-1:1).'...
    +data(last).dep(b2-ns+1:b2,:).*(1:ns).')/(ns+1); ...
    data(first).dep(e1+ns+1:end,:); data(last).dep(b2+1:end,:)];
data(last)=[];

% set new ab, ae, npts
ab=fixmodserial([ab(first,1) min(ab(first,2),sab(2))]);
ae=fixmodserial([ab(first,1) max(ae(first,2),sae(2))]);
npts=npts(1)+npts(2)-ns;

end


function [flags]=flagexactdupes(b,e)
% flags records that start/end at the same time as another as duplicates

% check inputs
if(~isequal(size(b),size(e)))
    error('seizmo:meld:badInputs',['Array sizes do not match:' ...
        'B: %s\nE: %s'],size(b),size(e))
end

% number of inputs
nrecs=size(b,1);

% assure times are valid
b=fixmodserial(b);
e=fixmodserial(e);

% find exact duplicates (note datetime must be valid)
dups=(b(:,ones(nrecs,1))==b(:,ones(nrecs,1)).' ...
    & b(:,2*ones(nrecs,1))==b(:,2*ones(nrecs,1)).') ...
    & (e(:,ones(nrecs,1))==e(:,ones(nrecs,1)).' ...
    & e(:,2*ones(nrecs,1))==e(:,2*ones(nrecs,1)).');
dupsu=(dups-dups.')>0;      % only delete if other not deleted
dupsu(tril(true(nrecs)))=0; % upper triangle
dups(triu(true(nrecs)))=0;  % lower triangle
flags=sum(dups+dupsu,2)>0;  % logical indices in group
        
end


function [flags]=flagdupes(b,e)
% flags records that are exact duplicates or are partial pieces

% check inputs
if(~isequal(size(b),size(e)))
    error('seizmo:meld:badInputs',['Array sizes do not match:' ...
        'B: %s\nE: %s'],size(b),size(e))
end

% number of inputs
nrecs=size(b,1);

% assure times are valid
b=fixmodserial(b);
e=fixmodserial(e);

% find duplicate data (note datetime must be valid)
dups=(b(:,ones(nrecs,1))>b(:,ones(nrecs,1)).' ...
    | (b(:,ones(nrecs,1))==b(:,ones(nrecs,1)).' ...
    & b(:,2*ones(nrecs,1))>=b(:,2*ones(nrecs,1)).')) ...
    & (e(:,ones(nrecs,1))<e(:,ones(nrecs,1)).' ...
    | (e(:,ones(nrecs,1))==e(:,ones(nrecs,1)).' ...
    & e(:,2*ones(nrecs,1))<=e(:,2*ones(nrecs,1)).'));
dupsu=(dups-dups.')>0;      % only delete if other not deleted
dupsu(tril(true(nrecs)))=0; % upper triangle
dups(triu(true(nrecs)))=0;  % lower triangle
flags=sum(dups+dupsu,2)>0;  % logical indices in group

end


function [flags]=flaghistorydupes(history)
% flags records that are the subset construct

% number of records
nrecs=numel(history);

% find history dupes
flags=false(nrecs,1);
for i=1:nrecs-1
    if(flags(i)); continue; end
    for j=i+1:nrecs
        if(flags(i)); break; end
        if(flags(j)); continue; end
        if(isempty(setdiff(history{j},history{i})))
            flags(j)=true;
        elseif(isempty(setdiff(history{i},history{j})))
            flags(i)=true;
        end
    end
end

end


function [option]=parse_meld_parameters(varargin)
% parses/checks meld p/v pairs

% valid values for string options
valid.TOLERANCEUNITS={'seconds' 'intervals'};
valid.OVERLAP={'sequential' 'truncate' 'add' 'average' 'blend'};
valid.UNEVENOVERLAP={'truncate' 'add' 'average' 'ignore'};
valid.GAP={'sequential' 'interpolate' 'fill'};
valid.INTERPOLATE={'spline' 'pchip' 'linear' 'nearest'};
valid.ADJUST={'longer' 'shorter' 'first' 'last' 'one' 'two'};
valid.SHIFTUNITS={'seconds' 'intervals'};
valid.TIMING={'utc' 'tai'};
valid.ONLY={'sequential' 'gaps' 'overlaps' 'seq' 's' 'g' 'o' 'ov'};
valid.SKIP={'sequential' 'gaps' 'overlaps' 'seq' 's' 'g' 'o' 'ov'};

% defaults
option.TOLERANCE=0.02; % seconds, any positive number
option.TOLERANCEUNITS='seconds'; % seconds/intervals
option.OVERLAP='sequential'; % sequential/truncate/add/average/blend
option.UNEVENOVERLAP='truncate'; % truncate/add/average/ignore
option.GAP='sequential'; % sequential/interpolate/fill
option.INTERPOLATE='spline'; % spline/pchip/linear/nearest
option.ADJUST='shorter'; % longer/shorter/first/last
option.SHIFTUNITS='intervals'; % seconds/intervals
option.SHIFTMAX=0.01; % interval: 0-0.5 , seconds: 0+
option.FILLER=0; % any number
option.TIMING='utc'; % utc/tai
option.USEABSOLUTETIMING=true; % true/false
option.REQUIREDCHARFIELDS={'knetwk' 'kstnm' 'khole' 'kcmpnm'};
option.REQUIREDREALFIELDS={'cmpinc' 'cmpaz'};
option.ALLOCATE=10; % size of temp space
option.ONLY=[]; % sequential/gaps/overlaps
option.SKIP=[]; % sequential/gaps/overlaps
option.VERBOSE=seizmoverbose; % default to seizmoverbose state
option.DEBUG=seizmodebug; % default to seizmodebug state

% require options specified by strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:meld:badInput',...
        'Not all options are specified with a string!');
end

% parse options
for i=1:2:nargin
    switch lower(varargin{i})
        case {'tolerance' 'tol' 't'}
            option.TOLERANCE=varargin{i+1};
        case {'toleranceunits' 'tunit' 'tu' 'tolu'}
            option.TOLERANCEUNITS=varargin{i+1};
        case {'overlap' 'over' 'lap' 'ov'}
            option.OVERLAP=varargin{i+1};
        case {'unevenoverlap' 'uneven' 'unover' 'unlap' 'unov'}
            option.UNEVENOVERLAP=varargin{i+1};
        case {'gap' 'ga' 'g'}
            option.GAP=varargin{i+1};
        case {'interpolate' 'interp' 'int' 'i'}
            option.INTERPOLATE=varargin{i+1};
        case {'adjust' 'adj' 'a'}
            option.ADJUST=varargin{i+1};
        case {'shiftunits' 'shiftu' 'sunit' 'su'}
            option.SHIFTUNITS=varargin{i+1};
        case {'shiftmax' 'shift' 'smax' 'max' 'sh'}
            option.SHIFTMAX=varargin{i+1};
        case {'filler' 'fill' 'f'}
            option.FILLER=varargin{i+1};
        case {'timing' 'time'}
            option.TIMING=varargin{i+1};
        case {'useabsolutetiming' 'useabs' 'abs' 'abstiming'}
            option.USEABSOLUTETIMING=varargin{i+1};
        case {'requiredcharfields' 'reqchar' 'charfields' 'cfields'}
            option.REQUIREDCHARFIELDS=varargin{i+1};
        case {'requiredrealfields' 'reqreal' 'realfields' 'rfields'}
            option.REQUIREDREALFIELDS=varargin{i+1};
        case {'allocate' 'all'}
            option.ALLOCATE=varargin{i+1};
        case 'only'
            option.ONLY=varargin{i+1};
        case 'skip'
            option.SKIP=varargin{i+1};
        case {'verbose' 'v'}
            option.VERBOSE=varargin{i+1};
        case {'debug' 'd'}
            option.DEBUG=varargin{i+1};
        otherwise
            error('seizmo:meld:badInput',...
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
        case 'tolerance'
            if(~isnumeric(value))
                error('seizmo:meld:badInput',...
                    'TOLERANCE must be a 1 or 2 element real array!');
            end
        case {'shiftmax' 'filler'}
            if(~isnumeric(value) || ~isscalar(value))
                error('seizmo:meld:badInput',...
                    '%s must be a scalar real number!',fields{i});
            end
        case {'overlap' 'gap' 'interpolate' 'toleranceunits' ...
                'adjust' 'shiftunits' 'timing' 'unevenoverlap'}
            if(~ischar(value) || size(value,1)~=1 ...
                    || ~any(strcmpi(value,valid.(fields{i}))))
                error('seizmo:meld:badInput',...
                    ['%s option must be one of the following:\n'...
                    sprintf('%s ',valid.(fields{i}){:})],fields{i});
            end
        case {'skip' 'only'}
            if(~isempty(value) && (~ischar(value) || size(value,1)~=1 ...
                    || ~any(strcmpi(value,valid.(fields{i})))))
                error('seizmo:meld:badInput',...
                    ['%s option must be one of the following:\n'...
                    sprintf('%s ',valid.(fields{i}){:})],fields{i});
            end
        case {'requiredcharfields' 'requiredrealfields'}
            % fix char arrays
            if(ischar(value))
                value=cellstr(value);
                option.(fields{i})=value;
            end
            if(~iscellstr(value))
                error('seizmo:meld:badInput',...
                    '%s option must be a cellstr of header fields!',...
                    fields{i});
            end
        case 'allocate'
            if(~isnumeric(value) || fix(value)~=value)
                error('seizmo:meld:badInput',...
                    'ALLOCATE must be a scalar integer!');
            end
        case {'useabsolutetiming' 'verbose' 'debug'}
            if(~islogical(value) || ~isscalar(value))
                error('seizmo:meld:badInput',...
                    '%s option must be a logical!',fields{i});
            end
        otherwise
            error('seizmo:meld:badInput',...
                'Unknown option: %s !',fields{i});
    end
end

% turn off verbose if debugging
if(option.DEBUG); option.VERBOSE=false; end

% handle tolerance
if(numel(option.TOLERANCE)==1)
    option.TOLERANCE=[-1 option.TOLERANCE];
end
option.TOLERANCE=sort(option.TOLERANCE);

% only 1 of skip/only can be set
if(~isempty(option.SKIP) && ~isempty(option.ONLY))
    error('seizmo:meld:badInput',...
        'Use options SKIP or ONLY, but not both!');
end

% set logical equivalent to skip/only string
option.SGO=true(1,3);
if(~isempty(option.ONLY))
    if(option.ONLY(1)=='s')
        option.SGO(2:3)=false;
    elseif(option.ONLY(1)=='g')
        option.SGO(1:2:3)=false;
    else % overlap
        option.SGO(1:2)=false;
    end
elseif(~isempty(option.SKIP))
    if(option.SKIP(1)=='s')
        option.SGO(1)=false;
    elseif(option.SKIP(1)=='g')
        option.SGO(2)=false;
    else % overlap
        option.SGO(3)=false;
    end
end

end

