function [data]=merge(data,varargin)
%MERGE    Merge SEIZMO records
%
%    Description: DATA=MERGE(DATA) will take all records in DATA and merge
%     any pairs that are within +/-0.02 seconds of being continuous.  The
%     output dataset will contain the merged records with all duplicate or
%     partial records removed.  A 'pair' must have a number of fields that
%     are identical to one another (see options REQUIREDCHARFIELDS and
%     REQUIREDREALFIELDS for a list of those).  By default the records are
%     merged end-to-end so no data is added or deleted -- just shifted to
%     be continuous in time.  Basically this is for eliminating small
%     'time tears'.
%
%     DATA=MERGE(DATA,...,'TOLERANCE',TOLERANCE,...) allows changing the
%     magnitude of the time tears that can be merged.  Setting TOLERANCE to
%     0.1 will therefore allow merging records with a gap or overlap within
%     +/-0.1 seconds.  TOLERANCE can also be a two-element vector, so that
%     gaps or overlaps around a certain magnitude can be targeted.  This is
%     particularly suited for removing leap seconds that have been stupidly
%     inserted into continuous data.
%
%     DATA=MERGE(DATA,...,'ADJUST',METHOD,...) allows changing which record
%     out of a mergible pair is shifted/interpolated to time-align with the
%     other.  There are six choices: 'FIRST' 'LAST' 'LONGER' 'SHORTER'
%     'ONE' & 'TWO.  The default is 'SHORTER' (which adjusts the shorter
%     record to time-align with the longer).  Method ONE adjust the record
%     with a lower index, while TWO adjust the higher.
%
%     DATA=MERGE(DATA,...,'OVERLAP',METHOD,...) allows changing how
%     overlaps are merged.  There are two choices: 'SEQUENTIAL' and
%     'TRUNCATE'.  The default is 'SEQUENTIAL', which basically just shifts
%     the timing of one of the records (as chosen by the ADJUST option) so
%     they no longer overlap and then combines just the two records.  This
%     is useful for cases where the data is actually continuous but a time
%     tear was inserted to deal with accrued time error.  The 'TRUNCATE'
%     option allows for deleting overlapping data from one of the records
%     (as chosen by the ADJUST option).  This is useful for merging records
%     that do have some redundant data.
%
%     DATA=MERGE(DATA,...,'GAP',METHOD,...) allows changing how gaps are
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
%     DATA=MERGE(DATA,...,'SHIFTMAX',MAXVALUE,...) allows changing the cap
%     on when the record-to-be-adjusted (see ADJUST option) is interpolated
%     or shifted to time-align with the other record.  This option only
%     applies to the overlaps and gaps that ARE NOT to be made sequential.
%     This means that if you are doing a truncation, gap interpolation or
%     gap filling the option is in effect.  A record in these cases can be
%     shifted at most one half the sample interval.  By default the
%     MAXVALUE is 0.01 intervals -- meaning the adjusted record is only
%     shifted if the time change is very minor, otherwise the data will be
%     interpolated at the aligned sample times.  A MAXVALUE of 0.5
%     intervals will always shift the data to the new times without
%     interpolating new values.  Really the choice depends on how sensitive
%     you think your data is to time shifts and/or how much you trust the
%     timing of the adjusted record.  If your trying to get relative
%     arrival times of P recordings then you probably are worried about
%     minor shifts in timing (which begs the question of why you are
%     dealing with such crappy data in the first place).  The default is
%     basically saying the data timing is pretty accurate and any new data
%     from new time points should be interpolated unless the time shift is
%     damn small and really would not change anything.
%
%     DATA=MERGE(DATA,...,'SHIFTUNITS',UNITS,...) allows changing the units
%     of the SHIFTMAX option.  By default UNITS is 'INTERVALS'.  This can
%     be changed to 'SECONDS' if that is more useful.
%
%     DATA=MERGE(DATA,...,'INTERPOLATE',METHOD,...) allows changing the
%     interpolation method.  The choices are basically those allowed in
%     Matlab's INTERP1 command: 'spline' 'pchip' 'linear' and 'nearest'.
%     The default is 'spline', which is continuous in the 1st and 2nd
%     derivatives.  Look out for artifacting if you use one of the other
%     options and are going to differentiate the data later.
%
%     DATA=MERGE(DATA,...,'FILLER',FILLER,...) allows changing the filler
%     when the GAP option is set to 'FILL'.  The default is zero.  Can be
%     any real number.
%
%     DATA=MERGE(DATA,...,'MERGESEQUENTIAL',LOGICAL,...) allows turning the
%     merging of sequential records (time tear of zero) on or off.  Will
%     not turn off making gaps or overlaps sequential (use GAP or OVERLAP
%     options above or see MERGEGAPS and MERGEOVERLAPS below).
%
%     DATA=MERGE(DATA,...,'MERGEOVERLAPS',LOGICAL,...) allows turning on/off
%     the merging of overlapping data that is within the time tear
%     tolerance.  Useful for just working on gaps.
%
%     DATA=MERGE(DATA,...,'MERGEGAPS',LOGICAL,...) allows turning on/off
%     the merging of data with gaps that are within the time tear
%     tolerance.  Useful for just working on overlaps.
%
%     DATA=MERGE(DATA,...,'USEABSOLUTETIMING',LOGICAL,...) allows turning
%     on/off the usage of the reference time fields to figure out the
%     timing of data.  This can be safely turned off if all the data share
%     the same reference time.  Leave it on if your reference times vary
%     with each record.
%
%     DATA=MERGE(DATA,...,'TIMING',STANDARD,...) allows changing the timing
%     standard assumed for the reference time.  The choices are: 'UTC' and
%     'TAI'.  The default is 'UTC', which has leap second awareness.  This
%     is useful for dealing with data that have had UTC leap seconds
%     properly inserted (basically MERGE won't even see the data overlap
%     because UTC times are converted to a leapless standard).  Proper
%     handling of leap seconds requires that the records' have their
%     reference time at the actual UTC time.  If the recording equipment
%     doesn't actually handle leap seconds then some time adjustment is/was
%     needed for the data.  See LEAPSECONDS for more info.  The 'TAI'
%     option is useful for data without any leap second concerns.
%
%     DATA=MERGE(DATA,...,'REQUIREDCHARFIELDS',FIELDS,...) allows changing
%     the character fields required to be equal between records before
%     checking if they can be merged.  The list is a cellstring array.  The
%     default is: {'knetwk' 'kstnm' 'khole' 'kcmpnm'}.
%
%     DATA=MERGE(DATA,...,'REQUIREDREALFIELDS',FIELDS,...) allows changing
%     the numerical fields required to be equal between records before
%     checking if they can be merged.  The list must be a cellstring array.
%     The default is: {'delta' 'cmpinc' 'cmpaz'}.  Note that LEVEN and NCMP
%     are also required but cannot be removed from the list.  Removing
%     DELTA from the list will allow creation of unevenly sampled records.
%
%     DATA=MERGE(DATA,...,'ALLOCATE',SIZE,...) sets the temporary space
%     initially allocated for merged records.  This is just a guess of the
%     maximum number of merged records created for a merge group.  The
%     default value is 10.  Not really worth changing.
%
%     DATA=MERGE(DATA,...,'VERBOSE',LOGICAL,...) turns on/off detailed
%     messages about the merges.  Useful for seeing what is happening.
%     Default is FALSE (off).
%
%     DATA=MERGE(DATA,...,'DEBUG',LOGICAL,...) turns on/off detailed
%     messages and some debugging messages.  Not really useful for anyone
%     but me.  Default is FALSE (off).
%
%    Notes:
%     - MERGE is a rather complicated function 'under the hood' as it tries
%       to allow for most of the sane ways I could come up with to combine
%       records.  The defaults should be good for the most common cases but
%       if you have a more complicated situation, MERGE (hopefully) can
%       save you some time. Good luck!
%     - Biggest caveat -- merging multiple records together can return with
%       different solutions depending on the order of the records.  This is
%       because each merge operation is separate from the next.  So which
%       records gets truncated, shifted, interpolated, or left alone
%       depends on the order in which the records are processed (starts at
%       the lower indices) and the ADJUST option.  So you may want to
%       consider preparing the data a bit before using merge if you are
%       merging 3+ into 1+ records.  Sorting by start time would probably
%       be enough.
%     - Run FIXDELTA first to take care of small differences in sample
%       rates caused by floating point inaccuracies!
%     - If you find an error you don't understand or a bug let me know!
%       This function hasn't been vetted as much as some of the others.
%     - How to speed things up?
%       - Are your reference times all the same? If yes, set
%         USEABSOLUTETIMING to FALSE.
%       - Do you care about leap seconds? If no, set TIMING to TAI.
%       - Do you care about timing accuracy? If not too much, then consider
%         setting SHIFTMAX to 0.5 (with SHIFTUNITS set to INTERVALS).  This
%         allows nudging the timing of records by half an interval so that
%         they time-align without interpolating the data.  BIG speed jump.
%
%    Tested on: Matlab r2007b
%
%    Header changes: B, E, NPTS, DEPMEN, DEPMIN, DEPMAX
%                    (see CHECKHEADER for more)
%
%    Usage:    data=merge(data)
%              data=merge(data,...,'tolerance',tolerance,...)
%              data=merge(data,...,'adjust',method,...)
%              data=merge(data,...,'overlap',method,...)
%              data=merge(data,...,'gap',method,...)
%              data=merge(data,...,'shiftmax',value,...)
%              data=merge(data,...,'shiftunits',units,...)
%              data=merge(data,...,'interpolate',method,...)
%              data=merge(data,...,'filler',filler,...)
%              data=merge(data,...,'mergesequential',logical,...)
%              data=merge(data,...,'mergeoverlaps',logical,...)
%              data=merge(data,...,'mergegaps',logical,...)
%              data=merge(data,...,'useabsolutetiming',logical,...)
%              data=merge(data,...,'timing',standard,...)
%              data=merge(data,...,'requiredcharfields',fields,...)
%              data=merge(data,...,'requiredrealfields',fields,...)
%              data=merge(data,...,'allocate',size,...)
%              data=merge(data,...,'verbose',logical,...)
%              data=merge(data,...,'debug',logical,...)
%
%    Examples:
%     Just friggin merge already!
%      data=merge(data)
%
%     Merge roughly 1 second gaps/overlaps:
%      data=merge(data,'tolerance',[0.99 1.01])
%
%     Merge just gaps:
%      data=merge(data,'mergesequential',false,'mergeoverlaps',false)
%
%     Merge gaps by inserting NaNs:
%      data=merge(data,'gap','fill','filler',nan)
%
%     Try it before you buy it:
%      merge(data,'verbose',true);
%
%    See also: removeduplicates, cut

%     Version History:
%        Dec.  6, 2008 - initial version
%        Dec.  8, 2008 - more options
%        Mar. 30, 2009 - major update: description added, tons of options, 
%                        handles the 'South Pole Case'
%        Apr.  1, 2009 - major update again: speedy multi-merge
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  1, 2009 at 08:20 GMT

% todo:
%   - uneven support - just toss together and sort after?
%   - debug!!

% check nargin
if(mod(nargin-1,2))
    error('seizmo:merge:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
error(seizmocheck(data,'dep'))

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% valid values for string options
valid.OVERLAP={'sequential' 'truncate'};
valid.GAP={'sequential' 'interpolate' 'fill'};
valid.INTERPOLATE={'spline' 'pchip' 'linear' 'nearest'};
valid.ADJUST={'longer' 'shorter' 'first' 'last' 'one' 'two'};
valid.SHIFTUNITS={'seconds' 'intervals'};
valid.TIMING={'utc' 'tai'};

% defaults
option.TOLERANCE=0.02; % seconds, any positive number
option.OVERLAP='sequential'; % sequential/truncate
option.GAP='sequential'; % sequential/interpolate/fill
option.INTERPOLATE='spline'; % spline/pchip/linear/nearest
option.ADJUST='shorter'; % longer/shorter/first/last
option.SHIFTUNITS='intervals'; % seconds/intervals
option.SHIFTMAX=0.01; % interval: 0-0.5 , seconds: 0+
option.FILLER=0; % any number
option.TIMING='utc'; % utc/tai
option.USEABSOLUTETIMING=true; % true/false
option.REQUIREDCHARFIELDS={'knetwk' 'kstnm' 'khole' 'kcmpnm'};
option.REQUIREDREALFIELDS={'delta' 'cmpinc' 'cmpaz'};
option.ALLOCATE=10; % size of temp space
option.MERGESEQUENTIAL=true; % on/off switch for merging sequential
option.MERGEGAPS=true; % on/off switch for merging gaps
option.MERGEOVERLAPS=true; % on/off switch for merging overlaps
option.VERBOSE=false; % turn on/off verbose messages
option.DEBUG=false; % turn on/off debugging messages

% get options from SEIZMO global
global SEIZMO
try
    fields=fieldnames(SEIZMO.MERGE);
    for i=1:numel(fields)
        option.(fields{i})=SEIZMO.MERGE.(fields{i});
    end
catch
end

% get options from command line
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error('seizmo:merge:badInput',...
            'Options must be specified as a strings!');
    end
    option.(upper(varargin{i}))=varargin{i+1};
end

% check options
fields=fieldnames(option);
for i=1:numel(fields)
    % get value of field and do a basic check
    value=option.(fields{i});
    if((ischar(value) && size(value,1)~=1))
        error('seizmo:merge:badInput',...
            'Bad value for option %s !',value);
    end
    
    % specific checks
    switch lower(fields{i})
        case 'tolerance'
            if(~isnumeric(value))
                error('seizmo:merge:badInput',...
                    'TOLERANCE must be a 1 or 2 element real array!');
            end
        case {'shiftmax' 'filler'}
            if(~isnumeric(value) || ~isscalar(value))
                error('seizmo:merge:badInput',...
                    '%s must be a scalar real number!',fields{i});
            end
        case {'overlap' 'gap' 'interpolate' 'adjust' 'shiftunits' 'timing'}
            if(~ischar(value) || size(value,1)~=1 || ~any(strcmpi(value,...
                    valid.(fields{i}))))
                error('seizmo:merge:badInput',...
                    ['%s option must be one of the following:\n'...
                    sprintf('%s ',valid.(fields{i}){:})],fields{i});
            end
        case {'requiredcharfields' 'requiredrealfields'}
            if(~iscellstr(value))
                error('seizmo:merge:badInput',...
                    '%s option must be a cellstr array of header fields!',...
                    fields{i});
            end
        case 'allocate'
            if(~isnumeric(value) || fix(value)~=value)
                error('seizmo:merge:badInput',...
                    'ALLOCATE must be a scalar integer!');
            end
        case {'useabsolutetiming' 'mergesequential' 'mergegaps'...
                'mergeoverlaps' 'verbose' 'debug'}
            if(~islogical(value) || ~isscalar(value))
                error('seizmo:merge:badInput',...
                    '%s option must be a logical!',fields{i});
            end
        otherwise
            error('seizmo:merge:badInput',...
                'Unknown option: %s !',fields{i});
    end
end

% handle tolerance
if(numel(option.TOLERANCE)==1)
    option.TOLERANCE=[-1 option.TOLERANCE];
end
option.TOLERANCE=sort(option.TOLERANCE);

% get full filenames (for verbose or debugging output)
nrecs=numel(data);
if(option.VERBOSE || option.DEBUG)
    fullname=strcat({data.path}.',{data.name}.');
end

% get header fields
ncmp=getncmp(data);
if(option.USEABSOLUTETIMING)
    [b,e,delta,npts,depmin,depmax,depmen,...
        nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec]=getheader(data,...
        'b','e','delta','npts','depmin','depmax','depmen',...
        'nzyear','nzjday','nzhour','nzmin','nzsec','nzmsec');
    dt=[nzyear nzjday nzhour nzmin nzsec nzmsec];
else
    [b,e,delta,npts,depmin,depmax,depmen]=...
        getheader(data,'b','e','delta','npts','depmin','depmax','depmen');
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
iftype=getenumid(data,'iftype');
leven=strcmp(getlgc(data,'leven'),'true');

% require timeseries and general x vs y
if(any(~strcmp(iftype,'itime') & ~strcmp(iftype,'ixy')))
    error('seizmo:merge:badRecordType',...
        'Records must be Time Series or General X vs Y !');
end

% get start and end of records in absolute time
if(option.USEABSOLUTETIMING)
    if(strcmp(option.TIMING,'utc'))
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

% add temp space to arrays
alloc=(nrecs+1):(nrecs+option.ALLOCATE);
data(nrecs+option.ALLOCATE).dep=[];
ab(alloc,:)=nan; ae(alloc,:)=nan; dt(alloc,:)=nan;
delta(alloc,1)=nan; npts(alloc,1)=nan; fullname(alloc,1)={''};
depmen(alloc,1)=nan; depmin(alloc,1)=nan; depmax(alloc,1)=nan;

% debug
if(option.DEBUG)
    disp('Group IDs:')
    disp(f);
    disp(' ');
    disp('Lookup Table:')
    disp(sprintf('%d - %d\n',[1:nrecs; h.']));
end

% loop through each group
destroy=false(nrecs+option.ALLOCATE,1);
for i=1:size(f,1)
    % get group member indices
    gidx=find(h==i);
    ng=numel(gidx);
    
    % backup for later
    origng=ng;
    
    % detail message
    if(option.VERBOSE || option.DEBUG)
        disp(' '); disp(' ');
        disp(sprintf('Processing Group: %d',i));
        disp(['Members: ' sprintf('%d ',gidx)]);
        disp(sprintf('Number in Group: %d',ng));
    end
    
    % find any exact duplicates
    bad=flagexactdupes(ab(gidx,:),ae(gidx,:));
    
    % detail message
    if(option.VERBOSE || option.DEBUG)
        disp(' ');
        disp('Deleting Duplicate(s):');
        disp(sprintf(' %d',gidx(bad)));
        disp(sprintf('Number Still in Group: %d',ng-sum(bad)));
    end
    
    % get rid of any exact duplicates
    destroy(gidx(bad))=true;
    ng=ng-sum(bad);
    gidx(bad)=[];
    
    % no records to merge with
    if(ng==1); continue; end
    
    % handle uneven / differing sample rates (NEEDS TO BE WRITTEN)
    if(any(~leven(gidx)) || numel(unique(delta(gidx)))~=1)
        % debug
        %any(~leven(gidx))
        %numel(unique(delta(gidx)))~=1
        %die
        
        % what to do?
        % - pass to mseq as if perfect
        % - merge .ind
        % - sort by .ind
        % - if diff .ind==0 whine
        error('seizmo:merge:unevenUnsupported',...
            ['Merging Uneven Records or Records with differing DELTA\n'...
             'is unsupported at the moment!']);
    end
    
    % all the same delta so share
    gdelta=delta(gidx(1));
    
    % separate arrays for adding new info
    history=num2cell(gidx).';
    newgidx=gidx;
    newng=ng;
    
    % loop until no files left unattempted
    attempted=false(ng,1); newidx=nrecs+1;
    while(any(~attempted))
        % get an unattempted file
        j=find(~attempted,1,'first');
        
        % go merge with other records
        [data(newidx),ab(newidx,:),ae(newidx,:),npts(newidx),...
            dt(newidx,:),fullname(newidx),newhistory]=gomerge(...
            data(gidx(j)),ab(gidx(j),:),ae(gidx(j),:),npts(gidx(j)),...
            dt(gidx(j),:),fullname(gidx(j)),history(j),newidx,...
            data(gidx),ab(gidx,:),ae(gidx,:),npts(gidx),dt(gidx,:),...
            fullname(gidx),ng,history(1:ng),gidx,gdelta,option);
        
        % check merge history
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
            depmen(newidx)=mean(data(newidx).dep(:));
            newidx=newidx+1;
        end
    end
    
    % detail message
    if(option.VERBOSE || option.DEBUG)
        disp(' ');
        disp(sprintf('Finished Merging Group: %d',i));
        disp(['Members: ' sprintf('%d ',newgidx)]);
        disp(sprintf('Number in Group: %d',newng));
    end
    
    % get longest records with unique time coverage
    good=~(flagdupes(ab(newgidx,:),ae(newgidx,:))...
        | flaghistorydupes(history));
    ngood=sum(good);
    goodidx=1:ngood;
    
    % detail message
    if(option.VERBOSE || option.DEBUG)
        disp('Deleting Duplicate(s) and/or Partial Piece(s):');
        disp(sprintf(' %d',newgidx(~good)));
        disp('Changing Indices Of Good Record(s):');
        disp(sprintf('%d ==> %d\n',[newgidx(good).'; newgidx(goodidx).']));
        disp('-------------------------------');
        disp(sprintf('%d kept / %d made / %d original',...
            ngood,newng-ng,origng));
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
    dt(newgidx(goodidx),:)=dt(newgidx(good),:);
    destroy(newgidx(ngood+1:end))=true;
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

% get relative times from absolute
if(option.USEABSOLUTETIMING)
    if(strcmp(option.TIMING,'utc'))
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
    'depmin',depmin,'depmax',depmax,'depmen',depmen);

% remove unwanted data
data(destroy)=[];

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end


function [mdata,mab,mae,mnpts,mdt,mname,mhistory]=gomerge(...
    mdata,mab,mae,mnpts,mdt,mname,mhistory,mabsidx,...
    data,ab,ae,npts,dt,name,nrecs,history,absidx,delta,option)

% if record is unmerged then identify as the original
if(isscalar([mhistory{:}])); idx=[mhistory{:}];
else idx=mabsidx;
end

for i=1:nrecs
    % skip records already merged in
    if(any(absidx(i)==[mhistory{:}])); continue; end
    
    % time tear (the shift required to make sequential)
    tdiff(1)=-delta+(ab(i,1)-mae(1))*86400+(ab(i,2)-mae(2));
    tdiff(2)=-delta+(mab(1)-ae(i,1))*86400+(mab(2)-ae(i,2));
    
    % use minimum shift (so we don't try to merge in strange ways)
    [lead,lead]=min(abs(tdiff));
    tdiff=tdiff(lead);
    
    % skip based on switches
    if((~option.MERGESEQUENTIAL && tdiff==0) ...
            || (~option.MERGEGAPS && tdiff>0) ...
            || (~option.MERGEOVERLAPS && tdiff<0))
        continue;
    end
    
    % check if within tolerance
    if(abs(tdiff)>=option.TOLERANCE(1)...
            && abs(tdiff)<=option.TOLERANCE(2))
        % merge type
        if(~tdiff)
            mergefunc=@mseq;
            type='SEQUENTIAL';
        elseif(tdiff>0)
            mergefunc=@mgap;
            type='GAP';
        else
            mergefunc=@mlap;
            type='OVERLAP';
        end
    else
        % data not mergible
        continue;
    end
    
    % detail message
    if(option.VERBOSE || option.DEBUG)
        disp(sprintf(...
            ['\nMerging Record:\n'...
             ' %d - %s\n'...
             ' begin: Day: %d Second: %f\n'...
             ' end:   Day: %d Second: %f\n'...
             ' npts:  %d\n'...
             ' merge history: '...
             sprintf('%d ',[mhistory{:}]) '\n'...
             'with\n'...
             ' %d - %s\n'...
             ' begin: Day: %d Second: %f\n'...
             ' end:   Day: %d Second: %f\n'...
             ' npts:  %d\n'...
             ' merge history: '...
             sprintf('%d ',absidx(i)) '\n'...
             'Time Tear:  %f seconds\n'...
             'Merge Type: %s'],...
            idx,mname{:},mab(1),mab(2),mae(1),mae(2),mnpts,...
            absidx(i),name{i},ab(i,1),ab(i,2),ae(i,1),ae(i,2),npts(i),...
            tdiff,type));
    end
    
    % merge the records
    [mdata,mab,mae,mnpts,mdt,mname]=mergefunc(...
        [mdata; data(i)],[mab; ab(i,:)],[mae; ae(i,:)],delta,...
        [mnpts npts(i)],[mdt; dt(i,:)],option,lead,[idx absidx(i)],...
        [mname name(i)],tdiff);
    
    % update merge history
    if(lead==1); mhistory={[mhistory{:} history{i}]};
    else mhistory={[history{i} mhistory{:}]};
    end
    
    % detail message
    if(option.VERBOSE || option.DEBUG)
        disp(sprintf(...
            ['Output Record:\n'...
             ' %d - %s\n'...
             ' begin: Day: %d Second: %f\n'...
             ' end:   Day: %d Second: %f\n'...
             ' npts: %d\n'...
             ' merge history: ' sprintf('%d ',[mhistory{:}])],...
            mabsidx,mname{:},mab(1),mab(2),mae(1),mae(2),mnpts));
    end
    
    % now recurse
    [mdata,mab,mae,mnpts,mdt,mname,mhistory]=gomerge(...
        mdata,mab,mae,mnpts,mdt,mname,mhistory,mabsidx,...
        data,ab,ae,npts,dt,name,nrecs,history,absidx,delta,option);
end

end


function [data,ab,ae,npts,dt,name]=mseq(...
    data,ab,ae,delta,npts,dt,option,first,idx,name,varargin)
%MSEQ    Merge sequential records

% who do we keep
last=3-first;
switch option.ADJUST
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
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Adjusting: %s (%d - %s)',...
        option.ADJUST,idx(kill),name{kill}));
end

% adjust
if(keep==first)
    ab=ab(first,:);
    ae=fixmodserial(ae(first,:)+[0 npts(last)*delta]);
else
    ab=fixmodserial(ab(last,:)-[0 npts(first)*delta]);
    ae=ae(last,:);
end

% merge
data(keep).dep=[data(first).dep; data(last).dep];
data(kill)=[];

% update name, npts, dt
name=name(keep);
npts=sum(npts);
dt=dt(keep,:);
    
end


function [data,ab,ae,npts,dt,name]=...
    mgap(data,ab,ae,delta,npts,dt,option,first,idx,name,diff)
%MGAP    Merge gaps

% how much do we need to shift the samples
shift=cmod(diff,delta);

% get max shift in seconds
if(isequal(option.SHIFTUNITS,'intervals'))
    maxshift=option.SHIFTMAX*delta;
else
    maxshift=option.SHIFTMAX;
end

% who do we keep
last=3-first;
switch option.ADJUST
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
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Gap Method: %s',option.GAP));
end

% how to handle gap?
switch option.GAP
    case 'sequential'
        % just shift the option.ADJUST record to make sequential
        [data,ab,ae,npts,dt,name]=mseq(...
            data,ab,ae,delta,npts,dt,option,first,idx,name);
    case 'interpolate'
        % detail message
        if(option.VERBOSE || option.DEBUG)
            disp(sprintf(...
                ['Interpolate Method: %s\n',...
                 'Adjusting: %s (%d - %s)'...
                 'Adjustment: %f seconds'],...
                option.INTERPOLATE,option.ADJUST,...
                idx(kill),name{kill},shift));
        end
        
        % update
        dt=dt(keep,:);
        name=name(keep);
        
        % shift or interpolate option.ADJUST record to align?
        if(abs(shift)<=maxshift)
            % detail message
            if(option.VERBOSE || option.DEBUG)
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
            if(option.VERBOSE || option.DEBUG)
                disp('Adjust Method: interpolate');
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
        if(option.VERBOSE || option.DEBUG)
            disp(sprintf(...
                ['Filler: %d\n',...
                 'Adjusting: %s (%d - %s)'...
                 'Adjustment: %f seconds'],...
                option.FILLER,option.ADJUST,idx(kill),name{kill},shift));
        end
        
        % update
        dt=dt(keep,:);
        name=name(keep);
        
        % shift or interpolate option.ADJUST record to align?
        if(abs(shift)<=maxshift)
            % detail message
            if(option.VERBOSE || option.DEBUG)
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
            if(option.VERBOSE || option.DEBUG)
                disp('Adjust Method: interpolate');
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
shift=cmod(diff,delta);

% get max shift allowed in seconds (otherwise interpolate)
if(isequal(option.SHIFTUNITS,'intervals'))
    maxshift=option.SHIFTMAX*delta;
else
    maxshift=option.SHIFTMAX;
end

% who do we keep
last=3-first;
switch option.ADJUST
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
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Overlap Method: %s',option.OVERLAP));
end

% how to handle overlap?
switch option.OVERLAP
    case 'sequential'
        % just shift the option.ADJUST record to make sequential
        [data,ab,ae,npts,dt,name]=mseq(...
            data,ab,ae,delta,npts,dt,option,first,idx,name);
    case 'truncate'
        % truncate overlap from option.ADJUST record
        % detail message
        if(option.VERBOSE || option.DEBUG)
            disp(sprintf(...
                ['Adjusting:  %s (%d - %s)\n'...
                 'Adjustment: %f seconds'],...
                option.ADJUST,idx(kill),name{kill},shift));
        end
        
        % update
        dt=dt(keep,:);
        name=name(keep);
        
        % shift or interpolate option.ADJUST record to align?
        if(abs(shift)<=maxshift)
            % detail message
            if(option.VERBOSE || option.DEBUG)
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
            if(option.VERBOSE || option.DEBUG)
                disp('Adjust Method: interpolate');
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
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Added: %d',nsamples));
end

% interpolate gap
gaptimes=(sae(2)+delta):delta:(sae(2)+nsamples*delta);
firsttimes=sab(2):delta:(sab(2)+(npts(first)-1)*delta);
lasttimes=ab(last,2):delta:(ab(last,2)+(npts(last)-1)*delta);
gapdata=interp1([firsttimes lasttimes],...
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
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Added: %d',nsamples));
end

% interpolate gap
gaptimes=(ae(first,2)+delta):delta:(ae(first,2)+nsamples*delta);
firsttimes=ab(first,2):delta:(ab(first,2)+(npts(first)-1)*delta);
lasttimes=sab(2):delta:(sab(2)+(npts(last)-1)*delta);
gapdata=interp1([firsttimes lasttimes],...
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
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Added: %d',nsamples));
end

% interpolate first record and gap
firsttimes=ab(first,2):delta:(ab(first,2)+(npts(first)-1)*delta);
lasttimes=ab(last,2):delta:(ab(last,2)+(npts(last)-1)*delta);
gaptimes=(sae(2)+delta):delta:(sae(2)+nsamples*delta);
newtimes=sab(2):delta:(sab(2)+(npts(first)-1)*delta);
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
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Added: %d',nsamples));
end

% interpolate last record and gap
firsttimes=ab(first,2):delta:(ab(first,2)+(npts(first)-1)*delta);
lasttimes=ab(last,2):delta:(ab(last,2)+(npts(last)-1)*delta);
gaptimes=(ae(first,2)+delta):delta:(ae(first,2)+nsamples*delta);
newtimes=sab(2):delta:(sab(2)+(npts(last)-1)*delta);
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
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Added: %d',nsamples));
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
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Added: %d',nsamples));
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
oldtimes=ab(first,2):delta:(ab(first,2)+(npts(first)-1)*delta);
newtimes=sab(2):delta:(sab(2)+(npts(first)-1)*delta);
data(first).dep=...
    interp1(oldtimes,data(first).dep,newtimes,option.INTERPOLATE,'extrap');
data(first).dep=data(first).dep.';

% figure out number of samples in gap
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);

% detail message
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Added: %d',nsamples));
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
oldtimes=ab(last,2):delta:(ab(last,2)+(npts(last)-1)*delta);
newtimes=sab(2):delta:(sab(2)+(npts(last)-1)*delta);
data(last).dep=...
    interp1(oldtimes,data(last).dep,newtimes,option.INTERPOLATE,'extrap');
data(last).dep=data(last).dep.';

% figure out number of samples in gap
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);

% detail message
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Added: %d',nsamples));
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
%ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% figure out number of samples dropped
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);

% detail message
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Truncated: %d',nsamples));
end

% combine data and drop overlap from first
data(last).dep=[data(first).dep(1:end-nsamples,:); data(last).dep];
data(first)=[];

% set new ab,ae,npts,dt
ab=fixmodserial(sab);
ae=ae(last,:);
npts=npts(1)+npts(2)-nsamples;

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

% detail message
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Truncated: %d',nsamples));
end

% combine data and drop overlap from last
data(first).dep=[data(first).dep; data(last).dep(nsamples+1:end,:)];
data(last)=[];

% set new ab, ae, npts
ab=ab(first,:);
ae=fixmodserial(sae);
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
%ae(last,:)=[ab(first,1) ae(last,2)+86400*(ae(last,1)-ab(first,1))];

% get shifted times of first
sab=ab(first,:)+[0 shift];
sae=ae(first,:)+[0 shift];

% interpolate first record
oldtimes=ab(first,2):delta:(ab(first,2)+(npts(first)-1)*delta);
newtimes=sab(2):delta:(sab(2)+(npts(first)-1)*delta);
data(first).dep=...
    interp1(oldtimes,data(first).dep,newtimes,option.INTERPOLATE,'extrap');
data(first).dep=data(first).dep.';

% figure out number of samples dropped
nsamples=round(abs(ab(last,2)-sae(2)-delta)/delta);

% detail message
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Truncated: %d',nsamples));
end

% combine data and drop overlap from first
data(last).dep=[data(first).dep(1:end-nsamples,:); data(last).dep];
data(first)=[];

% set new ab,ae,npts,dt
ab=fixmodserial(sab);
ae=ae(last,:);
npts=npts(1)+npts(2)-nsamples;

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
oldtimes=ab(last,2):delta:(ab(last,2)+(npts(last)-1)*delta);
newtimes=sab(2):delta:(sab(2)+(npts(last)-1)*delta);
data(last).dep=...
    interp1(oldtimes,data(last).dep,newtimes,option.INTERPOLATE,'extrap');
data(last).dep=data(last).dep.';

% figure out number of samples dropped
nsamples=round(abs(sab(2)-ae(first,2)-delta)/delta);

% detail message
if(option.VERBOSE || option.DEBUG)
    disp(sprintf('Samples Truncated: %d',nsamples));
end

% combine data and drop overlap from last
data(first).dep=[data(first).dep; data(last).dep(nsamples+1:end,:)];
data(last)=[];

% set new ab, ae, npts
ab=ab(first,:);
ae=fixmodserial(sae);
npts=npts(1)+npts(2)-nsamples;

end


function [flags]=flagexactdupes(b,e)
% flags records that start/end at the same time as another as duplicates

% check inputs
if(~isequal(size(b),size(e)))
    error('seizmo:merge:badInputs',['Array sizes do not match:' ...
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
    error('seizmo:merge:badInputs',['Array sizes do not match:' ...
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
dupsu(tril(true(nrecs)))=0; % upper trianrecsle
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
