function [data]=timeshift(data,shift,iztype,timing,option,varargin)
%TIMESHIFT    Shift timing of SEIZMO records
%
%    Usage:    data=timeshift(data,shift)
%              data=timeshift(data,shift,iztype)
%              data=timeshift(data,shift,iztype,timing)
%              data=timeshift(data,shift,iztype,timing,option)
%              data=timeshift(data,shift,iztype,timing,option,...
%                             field1,...,fieldN)
%
%    Description:
%     TIMESHIFT(DATA,SHIFT) adjusts the relative & reference timing of
%     SEIZMO records in DATA by SHIFT seconds.  This essentially shifts the
%     reference time of the records while preserving the actual timing of
%     the relative time fields.  SHIFT is added to all defined relative
%     time header fields (see Header changes section) and subtracted from
%     the reference time.  Basically equivalent to the SAC command
%     'chnhdr allt shift'.  Since reference time resolution is limited to
%     milliseconds, this operation is limited to such a resolution.  Giving
%     a shift with more precision than to the millisecond will not be
%     honored unless OPTION (see below) is 'RELATIVE' or 'USER'.
%
%     TIMESHIFT(DATA,SHIFT,IZTYPE) changes the output records' header field
%     'iztype' to IZTYPE.  This value is passed directly to CHANGEHEADER as
%     'changeheader(DATA,'iztype',IZTYPE)'.  The default IZTYPE is 'iunkn'.
%
%     TIMESHIFT(DATA,SHIFT,IZTYPE,TIMING) timeshifts while interpreting the
%     absolute timing as the specified type ('UTC' or 'TAI' - 'UTC' is
%     default).  Handle UTC leap seconds by setting this option to 'UTC'.
%     There are no leap seconds in the International Atomic Time standard
%     (TAI), so setting the absolute timing to 'tai' will basically skip
%     any leap second functions.
%
%     TIMESHIFT(DATA,SHIFT,IZTYPE,TIMING,OPTION) facilitates changing only
%     reference or relative time fields by setting OPTION to 'REFERENCE' or
%     'RELATIVE'.  The default is 'BOTH' ([] also works), which changes
%     both sets.  Note that if OPTION is 'REFERENCE', SHIFT is applied to
%     the reference time without a sign change.  All of the above options
%     will also apply SHIFT to any additional fields given (see next usage
%     format).  To apply a shift to only the additional fields, use option
%     'USER'.
%
%     TIMESHIFT(DATA,SHIFT,IZTYPE,TIMING,OPTION,FIELD1,...,FIELDN) adjusts
%     header fields FIELD1 TO FIELDN by SHIFT seconds.  Giving fields that
%     are already identified as timing fields (A, B, E, F, O, Tn) will
%     apply SHIFT twice unless OPTION is set to 'REFERENCE' or 'USER'!
%
%    Notes:
%     - DOES NOT WORK FOR XYZ RECORDS!
%     - Spectral records update the SB field instead of the B or E field!
%     - Since reference time resolution is limited to milliseconds, this
%       operation is limited to such a resolution.  Giving a shift with
%       more precision than to the millisecond will not be honored unless
%       OPTION is 'RELATIVE' or 'USER'.
%
%    Header changes:
%     Relative time fields:  A, B, E, F, O, Tn and any user-defined field
%     Reference time fields: NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC
%     Spectral records: SB instead of B or E
%
%    Examples:
%     % Shift the reference time to the origin time (note '-' sign):
%     data=timeshift(data,-gh(data,'o'),'io')
%
%     % Also useful for quickly plotting data aligned on a phase:
%     plot0(timeshift(data,-Parrivaltimes))
%
%    See also: CHANGEHEADER, GETHEADER, FIXTIMES, SYNCHRONIZE, CUT, MELD

%     Version History:
%        Dec. 13, 2008 - initial version
%        Mar. 12, 2009 - doc update
%        Mar. 29, 2009 - added OPTION input to allow for more flexibility
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 24, 2009 - added explaination about inaccurate shifting,
%                        improved checks on options
%        Sep. 25, 2009 - added iztype field parameter, better fixes for
%                        reftime millisec limit
%        Feb.  2, 2010 - proper SEIZMO handling, seizmoverbose support,
%                        versioninfo caching, fix bug disallowing xy data
%        Feb. 24, 2010 - require a finite shift
%        Mar.  8, 2010 - drop versioninfo caching (too difficult to debug)
%        Aug. 21, 2010 - nargchk fix, better checkheader usage, update
%                        undef checking
%        Nov.  1, 2011 - doc update
%        Jan. 30, 2012 - minor code comment
%        Feb.  4, 2012 - doc update
%        Feb.  5, 2012 - works with spectral records (updates sb not b/e)
%        Feb.  7, 2012 - merge to meld update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  7, 2012 at 12:40 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% check struct/header
data=checkheader(data,...
    'XYZ_IFTYPE','ERROR');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt timeshift
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);

    % detail message
    if(verbose)
        disp('Time-Shifting Record(s)');
        print_time_left(0,nrecs);
    end
    
    % valid option values
    valid.TIMING={'UTC' 'TAI'};
    valid.OPTION={'BOTH' 'REFERENCE' 'RELATIVE' 'USER'};

    % check shift is numeric
    if(~isreal(shift) || any(isnan(shift) | isinf(shift)))
        error('seizmo:timeshift:badShift',...
            'SHIFT option must be a finite real-valued array!');
    elseif(~any(numel(shift)==[1 nrecs]))
        error('seizmo:timeshift:badShift',...
            ['SHIFT option must be a scalar or have '...
            'the same number of elements as DATA!']);
    end

    % expand scalar shift, make column vector
    if(numel(shift)==1)
        shift(1:nrecs,1)=shift;
    else
        shift=shift(:);
    end

    % default iztype
    if(nargin<3 || isempty(iztype))
        iztype='iunkn';
    end

    % default/check timing
    if(nargin<4 || isempty(timing))
        timing='utc';
    elseif(~ischar(timing) || size(timing,1)~=1 ...
            || ~any(strcmpi(timing,valid.TIMING)))
        error('seizmo:timeshift:badTiming',...
            ['TIMING must be a string of one of the following:\n'...
            sprintf('%s ',valid.TIMING{:})]);
    end

    % default/check/implement shift option
    if(nargin<5 || isempty(option)); option='both'; end
    if(~ischar(option) || size(option,1)~=1 ...
            || ~any(strcmpi(option,valid.OPTION)))
        error('seizmo:timeshift:badOption',...
            ['OPTION must be a string of one of the following:\n'...
            sprintf('%s ',valid.OPTION{:})]);
    end
    switch lower(option)
        case 'both'
            % force shift to nearest millisecond
            shift=round(shift*1000)/1000;
            refshift=-shift;
            relshift=shift;
            usershift=shift;
        case 'reference'
            % force shift to nearest millisecond
            shift=round(shift*1000)/1000;
            refshift=shift;
            relshift=0.*shift;
            usershift=shift;
        case 'relative'
            refshift=0.*shift;
            relshift=shift;
            usershift=shift;
        case 'user'
            refshift=0.*shift;
            relshift=0.*shift;
            usershift=shift;
    end

    % get header fields
    nvararg=numel(varargin);
    user=cell(1,nvararg);
    [iftype,a,b,sb,e,f,o,t,...
        nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec,user{:}]=getheader(data,...
        'iftype id','a','b','sb','e','f','o','t',...
        'nzyear','nzjday','nzhour','nzmin','nzsec','nzmsec',varargin{:});

    % get new absolute timing
    times=fixtimes( ...
        [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000+refshift],timing);

    % deal with user fields
    for i=1:nvararg
        if(~isnumeric(user{i}))
            error('seizmo:timeshift:badUserField',...
                'User given fields must return numeric arrays!');
        end
        sz=size(user{i},2);
        user{i}=user{i}+usershift(:,ones(sz,1));
    end

    % combine field and values for changeheader call
    user=[varargin; user];
    
    % handle spectral records switch of b/sb, skip e
    spectral=ismember(iftype,{'irlim' 'iamph'});
    time=~spectral;
    b(time)=b(time)+relshift(time);
    e(time)=e(time)+relshift(time);
    sb(spectral)=sb(spectral)+relshift(spectral);

    % change header fields
    data=changeheader(data,'iztype',iztype,...
        'nzyear',times(:,1),'nzjday',times(:,2),'nzhour',times(:,3),...
        'nzmin',times(:,4),'nzsec',fix(times(:,5)+1e-9),... % handling
        'nzmsec',fix(1000*mod(times(:,5)+1e-9,1)),...       % precision
        'a',a+relshift,'b',b,'sb',sb,'e',e,'f',f+relshift,...
        'o',o+relshift,'t',t+relshift(:,ones(10,1)),user{:});
    
    % detail message
    if(verbose); print_time_left(nrecs,nrecs); end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
