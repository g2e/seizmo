function [data]=timeshift(data,shift,timing,option,varargin)
%TIMESHIFT    Shift timing of SEIZMO records
%
%    Usage:    data=timeshift(data,shift)
%              data=timeshift(data,shift,timing)
%              data=timeshift(data,shift,timing,option)
%              data=timeshift(data,shift,timing,option,field1,...,fieldN)
%
%    Description: TIMESHIFT(DATA,SHIFT) adjusts the relative timing of
%     SEIZMO records in DATA by SHIFT seconds.  This adjustment is added
%     to all defined header time fields (see Header changes section).  The
%     reference time fields are then adjusted by -SHIFT.  This preserves
%     the actual timing of the data and is basically equivalent to SAC's
%     'chnhdr allt shift' command.
%
%     TIMESHIFT(DATA,SHIFT,TIMING) timeshifts while interpreting the
%     absolute timing as the specified type ('UTC' or 'TAI' - 'UTC' is
%     default).  Handle UTC leap seconds by setting this option to 'UTC'.
%     There are no leap seconds in the International Atomic Time standard
%     (TAI), so setting the absolute timing to 'tai' will basically skip
%     any leap second functions.
%
%     TIMESHIFT(DATA,SHIFT,TIMING,OPTION) allows for changing just the
%     reference and relative time fields by setting OPTION to 'REFERENCE'
%     or 'RELATIVE'.  The default is 'BOTH' ([] also works), which changes
%     both.  Note that if the 'REFERENCE' option is given, the shift is
%     applied without a sign change.  All options will apply a shift to
%     additional fields given (see next usage format).  To apply a shift to
%     only the additional fields, use option 'USER'.
%
%     TIMESHIFT(DATA,SHIFT,TIMING,OPTION,FIELD1,...,FIELDN) also adjusts
%     header fields FIELD1 TO FIELDN by SHIFT seconds.
%
%    Notes:
%     - DOES NOT WORK FOR SPECTRAL OR XYZ RECORDS!
%
%    Tested on: Matlab r2007b
%
%    Header changes: B, E, A, F, O, T0-T9,
%                    NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC
%
%    Examples:
%     Shift the reference time to the origin time (note '-' sign):
%      data=timeshift(data,-gh(data,'o'))
%
%     Also useful for quickly plotting data aligned on a phase:
%      plot0(timeshift(data,-Parrivaltimes))
%
%    See also: changeheader, getheader, fixtimes

%     Version History:
%        Dec. 13, 2008 - initial version
%        Mar. 12, 2009 - doc update
%        Mar. 29, 2009 - added OPTION input to allow for more flexibility
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 21:10 GMT

% todo:

% check nargin
msg=nargchk(2,inf,nargin);
if(~isempty(msg)); error(msg); end

% get undefined value
[h,idx]=versioninfo(data);
undef=[h.undef.ntype].';
undef=undef(idx);

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% update header
data=checkheader(data);

% number of records
nrecs=numel(data);

% check shift is numeric
if(~isnumeric(shift))
    error('seizmo:timeshift:badShift',...
        'SHIFT option must be a numeric array!');
elseif(~any(numel(shift)==[1 nrecs]))
    error('seizmo:timeshift:badShift',...
        ['SHIFT option must be a scalar or have '...
        'the same number of elements as DATA!']);
end

% expand scalar shift
if(numel(shift)==1)
    shift(1:nrecs,1)=shift;
else
    shift=shift(:);
end

% default timing
if(nargin==2 || isempty(timing))
    timing='utc';
end

% shift option
if(nargin<=3 || isempty(option))
    refshift=-shift;
    relshift=shift;
    usershift=shift;
else
    if(~ischar(option))
        error('seizmo:timeshift:badOption',...
            ['OPTION option must be string of one of the following:\n'...
             'BOTH  REFERENCE  RELATIVE  USER']);
    end
    switch lower(option)
        case 'both'
            refshift=-shift;
            relshift=shift;
            usershift=shift;
        case 'reference'
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
        otherwise
            error('seizmo:timeshift:badOption',...
                ['OPTION option must be string of one of the following:'...
                 '\nBOTH  REFERENCE  RELATIVE  USER']);
    end
end

% get header fields
nvararg=numel(varargin);
user=cell(1,nvararg);
[a,b,e,f,o,t,nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec,user{:}]=...
    getheader(data,'a','b','e','f','o','t',...
    'nzyear','nzjday','nzhour','nzmin','nzsec','nzmsec',varargin{:});

% only itime, ixy
iftype=getenumid(data,'iftype');
if(any(~(strcmp(iftype,'itime') | strcmp(iftype,'ixy'))))
    error('seizmo:timeshift:badFiletype',...
        'Filetypes for records must be Times Series or General X vs Y !');
end

% get new absolute timing
switch lower(timing)
    case 'tai'
        times=[nzyear nzjday nzhour nzmin nzsec+refshift+nzmsec/1000];
        times=fixtimes(times);
    case 'utc'
        times=[nzyear nzjday nzhour nzmin nzsec+refshift+nzmsec/1000];
        times=fixtimes(times,'utc');
    otherwise
        error('seizmo:timeshift:badTimes',...
            'TIMING option must be ''TAI'' or ''UTC''!');
end

% undefined to NaN
a(a==undef)=nan; b(b==undef)=nan;
e(e==undef)=nan; f(f==undef)=nan;
o(o==undef)=nan; t(t==undef(:,ones(10,1)))=nan;

% deal with user fields
for i=1:nvararg
    if(~isnumeric(user{i}))
        error('seizmo:timeshift:badUserField',...
            'User given fields must return numeric arrays!');
    end
    sz=size(user{i},2);
    user{i}(user{i}==undef(:,ones(sz,1)))=nan;
    user{i}=user{i}+usershift(:,ones(sz,1));
end
user=[varargin; user];

% change header fields
data=changeheader(data,...
    'nzyear',times(:,1),'nzjday',times(:,2),'nzhour',times(:,3),...
    'nzmin',times(:,4),'nzsec',floor(times(:,5)),...
    'nzmsec',floor(1000*mod(times(:,5),1)),...
    'a',a+relshift,'b',b+relshift,'e',e+relshift,'f',f+relshift,...
    'o',o+relshift,'t',t+relshift(:,ones(10,1)),user{:});

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
