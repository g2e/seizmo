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
%    Description: TIMESHIFT(DATA,SHIFT) adjusts the relative timing of
%     SEIZMO records in DATA by SHIFT seconds.  This adjustment is added
%     to all defined header time fields (see Header changes section).  The
%     reference time fields are then adjusted by -SHIFT.  This preserves
%     the actual timing of the data and is basically equivalent to SAC's
%     'chnhdr allt shift' command.
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
%     - DOES NOT WORK FOR SPECTRAL OR XYZ RECORDS!
%     - Since reference time resolution is limited to 1 millisecond, this
%       operation is limited to such a resolution.  Giving a shift with
%       more precision than to the millisecond will not be honored unless
%       OPTION is 'RELATIVE' or 'USER'.
%
%    Header changes: NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC
%                    A, B, E, F, O, Tn, and any user-defined field
%
%    Examples:
%     Shift the reference time to the origin time (note '-' sign):
%      data=timeshift(data,-gh(data,'o'),'io')
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
%        June 24, 2009 - added explaination about inaccurate shifting,
%                        improved checks on options
%        Sep. 25, 2009 - added iztype field parameter, better fixes for
%                        reftime millisec limit
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 25, 2009 at 05:30 GMT

% todo:

% check nargin
msg=nargchk(2,inf,nargin);
if(~isempty(msg)); error(msg); end

% get undefined value
[h,idx]=versioninfo(data);
undef=getsubfield(h,'undef','ntype').';
undef=undef(idx);

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% update header
data=checkheader(data);

% number of records
nrecs=numel(data);

% valid option values
valid.TIMING={'UTC' 'TAI'};
valid.OPTION={'BOTH' 'REFERENCE' 'RELATIVE' 'USER'};

% check shift is numeric
if(~isreal(shift))
    error('seizmo:timeshift:badShift',...
        'SHIFT option must be a numeric array!');
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
times=fixtimes( ...
    [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000+refshift],timing);

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

% combine field and values for changeheader call
user=[varargin; user];

% change header fields
data=changeheader(data,'iztype',iztype,...
    'nzyear',times(:,1),'nzjday',times(:,2),'nzhour',times(:,3),...
    'nzmin',times(:,4),'nzsec',fix(times(:,5)+1e-9),...
    'nzmsec',fix(1000*mod(times(:,5)+1e-9,1)),...
    'a',a+relshift,'b',b+relshift,'e',e+relshift,'f',f+relshift,...
    'o',o+relshift,'t',t+relshift(:,ones(10,1)),user{:});

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
