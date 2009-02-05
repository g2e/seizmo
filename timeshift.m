function [data]=timeshift(data,shift,timing,varargin)
%TIMESHIFT    Shift timing of SEIZMO records
%
%    Description: TIMESHIFT(DATA,SHIFT) adjusts the relative timing of
%     SEIZMO records in DATA by SHIFT seconds.  This adjustment is added
%     to all defined header time fields.  The reference time fields are
%     then adjusted by -SHIFT.
%
%     TIMESHIFT(DATA,SHIFT,TIMING) sets the absolute timing to a specified
%     type ('UTC' or 'TAI' - 'TAI' is default).  Handle leapseconds by
%     setting this option to 'UTC'.
%
%     TIMESHIFT(DATA,SHIFT,TIMING,FIELD1,...,FIELDN) also adjusts header
%     fields FIELD1 TO FIELDN by SHIFT seconds.
%
%    Notes:
%     - DOES NOT WORK FOR SPECTRAL OR XYZ RECORDS!
%
%    Tested on: Matlab r2007b
%
%    Header changes: B, E, A, F, O, T0-T9,
%                    NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC
%
%    Usage:    data=timeshift(data,shift)
%              data=timeshift(data,shift,field1,...,fieldN)
%
%    Examples:
%
%    See also: changeheader, getheader

%     Version History:
%        Dec. 13, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec. 13, 2008 at 07:30 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

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
    timing='tai';
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
        times=[nzyear nzjday nzhour nzmin nzsec-shift+nzmsec/1000];
        times=fixtimes(times);
    case 'utc'
        times=[nzyear nzjday nzhour nzmin nzsec-shift+nzmsec/1000];
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
    user{i}=user{i}+shift(:,ones(sz,1));
end
user=[varargin; user];

% change header fields
data=changeheader(data,...
    'nzyear',times(:,1),'nzjday',times(:,2),'nzhour',times(:,3),...
    'nzmin',times(:,4),'nzsec',floor(times(:,5)),...
    'nzmsec',floor(1000*mod(times(:,5),1)),...
    'a',a+shift,'b',b+shift,'e',e+shift,'f',f+shift,'o',o+shift,...
    't',t+shift(:,ones(10,1)),user{:});

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
