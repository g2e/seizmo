function [data]=synchronize(data,field,option,timing,varargin)
%SYNCHRONIZE    Synchronizes the reference times of SEIZMO records
%
%    Usage:    data=synchronize(data)
%              data=synchronize(data,field)
%              data=synchronize(data,field,option)
%              data=synchronize(data,field,option,timing)
%              data=synchronize(data,field,option,timing,field1,...,fieldN)
%
%    Description: SYNCHRONIZE(DATA) synchronizes the reference times for
%     all records in DATA.  It sets the unified reference time to the
%     starting time of the latest starting record.  All timing fields (A,
%     B, E, F, O, Tn) are correspondingly shifted so that their timing is
%     preserved.  This is particularly useful when combined with CUT to get
%     a time window corresponding to a specific absolute time range.
%
%     SYNCHRONIZE(DATA,FIELD) allows changing the field for which the
%     records in DATA are synced to FIELD.  FIELD may be any header field
%     which returns one numeric value per record.  By default FIELD is 'b'.
%     If FIELD is not a recognized timing field, it will not be shifted
%     unless it is supplied as a user field (see final usage description).
%     Field may also be 'Z' to align on a reference time.
%
%     SYNCHRONIZE(DATA,FIELD,OPTION) sets whether reference times are
%     synced to the first or last occuring FIELD (in absolute time sense).
%     OPTION may be either 'FIRST' or 'LAST' and is 'LAST' by default.
%
%     SYNCHRONIZE(DATA,FIELD,OPTION,TIMING) sets the absolute timing
%     standard to TIMING.  TIMING must be either 'UTC' or 'TAI'.  The
%     default is 'UTC' and takes UTC leap seconds into account.
%
%     SYNCHRONIZE(DATA,FIELD,OPTION,TIMING,FIELD1,...,FIELDN) also shifts
%     the user-defined timing fields FIELD1 through FIELDN so that their
%     timing is preserved in the sync.  If a non-standard timing field is
%     used for syncing, then that field must be included in this list too.
%     See the examples below for further clarity.
%
%    Notes:
%     - Due to the reference time resolution of 1 millisecond, the field
%       being aligned on may not be zero but will be less than 1
%       millisecond.  Timing is preserved.
%
%    Header changes: NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC,
%                    A, B, E, F, O, Tn, and any user-defined field.
%
%    Examples:
%     Just like SAC:
%      data=synchronize(data);
%
%     Sync to the first occurance of USER2, making sure to update USERN:
%      data=synchronize(data,'user2','first',[],'user');
%
%    See also: timeshift, cut

%     Version History:
%        June 24, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:30 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% get undefined value
[h,idx]=versioninfo(data);
undef=[h.undef.ntype].';
undef=undef(idx);

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% number of records
nrecs=numel(data);

% valid option values
valid.TIMING={'UTC' 'TAI'};
valid.OPTION={'LAST' 'FIRST'};

% check field (catching warnings)
if(nargin<2 || isempty(field))
    field='b';
    values=getheader(data,field);
elseif(~ischar(field) || size(field,1)~=1)
    error('seizmo:synchronize:badField',...
        'FIELD must be a character string!');
else
    if(strcmpi(field,'z'))
        values=zeros(nrecs,1);
    else
        warning('off','seizmo:getheader:fieldInvalid')
        values=getheader(data,field);
        warning('on','seizmo:getheader:fieldInvalid')
        if(isnan(values) || ~isnumeric(values) || size(values,2)>1)
            error('seizmo:synchronize:badOption',...
                'FIELD must be a valid numeric header field string!');
        end
    end
end

% check option
if(nargin<3 || isempty(option))
    option='last';
elseif(~ischar(option) || size(option,1)~=1 ...
        || ~any(strcmpi(option,valid.OPTION)))
    error('seizmo:synchronize:badOption',...
        ['OPTION must be a string of one of the following:\n'...
         sprintf('%s ',valid.OPTION{:})]);
end

% check timing
if(nargin<4 || isempty(timing))
    timing='utc';
elseif(~ischar(timing) || size(timing,1)~=1 ...
        || ~any(strcmpi(timing,valid.TIMING)))
    error('seizmo:synchronize:badTiming',...
        ['TIMING must be a string of one of the following:\n'...
        sprintf('%s ',valid.TIMING{:})]);
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
    error('seizmo:synchronize:badFiletype',...
        'Filetypes for records must be Times Series or General X vs Y !');
end

% get shift
reftimes=[nzyear nzjday nzhour nzmin nzsec+nzmsec/1000];
times=sortrows(fixtimes([reftimes(:,1:4) reftimes(:,5)+values],timing));
times(:,5)=round(1000*times(:,5))/1000; % account for millisecond limit
switch lower(option)
    case 'last'
        synctime=times(end,:);
    case 'first'
        synctime=times(1,:);
end

% get shift
shift=timediff(synctime,reftimes,timing);

% shift reftimes
reftimes=fixtimes([reftimes(:,1:4) reftimes(:,5)-shift],timing);

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
    'nzyear',reftimes(:,1),'nzjday',reftimes(:,2),...
    'nzhour',reftimes(:,3),'nzmin',reftimes(:,4),...
    'nzsec',floor(reftimes(:,5)),...
    'nzmsec',round(1000*mod(reftimes(:,5),1)),...
    'a',a+shift,'b',b+shift,'e',e+shift,'f',f+shift,...
    'o',o+shift,'t',t+shift(:,ones(10,1)),user{:});

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
