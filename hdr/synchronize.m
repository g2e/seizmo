function [data]=synchronize(data,field,option,iztype,timing,varargin)
%SYNCHRONIZE    Synchronizes the reference times of SEIZMO records
%
%    Usage:    data=synchronize(data)
%              data=synchronize(data,field)
%              data=synchronize(data,field,option)
%              data=synchronize(data,field,option,iztype)
%              data=synchronize(data,field,option,iztype,timing)
%              data=synchronize(data,field,option,iztype,timing,...
%                               field1,...,fieldN)
%
%    Description:
%     SYNCHRONIZE(DATA) synchronizes the reference times for all records in
%     DATA.  It sets the unified reference time to the starting time of the
%     latest starting record.  All timing fields (A, B, E, F, O, Tn) are
%     correspondingly shifted so that their timing is preserved.  This is
%     particularly useful when combined with CUT to get a time window
%     corresponding to a specific absolute time range.
%
%     SYNCHRONIZE(DATA,FIELD) changes the header field that the
%     synchronization is based on.  FIELD may be any header field that
%     returns one numeric value per record.  By default, FIELD is 'b'.  If
%     FIELD is not a recognized timing field (A, B, E, F, O, Tn), it will
%     not be shifted in the resulting output unless it is also supplied as
%     a user field (see final usage description).  Field may also be 'Z'
%     to align on a reference time.  Field can be a 1x5 or 1x6 numeric
%     vector or a KZDTTM string to directly set the absolute time to
%     synchronize on.
%
%     SYNCHRONIZE(DATA,FIELD,OPTION) sets whether reference times are
%     synced to the first or last occuring FIELD (in absolute time sense).
%     OPTION may be either 'FIRST' or 'LAST' (default).
%
%     SYNCHRONIZE(DATA,FIELD,OPTION,IZTYPE) changes the output records'
%     header field 'iztype' to IZTYPE.  This value is passed directly to
%     CHANGEHEADER as 'changeheader(DATA,'iztype',IZTYPE)'.  The default
%     IZTYPE is 'iunkn'.
%
%     SYNCHRONIZE(DATA,FIELD,OPTION,IZTYPE,TIMING) sets the absolute timing
%     standard to TIMING.  TIMING must be either 'UTC' or 'TAI'.  The
%     default is 'UTC' and takes UTC leap seconds into account.  Using
%     'TAI' assumes that reference times are in TAI time.
%
%     SYNCHRONIZE(DATA,FIELD,OPTION,IZTYPE,TIMING,FIELD1,...,FIELDN) shifts
%     the user-defined timing fields FIELD1 thru FIELDN so their timing is
%     also preserved in the sync.  If a non-standard timing field is used
%     for syncing (FIELD parameter), then that field must be included in
%     this list as well.  See the examples below for further clarity.
%
%    Notes:
%     - Due to the reference time resolution of 1 millisecond, the field
%       being aligned on may not be zero but will be within a millisecond.
%     - Spectral records uses SB instead of B as a timing field!  The E
%       values are temporarily replaced with a valid E timing field and
%       then set back to the frequency counterpart.
%
%    Header changes: NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC,
%                    A, B, E, F, O, Tn, and any user-defined field
%                    SB (if spectral)
%
%    Examples:
%     % Just like SAC:
%     data=synchronize(data);
%
%     % Sync to the first occurance of USER2, making sure to update USERN:
%     data=synchronize(data,'user2','first',[],[],'user');
%
%     % Cut out the first day of 2009 from records:
%     data=cut(synchronize(data,[2009 1 0 0 0]),0,86400);
%
%    See also: TIMESHIFT, CUT

%     Version History:
%        June 24, 2009 - initial version
%        Sep. 24, 2009 - added iztype field parameter
%        Sep. 25, 2009 - better fix for reftime millisec limit, true sync
%        Oct. 17, 2009 - added direct absolute time input
%        Feb.  3, 2010 - proper SEIZMO handling, seizmoverbose support, fix
%                        xy datatype bug, fix checking order bug
%        Aug. 21, 2010 - nargchk fix, better checkheader usage, update
%                        undef checking, drop versioninfo caching
%        Nov.  1, 2011 - doc update
%        Jan. 28, 2012 - doc update, allow date vector input, comments
%        Feb.  5, 2012 - allow spectral records (use SB not B/E)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  5, 2012 at 16:20 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check headers
data=checkheader(data,...
    'XYZ_IFTYPE','ERROR');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt time synchronization
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % detail message
    if(verbose)
        disp('Sychronizing Record(s)');
        print_time_left(0,nrecs);
    end

    % valid option values
    valid.TIMING={'UTC' 'TAI'};
    valid.OPTION={'LAST' 'FIRST'};

    % default/check field (catching warnings)
    abstime=false;
    if(nargin>1 && iscell(field) && isscalar(field)); field=field{1}; end
    if(nargin<2 || isempty(field))
        field='b';
        values=getheader(data,field);
    elseif(isnumeric(field))
        szf=size(field);
        if(szf(1)~=1 || ~any(szf(2)==[2 3 5 6]) ...
                || any(isnan(field) | isinf(field)))
            error('seizmo:synchronize:badField',...
                'FIELD must be 1x2/3/5/6 and not be NaN or Inf!');
        end
        if(any(szf(2)==[2 3]) && ~isequal(field,round(field)))
            error('seizmo:synchronize:badField',...
                'FIELD date vector must be integer!');
        elseif(any(szf(2)==[5 6]) ...
                && ~isequal(field(1:end-1),round(field(1:end-1))))
            error('seizmo:synchronize:badField',...
                'FIELD vector must be integer except for seconds!');
        end
        % expand date w/o time to time at 12am on that day
        if(any(szf(2)==[2 3])); field=[field 0 0 0]; end
        abstime=true;
    elseif(~ischar(field) || size(field,1)~=1)
        error('seizmo:synchronize:badField',...
            'FIELD must be a character string!');
    else
        if(isequal(size(field),[1 29]))
            num=true(1,29);
            num(1,[5 8 11 12 16 17 20 23 26])=false;
            if(~strcmp(field(~num),'-- () ::.') ...
                    || ~all(isstrprop(field(num),'digit')))
                error('seizmo:synchronize:badField',...
                    'FIELD kzdttm-style string improperly formatted!');
            end
            field=str2double([cellstr(field(1:4)) cellstr(field(13:15)) ...
                cellstr(field(18:19)) cellstr(field(21:22)) ...
                cellstr(field(24:25)) cellstr(field(27:29))]);
            abstime=true;
        elseif(strcmpi(field,'z'))
            values=zeros(nrecs,1);
        else
            warning('off','seizmo:getheader:fieldInvalid')
            values=getheader(data,field);
            warning('on','seizmo:getheader:fieldInvalid')
            if(any(isnan(values)) ...
                    || ~isnumeric(values) || size(values,2)>1)
                error('seizmo:synchronize:badOption',...
                    'FIELD must be a valid numeric header field string!');
            end
        end
    end

    % force values to nearest millisecond
    if(~abstime); values=round(values*1000)/1000; end

    % default/check option
    if(nargin<3 || isempty(option))
        option='last';
    elseif(~ischar(option) || size(option,1)~=1 ...
            || ~any(strcmpi(option,valid.OPTION)))
        error('seizmo:synchronize:badOption',...
            ['OPTION must be a string of one of the following:\n'...
            sprintf('%s ',valid.OPTION{:})]);
    end

    % default iztype
    if(nargin<4 || isempty(iztype))
        iztype='iunkn';
    end

    % default/check timing
    if(nargin<5 || isempty(timing))
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
    [iftype,a,b,sb,sdelta,nspts,e,f,o,t,...
        nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec,user{:}]=getheader(data,...
        'iftype id','a','b','sb','sdelta','nspts','e','f','o','t',...
        'nzyear','nzjday','nzhour','nzmin','nzsec','nzmsec',varargin{:});
    se=sb+sdelta.*(nspts-1);
    
    % fix sync field if spectral & b/e is synced on
    spectral=ismember(iftype,{'irlim' 'iamph'});
    if(~abstime && any(spectral))
        switch lower(field)
            case 'b'
                values(spectral)=sb(spectral);
            case 'e'
                values(spectral)=se(spectral);
        end
    end

    % get shift
    reftimes=[nzyear nzjday nzhour nzmin nzsec+nzmsec/1000];
    if(abstime)
        synctime=fixtimes(field,timing);
        % Note: timediff on 'utc' calls fixtimes but 'tai' does not
        %       so we call fixtimes here to assure reftimes is proper
        if(strcmpi(timing,'tai'))
            reftimes=fixtimes(reftimes);
        end
    else
        % sorting b/c we need the first or last occurance
        % - using a datenum conversion here would lose the leapsecond info
        times=sortrows(fixtimes(...
            [reftimes(:,1:4) reftimes(:,5)+values],timing));
        times(:,5)=fix(1000*(times(:,5)+1e-9))/1000; % force msec precision
        switch lower(option)
            case 'last'
                synctime=times(end,:);
            case 'first'
                synctime=times(1,:);
        end
    end

    % get shift
    shift=timediff(synctime,reftimes,timing);

    % shift user fields, but not undefined user fields
    for i=1:nvararg
        if(~isnumeric(user{i}))
            error('seizmo:synchronize:badUserField',...
                'User given fields must return numeric arrays!');
        end
        sz=size(user{i},2);
        user{i}=user{i}+shift(:,ones(sz,1));
    end

    % combine field and values for changeheader call
    user=[varargin; user];
    
    % handle spectral records switch of b/sb, skip e
    time=~spectral;
    b(time)=b(time)+shift(time);
    e(time)=e(time)+shift(time);
    sb(spectral)=sb(spectral)+shift(spectral);

    % change header fields
    data=changeheader(data,'iztype',iztype,...
        'nzyear',synctime(1),'nzjday',synctime(2),...
        'nzhour',synctime(3),'nzmin',synctime(4),...
        'nzsec',fix(synctime(5)+1e-9),...              % trying to handle
        'nzmsec',fix(1000*mod(synctime(5)+1e-9,1)),... % precision issues
        'a',a+shift,'b',b,'sb',sb,'e',e,'f',f+shift,...
        'o',o+shift,'t',t+shift(:,ones(10,1)),user{:});
    
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
