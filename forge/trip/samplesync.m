function [varargout]=samplesync(varargin)
%SAMPLESYNC    
%
%    Usage:    data=samplesync(data)
%              [data1,...,dataN]=samplesync(data1,...,dataN)
%              [...]=samplesync(...,'adjust',method,...)
%              [...]=samplesync(...,'shiftmax',value,...)
%              [...]=samplesync(...,'shiftunits',units,...)
%              [...]=samplesync(...,'interpolate',method,...)
%              [...]=samplesync(...,'extrapval',value,...)
%              [...]=samplesync(...,'span',method,...)
%              [...]=samplesync(...,'abstiming',logical,...)
%              [...]=samplesync(...,'timestandard',standard,...)
%
%    Description:
%     DATA=SAMPLESYNC(DATA)
%
%     [DATA1,...,DATAN]=SAMPLESYNC(DATA1,...,DATAN)
%
%     [...]=SAMPLESYNC(...,'ADJUST',METHOD,...)
%
%     [...]=SAMPLESYNC(...,'SHIFTMAX',VALUE,...)
%
%     [...]=SAMPLESYNC(...,'SHIFTUNITS',UNITS,...)
%
%     [...]=SAMPLESYNC(...,'INTERPOLATE',METHOD,...)
%
%     [...]=SAMPLESYNC(...,'EXTRAPVAL',VALUE,...)
%
%     [...]=SAMPLESYNC(...,'SPAN',METHOD,...)
%
%     [...]=SAMPLESYNC(...,'ABSTIMING',LOGICAL,...)
%
%     [...]=SAMPLESYNC(...,'TIMESTANDARD',STANDARD,...)
%
%    Notes:
%
%    Examples:
%
%    See also: MELD, ROTATE, ROTATE_CORRELATIONS

%     Version History:
%        Mar. 21, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 21, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% parse pv pairs
[data,pv]=parse_samplesync_inputs(varargin{:});

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

try
    % work based on number of datasets
    if(iscell(data)) % 2+ datasets
        ndata=numel(data);
        tmpdata=data{1}(ones(ndata,1));
        for i=1:numel(data{1})
            for j=1:ndata; tmpdata(j)=data{j}(i); end
            tmpdata=getinstep(tmpdata,pv);
            for j=1:ndata; data{j}(i)=tmpdata(j); end
        end
    else % 1 dataset
        data=getinstep(data,pv);
    end
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% output
varargout=data;

end


function [data]=getinstep(data,pv)

% number of records
nrecs=numel(data);

% get header info
if(pv.ABSTIMING)
    [babs,eabs,npts,delta]=getheader(data,...
        ['b ' pv.TIMESTANDARD],['e ' pv.TIMESTANDARD],'npts','delta');
    babs=cell2mat(babs); eabs=cell2mat(eabs);
    b=timediff(babs(1,:),babs,pv.TIMESTANDARD);
    e=timediff(babs(1,:),eabs,pv.TIMESTANDARD);
else
    [b,e,npts,delta]=getheader(data,'b','e','npts','delta');
end
delta=delta(1);

% use adjust parameter to figure out which record will not be adjusted
switch lower(pv.ADJUST)
    case 'first'
        idx=find(b==max(b),1);
    case 'last'
        idx=find(b==min(b),1);
    case 'longer'
        idx=find(npts==min(npts),1);
    case 'shorter'
        idx=find(npts==max(npts),1);
    case 'one'
        idx=nrecs;
    otherwise % two
        idx=1;
end

% what is the necessary shift?
shift=mod(b(idx)-b,delta);
shift(shift>delta/2)=-delta+shift(shift>delta/2);

% get shift max
switch lower(pv.SHIFTUNITS)
    case 'seconds'
        shiftmax=pv.SHIFTMAX;
    case 'intervals'
        shiftmax=pv.SHIFTMAX*delta;
end

% new times
newb=b+shift;
newe=e+shift;

% pad or truncate?
switch lower(pv.SPAN)
    case 'pad'
        % how many points to add?
        bpad=round((newb-min(newb))/delta);
        epad=round((max(newe)-newe)/delta);
        newb=min(newb); newe=max(newe);
    case 'trim'
        % how many points to remove?
        bcut=round((max(newb)-newb)/delta);
        ecut=round((newe-min(newe))/delta);
        newb=max(newb); newe=min(newe);
end

% handle one record at a time
for i=1:nrecs
    % skip unadjusted
    if(i==idx); continue; end
    
    % interpolate records that need to be adjusted more than shiftmax
    if(shift>shiftmax)
        % interpolating to time-align samples
        data(i).dep=interp1(b(i)+(0:npts(i)-1)*delta,data(i).dep,...
            b(i)+(0:npts(i)-1)*delta+shift(i),pv.INTERPOLATE,'extrap');
    end
    
    % pad or truncate?
    ncmp=size(data(i).dep,2);
    switch lower(pv.SPAN)
        case 'pad'
            data(i).dep=[zeros(bpad(i),ncmp);
                         data(i).dep;
                         zeros(epad(i),ncmp)];
        case 'trim'
            data(i).dep=data(i).dep(1+bcut(i):end-ecut(i),:);
    end
end

% update timing
% - need to fix xabs & npts
if(pv.ABSTIMING)
    data=changeheader(data,['b ' pv.TIMESTANDARD],babs,...
        ['e ' pv.TIMESTANDARD],eabs,'npts',npts);
else
    data=changeheader(data,'b',b,'e',e,'npts',npts);
end

end


function [data,pv]=parse_samplesync_inputs(varargin)

% strip datasets from inputs
data=cell(0,0);
ndata=0;
delete=false(nargin,1);
for i=1:nargin
    if(isseizmo(varargin{i}))
        ndata=ndata+1;
        data{ndata}=varargin{i};
        varargin{i}=[];
        delete(i)=true;
    end
end
varargin(delete)=[];
nargs=numel(varargin);

% require 1+ datasets
if(isempty(data))
    error('seizmo:samplesync:badInput',...
        '1 or more datasets required!');
end

% check dataset sizes
nrecs=numel(data{1});
if(ndata>1)
    for i=2:ndata
        if(numel(data{i})~=nrecs)
            error('seizmo:samplesync:badInput',...
                'Datasets must be equally sized!');
        end
    end
end

% uncell single dataset
if(ndata==1); data=data{1}; end

% require equal samplerates
if(ndata==1) % all records in single dataset
    error(seizmocheck(data,'dep'));
    oldseizmocheckstate=seizmocheck_state(false);
    try
        data=checkheader(data,...
            'MULTIPLE_DELTA','ERROR',...
            'FALSE_LEVEN','ERROR',...
            'NONTIME_IFTYPE','ERROR');
    catch
        % toggle checking back
        seizmocheck_state(oldseizmocheckstate);
        
        % rethrow error
        error(lasterror);
    end
    seizmocheck_state(oldseizmocheckstate);
else % records with same indices
    delta=cell(ndata,1);
    for i=1:ndata
        error(seizmocheck(data{i},'dep'));
        oldseizmocheckstate=seizmocheck_state(false);
        try
            data{i}=checkheader(data{i},...
                'FALSE_LEVEN','ERROR',...
                'NONTIME_IFTYPE','ERROR');
            delta{i}=getheader(data{i},'delta');
        catch
            % toggle checking back
            seizmocheck_state(oldseizmocheckstate);
            
            % rethrow error
            error(lasterror);
        end
        seizmocheck_state(oldseizmocheckstate);
    end
    if(~isequal(delta{:}))
        error('seizmo:samplesync:badInput',...
            'Records to be synced must have equal sample spacing!');
    end
end

% valid values for string options
valid.INTERPOLATE={'spline' 'pchip' 'linear' 'nearest'};
valid.ADJUST={'longer' 'shorter' 'first' 'last' 'one' 'two'};
valid.SHIFTUNITS={'seconds' 'intervals'};
valid.TIMESTANDARD={'utc' 'tai'};
valid.SPAN={'none' 'pad' 'trim'};

% defaults
pv.INTERPOLATE='spline'; % spline/pchip/linear/nearest
pv.ADJUST='shorter'; % longer/shorter/first/last
pv.SHIFTUNITS='intervals'; % seconds/intervals
pv.SHIFTMAX=0.01; % interval: 0-0.5 , seconds: 0+
pv.TIMESTANDARD='utc'; % utc/tai
pv.ABSTIMING=true; % true/false
pv.SPAN='none'; % none/trim/pad
pv.EXTRAP='extrap'; % 'extrap' or some scalar number (nan,0,etc)

% require parameters are specified by strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:samplesync:badInput',...
        'Parameters must be specified with strings!');
end

% parse parameters
for i=1:2:nargs
    switch lower(varargin{i})
        case {'interpolate' 'interp' 'int' 'i'}
            pv.INTERPOLATE=varargin{i+1};
        case {'adjust' 'adj' 'a'}
            pv.ADJUST=varargin{i+1};
        case {'shiftunits' 'shiftu' 'sunit' 'su'}
            pv.SHIFTUNITS=varargin{i+1};
        case {'shiftmax' 'shift' 'smax' 'max' 'sh'}
            pv.SHIFTMAX=varargin{i+1};
        case {'timestandard' 'standard' 'timing' 'time'}
            pv.TIMESTANDARD=varargin{i+1};
        case {'useabsolutetiming' 'useabs' 'abs' 'abstiming'}
            pv.ABSTIMING=varargin{i+1};
        case {'span' 'length'}
            pv.SPAN=varargin{i+1};
        case {'extrapval' 'ex' 'extrap' 'exval' 'extrapolate'}
            pv.EXTRAP=varargin{i+1};
        otherwise
            error('seizmo:samplesync:badInput',...
                'Unknown option: %s',varargin{i});
    end
end

% check parameters
fields=fieldnames(pv);
for i=1:numel(fields)
    % get value of field and do a basic check
    value=pv.(fields{i});
    
    % specific checks
    switch lower(fields{i})
        case 'shiftmax'
            if(~isnumeric(value) || ~isscalar(value))
                error('seizmo:samplesync:badInput',...
                    '%s must be a scalar real number!',fields{i});
            end
        case {'interpolate' 'adjust' 'shiftunits' 'timestandard' 'span'}
            if(~ischar(value) || size(value,1)~=1 ...
                    || ~any(strcmpi(value,valid.(fields{i}))))
                error('seizmo:samplesync:badInput',...
                    ['%s parameter must be one of the following:\n'...
                    sprintf('%s ',valid.(fields{i}){:})],fields{i});
            end
        case 'abstiming'
            if(~islogical(value) || ~isscalar(value))
                error('seizmo:samplesync:badInput',...
                    '%s parameter must be a logical!',fields{i});
            end
        case 'extrap'
            % extrap or a scalar value
            if((ischar(value) && ~strcmpi(value,'extrap')) ...
                    || ~isscalar(value) || ~isnumeric(value))
                error('seizmo:samplesync:badInput',...
                    'EXTRAP parameter must be ''extrap'' or a number!');
            end
        otherwise
            error('seizmo:samplesync:badInput',...
                'Unknown option: %s !',fields{i});
    end
end

end
