function [data]=convolve(data,varargin)
%CONVOLVE    Convolve SEIZMO records with a function
%
% COMING SOON!!!

% todo:
% - need options to create convolve function
%   - type (all windows?)
%     - triangle
%     - ramp
%     - square
%     - gaussian (other windows?)
%     - custom
%       - gauspuls
%       - chirp
%       - sinc
%       - sawtooth
%   - rampfraction (for ramp type) (0 to 0.5)
%   - nstddev (for gaussian type)
%   - halfduration
%   - maxamplitude
%   - offset
%
% - notes: 
%   - CONV is causal - need to offset by halfduration (just timing)
%   - study specfem3d - how do they work with kinematic sources
%   - 

% check nargin
if(mod(nargin-1,2))
    error('seizmo:convolve:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% toggle off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% valid option values and ranges
valid.TYPE={'triangle' 'ramp' 'square' 'gaussian' 'custom'};
valid.SKEW=[-1 1];
valid.RAMPFRACTION=[0 1];
valid.HALFDURATION=[0 inf];

% default options
option.TYPE='gaussian';
option.HALFDURATION=10;
option.MAXAMPLITUDE=1;
option.NSTDDEV=2.5;
option.SKEW=0;
option.RAMPFRACTION=0.25;
option.OFFSET=0;
option.CUSTOM=[];

% get options from SEIZMO global
global SEIZMO
try
    fields=fieldnames(SEIZMO.CONVOLVE);
    for i=1:numel(fields)
        if(~isempty(SEIZMO.CONVOLVE.(fields{i})))
            option.(fields{i})=SEIZMO.CONVOLVE.(fields{i});
        end
    end
catch
end

% get options from command line
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error('seizmo:convolve:badInput',...
            'Options must be specified as a strings!');
    end
    if(~isempty(varargin{i+1}))
        option.(upper(varargin{i}))=varargin{i+1};
    end
end

% check options
fields=lower(fieldnames(option));
for i=1:numel(fields)
    % get value of field and do a basic check
    value=option.(fields{i});
    
    % specific checks
    switch lower(fields{i})
        case 'type'
            if(iscellstr(value))
                value=char(value);
                option.(fields{i})=char(option.(fields{i}));
            end
            if(isempty(value) || ~ischar(value) ...
                    || ~any(size(value,1)==[1 nrecs]) ...
                    || ~isempty(setdiff(value,valid.(fields{i}))))
                error('seizmo:convolve:badInput',...
                    ['%s option must be one of the following:\n'...
                    sprintf('%s ',valid.(fields{i}){:})],fields{i});
            end
            if(size(value,1)==1)
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),:);
            end
        case {'maxamplitude' 'offset' 'nstddev'}
            if(isempty(value) || ~isreal(value) ...
                    || ~any(numel(value)==[1 nrecs]))
                error('seizmo:convolve:badInput',...
                    ['%s option must be a real number or \n'...
                    'a real array with one option per record.'],fields{i});
            end
            if(isscalar(value))
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),1);
            end
        case {'halfduration' 'skew'}
            if(isempty(value) || ~isreal(value) ...
                    || any(value<valid.(fields{i})(1) ...
                    | value>valid.(fields{i})(2)) ...
                    || ~any(numel(value)==[1 nrecs]))
                error('seizmo:convolve:badInput',...
                    ['%s option must be a real number within %d and %d\n'...
                    'or a real array with one option per record,\n'...
                    'all within %d and %d.'],fields{i},...
                    valid.(fields{i})(1),valid.(fields{i})(2),...
                    valid.(fields{i})(1),valid.(fields{i})(2));
            end
            if(isscalar(value))
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),1);
            end
        case 'rampfraction'
            sz=size(value);
            if(isempty(value) || ~isreal(value) ...
                    || any(value<valid.(fields{i})(1) ...
                    | value>valid.(fields{i})(2)) ...
                    || numel(sz)>2 || ~any(sz(2)==1:2) ...
                    || ~any(sz(1)==[1 nrecs]))
                error('seizmo:convolve:badInput',...
                    ['%s option must be a real number within %d and %d\n'...
                    'or a real array with one option per record,\n'...
                    'all within %d and %d.'],fields{i},...
                    valid.(fields{i})(1),valid.(fields{i})(2),...
                    valid.(fields{i})(1),valid.(fields{i})(2));
            end
            if(sz(1)==1)
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),:);
            end
            if(sz(2)==1)
                option.(fields{i})=[option.(fields{i}) option.(fields{i})];
            end
        case 'custom'
            csz=size(value);
            if(~isempty(value) && (~isreal(value) ...
                    || numel(csz)>2 || csz(1)~=2))
                error('seizmo:convolve:badInput',...
                    '%s option must be a 2xN real array!',fields{i});
            end
        otherwise
            error('seizmo:convolve:badInput',...
                'Unknown option: %s !',fields{i});
    end
end

% get header fields
iftype=getenumid(data,'iftype');
leven=getlgc(data,'leven');
[delta,npts,e]=getheader(data,'delta','npts','e');

% require timeseries and general x vs y
if(any(~strcmp(iftype,'itime') & ~strcmp(iftype,'ixy')))
    bad=find(~strcmp(iftype,'itime') & ~strcmp(iftype,'ixy'));
    error('seizmo:convolve:badRecordType',...
        ['Records:' sprintf('%d ',bad) '\n'...
        'Illegal operation on spectral/xyz record!']);
end

% require evenly spaced
if(any(~strcmpi(leven,'true')))
    error('seizmo:convolve:badRecordSpacing',...
        'Illegal operation on unevenly spaced record!')
end

% convolve and update
[depmen,depmin,depmax]=deal(nan(nrecs,1));
for i=1:nrecs
    % dataless support
    if(~npts(i)); continue; end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % build convolution function
    % - all non-custom options are symmetric
    last=fix(option.HALFDURATION(i)/delta(i));
    switch lower(option.TYPE)
        case 'gaussian'
            
        case 'ramp'
            
        case 'triangle'
            tripuls(-last:delta(i):last,option.HALFDURATION(i),option.SKEW(i));
        case 'square'
            stf=option.MAXAMPLITUDE(i).*...
                ones(1,2*fix(option.HALFDURATION(i)/delta(i))+1);
        case 'custom'
            option.HALFDURATION(i)=0;
            option.OFFSET(i)=0;
    end
    
    % convolve
    data(i).dep=conv(data(i).dep,stf);
    
    % change class back
    data(i).dep=oclass(data(i).dep);
    
    % update npts, e
    addnpts=numel(stf)-1;
    npts(i)=npts(i)+addnpts;
    e(i)=e(i)+addnpts*delta(i);
    
    % get dep* info
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
end

% update headers
data=changeheader(data,'depmen',depmen,'depmin',depmin,'depmax',depmax,...
    'npts',npts,'e',e);

% fix timing
data=timeshift(data,-option.HALFDURATION+option.OFFSET);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end