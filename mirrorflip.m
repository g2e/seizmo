function [data]=mirrorflip(data,varargin)
%MIRRORFLIP    Returns a mirror-flip of SEIZMO records
%
%    Usage:    data=mirrorflip(data)
%              data=mirrorflip(data,...,'type',type)
%              data=mirrorflip(data,...,'length',N)
%
%    Description: MIRRORFLIP(DATA) returns a mirror-flipped version of
%     SEIZMO records in DATA prepended to the original records.  The output
%     is thus 2*NPTS-1 points long, where NPTS is the number of points in
%     the original record.  The begin time is shifted forward to B-(E-B),
%     where B is the original begin time and E is the end time of the
%     record.
%
%     MIRRORFLIP(DATA,...,'TYPE',TYPE) allows changing the direction and
%     output.  Valid TYPE strings are the following:
%      'front'      - returns just the front-sided mirror-flip 
%      'back'       - returns just the back-sided mirror-flip
%      'prepend'    - returns the front-sided mirror-flip prepended to the
%                     original record (DEFAULT)
%      'append'     - returns the back-sided mirror-flip appended to the
%                     original record
%      'both'       - returns both the back and front-sided mirror-flips
%                     attached to the original record
%
%     MIRRORFLIP(DATA,...,'LENGTH',N) changes the number of samples that
%     are mirror-flipped.  By default N is NPTS-1 (the maximum), where NPTS
%     is the number of points in each record.  Setting N to NaN will force
%     the default behavior.
%
%    Notes:
%     - Mirror-flipping is often useful for avoiding strong edge-effects
%       in filter-related processing.  This only limits such effects on the
%       interior portion of data.  Data near the edge should not be trusted
%       in any case.
%
%    Header changes: B, E, NPTS, DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Use mirror-flip to dodge edge-effects when decimating:
%      [b,e]=getheader(data,'b','e');
%      data=mirrorflip(data,'type','both');
%      data=squish(data,[2 2 2 5]);
%      data=cut(data,b,e);
%
%    See also: reverse

%     Version History:
%        May  14, 2009 - initial version
%        May  28, 2009 - bug fix: had to force order of operations due to
%                        concatination/arithmitic ambiguity
%        May  29, 2009 - make nan the default for length option
%        June  3, 2009 - fixes for options checking
%        June 12, 2009 - add testing table
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  3, 2009 at 16:20 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:mirrorflip:badNumInputs',...
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

% valid options
valid.TYPE={'front' 'back' 'prepend' 'append' 'both'};

% default options
option.TYPE='prepend';
option.LENGTH=nan;

% get options from SEIZMO global
global SEIZMO
try
    fields=fieldnames(SEIZMO.MIRRORFLIP);
    for i=1:numel(fields)
        if(~isempty(SEIZMO.MIRRORFLIP.(fields{i})))
            option.(fields{i})=SEIZMO.MIRRORFLIP.(fields{i});
        end
    end
catch
end

% get options from command line
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error('seizmo:mirrorflip:badInput',...
            'Options must be specified as a strings!');
    end
    if(~isempty(varargin{i+1}))
        option.(upper(varargin{i}))=varargin{i+1};
    end
end

% check options
nrecs=numel(data);
fields=fieldnames(option);
for i=1:numel(fields)
    % specific checks
    switch lower(fields{i})
        case 'type'
            if(iscellstr(option.(fields{i})))
                option.(fields{i})=char(option.(fields{i}));
            end
            if(isempty(option.(fields{i})) ...
                    || ~ischar(option.(fields{i})) ...
                    || ~any(size(option.(fields{i}),1)==[1 nrecs]) ...
                    || ~isempty(setdiff(lower(option.(fields{i})),...
                    valid.(fields{i}))))
                error('seizmo:mirrorflip:badInput',...
                    ['%s option must be one of the following:\n'...
                    sprintf('%s ',valid.(fields{i}){:})],fields{i});
            end
            if(size(option.(fields{i}),1)==1)
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),:);
            end
            option.(fields{i})=cellstr(option.(fields{i}));
        case 'length'
            if(~isempty(option.(fields{i})) && (...
                    any(fix(option.(fields{i}))~=option.(fields{i})) ...
                    || ~any(numel(option.(fields{i}))==[1 nrecs])))
                error('seizmo:mirrorflip:badInput',...
                    ['%s option must be an integer or an array\n'...
                    'of integers with one option per record.'],fields{i});
            end
            if(isscalar(option.(fields{i})))
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),1);
            end
    end
end

% get header fields
iftype=getenumid(data,'iftype');
leven=strcmp(getlgc(data,'leven'),'true');
[b,e,delta,npts]=getheader(data,'b','e','delta','npts');

% fix empty length option
if(isnan(option.LENGTH))
    option.LENGTH=npts-1;
end

% force length option to be no larger than npts-1
option.LENGTH=min(option.LENGTH,npts-1);

% require timeseries and general x vs y
if(any(~strcmp(iftype,'itime') & ~strcmp(iftype,'ixy')))
    bad=find(~strcmp(iftype,'itime') & ~strcmp(iftype,'ixy'));
    error('seizmo:mirrorflip:badRecordType',...
        ['Records:' sprintf('%d ',bad) '\n'...
        'Illegal operation on spectral/xyz record!']);
end

% work on records individually
[depmen,depmin,depmax]=deal(nan(nrecs,1));
for i=1:nrecs
    % dataless support
    if(~npts(i)); continue; end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % mirror-flip
    switch lower(option.TYPE{i})
        case 'front'
            data(i).dep=2*data(i).dep(ones(option.LENGTH(i),1),:) ...
                -data(i).dep(option.LENGTH(i)+1:-1:2,:);
            % handle special case (single point => dataless)
            if(isempty(data(i).dep))
                npts(i)=0;
                b(i)=nan;
                e(i)=nan;
                if(~leven(i))
                    data(i).ind=[];
                end
            else
                npts(i)=option.LENGTH(i);
                if(leven(i))
                    newb=b(i)-delta(i)*option.LENGTH(i);
                    e(i)=b(i)-delta(i);
                else
                    newb=2*b(i)-data(i).ind(option.LENGTH(i)+1);
                    e(i)=2*b(i)-data(i).ind(2);
                    data(i).ind=2*data(i).ind(1) ...
                        -data(i).ind(option.LENGTH(i)+1:-1:2);
                end
                b(i)=newb;
            end
        case 'back'
            data(i).dep=2*data(i).dep(end*ones(option.LENGTH(i),1),:) ...
                -data(i).dep(end-1:-1:end-option.LENGTH(i),:);
            % handle special case (single point => dataless)
            if(isempty(data(i).dep))
                npts(i)=0;
                b(i)=nan;
                e(i)=nan;
                if(~leven(i))
                    data(i).ind=[];
                end
            else
                npts(i)=option.LENGTH(i);
                if(leven(i))
                    newe=e(i)+delta(i)*option.LENGTH(i);
                    b(i)=e(i)+delta(i);
                else
                    newe=2*e(i)-data(i).ind(end-option.LENGTH(i));
                    b(i)=2*e(i)-data(i).ind(end-1);
                    data(i).ind=2*data(i).ind(end) ...
                        -data(i).ind(end-1:-1:end-option.LENGTH(i));
                end
                e(i)=newe;
            end
        case 'prepend'
            data(i).dep=[(2*data(i).dep(ones(option.LENGTH(i),1),:) ...
                -data(i).dep(option.LENGTH(i)+1:-1:2,:)); data(i).dep];
            npts(i)=npts(i)+option.LENGTH(i);
            if(leven(i))
                b(i)=b(i)-delta(i)*option.LENGTH(i);
            else
                b(i)=2*b(i)-data(i).ind(option.LENGTH(i)+1);
                data(i).ind=[(2*data(i).ind(1) ...
                    -data(i).ind(option.LENGTH(i)+1:-1:2)); ...
                    data(i).ind];
            end
        case 'append'
            data(i).dep=[data(i).dep; ...
                (2*data(i).dep(end*ones(option.LENGTH(i),1),:) ...
                -data(i).dep(end-1:-1:end-option.LENGTH(i),:))];
            npts(i)=npts(i)+option.LENGTH(i);
            if(leven(i))
                e(i)=e(i)+delta(i)*option.LENGTH(i);
            else
                e(i)=2*e(i)-data(i).ind(end-option.LENGTH(i));
                data(i).ind=[data(i).ind; ...
                    (2*data(i).ind(end) ...
                    -data(i).ind(end-1:-1:end-option.LENGTH(i)))];
            end
        case 'both'
            data(i).dep=[(2*data(i).dep(ones(option.LENGTH(i),1),:) ...
                -data(i).dep(option.LENGTH(i)+1:-1:2,:)); data(i).dep; ...
                (2*data(i).dep(end*ones(option.LENGTH(i),1),:) ...
                -data(i).dep(end-1:-1:end-option.LENGTH(i),:))];
            npts(i)=npts(i)+2*option.LENGTH(i);
            if(leven(i))
                b(i)=b(i)-delta(i)*option.LENGTH(i);
                e(i)=e(i)+delta(i)*option.LENGTH(i);
            else
                b(i)=2*b(i)-data(i).ind(option.LENGTH(i)+1);
                e(i)=2*e(i)-data(i).ind(end-option.LENGTH(i));
                data(i).ind=[(2*data(i).ind(1) ...
                    -data(i).ind(option.LENGTH(i)+1:-1:2)); ...
                    data(i).ind; ...
                    (2*data(i).ind(end) ...
                    -data(i).ind(end-1:-1:end-option.LENGTH(i)))];
            end
    end
    
    % change class back
    data(i).dep=oclass(data(i).dep);
    
    % get dep* info
    depmen(i)=mean(data(i).dep(:));
    if(~isempty(data(i).dep))
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
    else
        depmin(i)=nan;
        depmax(i)=nan;
    end
end

% update headers
data=changeheader(data,'depmen',depmen,'depmin',depmin,'depmax',depmax,...
    'npts',npts,'b',b,'e',e);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
