function [data]=mirrorflip(data,varargin)
%MIRRORFLIP    Returns a mirror-flip of SEIZMO records
%
%    Usage:    data=mirrorflip(data)
%              data=mirrorflip(data,...,'type',type,...)
%              data=mirrorflip(data,...,'length',N,...)
%
%    Description:
%     MIRRORFLIP(DATA) returns a mirror-flipped version of SEIZMO records
%     in DATA prepended to the original records.  The output is thus
%     2*NPTS-1 points long, where NPTS is the number of points in the
%     original record.  The begin time is shifted forward to B-(E-B), where
%     B is the original begin time and E is the end time of the record.
%
%     MIRRORFLIP(DATA,...,'TYPE',TYPE,...) allows changing the direction
%     and output.  Valid TYPE strings are the following:
%      'front'      - returns just the front-sided mirror-flip 
%      'back'       - returns just the back-sided mirror-flip
%      'prepend'    - returns the front-sided mirror-flip prepended to the
%                     original record (DEFAULT)
%      'append'     - returns the back-sided mirror-flip appended to the
%                     original record
%      'both'       - returns both the back and front-sided mirror-flips
%                     attached to the original record
%
%     MIRRORFLIP(DATA,...,'LENGTH',N,...) changes the number of samples
%     that are mirror-flipped.  By default N is NPTS-1 (the maximum), where
%     NPTS is the number of points in each record.  Setting N to NaN will
%     force the default behavior.
%
%    Notes:
%     - Mirror-flipping is often useful for reducing strong edge-effects
%       in filter-related processing by faking data ahead of the record
%       window being filtered on.  This is known in signal processing as
%       pseudo-initial conditioning the filter.
%
%    Header changes: B, E, NPTS, DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Use mirror-flip to dodge edge-effects when decimating:
%     [b,e]=getheader(data,'b','e');
%     data=mirrorflip(data,'type','both');
%     data=squish(data,[2 2 2 5]);
%     data=cut(data,b,e);
%
%    See also: REVERSE, IIRFILTER

%     Version History:
%        May  14, 2009 - initial version
%        May  28, 2009 - bug fix: had to force order of operations due to
%                        concatination/arithmitic ambiguity
%        May  29, 2009 - make nan the default for length option
%        June  3, 2009 - fixes for options checking
%        Sep. 11, 2009 - minor doc update
%        Jan. 30, 2010 - proper SEIZMO handling, seizmoverbose support
%        Feb. 11, 2011 - mass seizmocheck fix
%        Jan. 28, 2012 - doc update, drop SEIZMO global, better checkheader
%                        usage
%        Mar. 13, 2012 - use getheader improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 15:05 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:mirrorflip:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
error(seizmocheck(data,'dep'));

% toggle off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt mirror-flip
try
    % check headers
    data=checkheader(data,'NONTIME_IFTYPE','ERROR');
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);

    % valid options
    valid.TYPE={'front' 'back' 'prepend' 'append' 'both'};

    % default options
    option.TYPE='prepend';
    option.LENGTH=nan;

    % get options from command line
    for i=1:2:nargin-1
        if(~ischar(varargin{i}))
            error('seizmo:mirrorflip:badInput',...
                'Options must be specified as strings!');
        end
        if(~isempty(varargin{i+1}))
            option.(upper(varargin{i}))=varargin{i+1};
        end
    end

    % check options
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

                if(~isempty(option.(fields{i})) ...
                        && (~isequalwithequalnans( ...
                        fix(option.(fields{i})),option.(fields{i})) ...
                        || ~any(numel(option.(fields{i}))==[1 nrecs])))
                    error('seizmo:mirrorflip:badInput',...
                        ['%s option must be an integer or an array\n'...
                        'of integers with one option per record!'],...
                        fields{i});
                end
                if(isscalar(option.(fields{i})))
                    option.(fields{i})=option.(fields{i})(ones(nrecs,1),1);
                end
        end
    end

    % get header fields
    [b,e,delta,npts,leven]=getheader(data,...
        'b','e','delta','npts','leven lgc');
    leven=~strcmpi(leven,'false');

    % fix empty length option
    if(isnan(option.LENGTH))
        option.LENGTH=npts-1;
    end

    % force length option to be no larger than npts-1
    option.LENGTH=min(option.LENGTH,npts-1);
    
    % detail message
    if(verbose)
        disp('Mirror-Flipping Record(s)');
        print_time_left(0,nrecs);
    end

    % work on records individually
    [depmen,depmin,depmax]=deal(nan(nrecs,1));
    for i=1:nrecs
        % dataless support
        if(~npts(i))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);

        % mirror-flip
        switch lower(option.TYPE{i})
            case 'front'
                data(i).dep=2*data(i).dep(ones(option.LENGTH(i),1),:) ...
                    - data(i).dep(option.LENGTH(i)+1:-1:2,:);
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
                data(i).dep=...
                    2*data(i).dep(end*ones(option.LENGTH(i),1),:) ...
                    - data(i).dep(end-1:-1:end-option.LENGTH(i),:);
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
                    - data(i).dep(option.LENGTH(i)+1:-1:2,:)); ...
                    data(i).dep];
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
                    - data(i).dep(end-1:-1:end-option.LENGTH(i),:))];
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
                    - data(i).dep(option.LENGTH(i)+1:-1:2,:)); ...
                    data(i).dep; ...
                    (2*data(i).dep(end*ones(option.LENGTH(i),1),:) ...
                    - data(i).dep(end-1:-1:end-option.LENGTH(i),:))];
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
        depmen(i)=nanmean(data(i).dep(:));
        if(~isempty(data(i).dep))
            depmin(i)=min(data(i).dep(:));
            depmax(i)=max(data(i).dep(:));
        else
            depmin(i)=nan;
            depmax(i)=nan;
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update headers
    data=changeheader(data,'npts',npts,'b',b,'e',e,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
