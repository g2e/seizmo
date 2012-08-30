function []=seizmo_template(data,varargin)
%SEIZMO_TEMPLATE    Template of SEIZMO functions
%
%    Usage:    []=seizmo_template(data)
%
%    Description:
%     []=SEIZMO_TEMPLATE(DATA) is a start point for a SEIZMO function.
%
%    Notes:
%     - something not mentioned above but still important
%
%    Examples:
%     % Blah:
%     seizmo_template
%
%    See also: SOME_FUNCTION, ANOTHER FUNCTION

%     Version History:
%        Jan.  1, 2009 - initial version
%        Feb. 11, 2011 - mass nargchk fix, mass seizmocheck fix, formatting
%                        updates
%        Jan. 28, 2012 - drop SEIZMO global, some code format changes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2012

% todo:
% - stuff to do

% check nargin
error(nargchk(1,4,nargin));

% check data structure (require data points)
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data,...
        'NONTIME_IFTYPE','ERROR',... % require timeseries
        'FALSE_LEVEN','ERROR');      % require evenly spaced
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% attempt rest
try
    % number of records
    nrecs=numel(data);
    
    % valid values for options
    valid.THING2={'opt1' 'opt2' 'opt3'};
    
    % option defaults
    option.THING1=3; % units, range
    option.THING2='opt1'; % opt1/opt2/opt3
    option.THING3=false; % turn on/off something
    
    % id funtion
    me=mfilename;
    szme=['seizmo:' me ':'];
    
    % get options from command line
    if(~iscellstr(varargin(1:2:end)))
        error([szme 'badInput'],...
            'Options must be specified as a string!');
    end
    for i=1:2:nargin-1
        % assign
        if(~isempty(varargin{i+1}))
            option.(upper(varargin{i}))=varargin{i+1};
        end
        
        % check
        if(~ismember(option.THING2,valid.THING2))
            error([szme 'badInput'],...
                ['THING2 must be of the following:\n',...
                sprintf('%s ',valid.THING2{:})]);
        end
    end
    
    % get header info
    [b,e,delta,npts,ncmp]=getheader(data,'b','e','delta','npts','ncmp');
    
    % loop over records
    [depmin,depmen,depmax]=deal(nan(nrecs,1));
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep)); continue; end
        
        % convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        
        % convert back
        data(i).dep=oclass(data(i).dep);
        
        % get dep*
        depmin=min(data(i).dep(:));
        depmen=nanmean(data(i).dep(:));
        depmax=max(data(i).dep(:));
    end
    
    % update headers
    data=changeheader(data,...
        'depmin',depmin,'depmen',depmen,'depmax',depmax);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end
