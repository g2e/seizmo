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
%    See also: SOME_FUNCTION

%     Version History:
%        Jan.  1, 2009 - initial version
%        Feb. 11, 2011 - mass nargchk fix, mass seizmocheck fix, formatting
%                        updates
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

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
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% get global
global SEIZMO

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
    
    % get options from SEIZMO global
    ME=upper(me);
    try
        fields=fieldnames(SEIZMO.(ME));
        for i=1:numel(fields)
            if(~isempty(SEIZMO.(ME).(fields{i})))
                option.(fields{i})=SEIZMO.(ME).(fields{i});
            end
        end
    catch
    end
    
    % get options from command line
    for i=1:2:nargin-1
        if(~ischar(varargin{i}))
            error([szme 'badInput'],...
                'Options must be specified as a string!');
        end
        if(~isempty(varargin{i+1}))
            option.(upper(varargin{i}))=varargin{i+1};
        end
    end
    
    % get header info
    [b,e,delta,npts,ncmp]=getheader(data,'b','e','delta','npts','ncmp');
    
    % loop over records
    depmin=nan(nrecs,1); depmen=depmin; depmax=depmin;
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
        depmen=mean(data(i).dep(:));
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
