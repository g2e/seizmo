function [data]=stack(data,varargin)
%STACK    Simple Time-Domain Stacking of SEIZMO records
%
%    Usage:    data=stack(data)
%              data=stack(data,...)
%
%    Description:
%     DATA=STACK(DATA) will return the stack of records in SEIZMO struct
%     DATA.  Note that the records are essentially passed to ADDRECORDS
%     which performs a sample-by-sample addition without regard to timing.
%     MAKE SURE ALL SAMPLES ARE TIME ALIGNED IN THIS USAGE CASE or your
%     stack likely will not be useful.  You will likely need to call CUT
%     and/or INTERPOLATE first to make this happen.
%
%     DATA=STACK(DATA,...) will pass any additional arguments to STACK on
%     to INTERPOLATE.  This makes STACK call INTERPOLATE internally.  See
%     the Examples section below for additional clarity.  See INTERPOLATE
%     for argument details.  This is the usual usage form.
%
%    Notes:
%
%    Header changes: DEPMIN, DEPMAX, DEPMEN, NPTS, B, E, DELTA, LEVEN
%
%    Examples:
%     % Normalize then stack:
%     data=stack(normalize(data));
%
%     % Return a stack with a delta of 0.2 and a time range of 100 to 200s
%     % making sure that any samples outside the input records data are set
%     % to 0 so they have no influence:
%     data=stack(normalize(data),5,[],100,200,0);
%
%    See also: ADDRECORDS, INTERPOLATE

%     Version History:
%        Mar. 24, 2010 - initial version
%        Dec. 30, 2010 - improved 1-liner description
%        Feb. 11, 2011 - mass nargchk fix
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt interpolation
try
    % interpolate
    if(nargin>1)
        % check options
        for i=1:nargin-1
            if(~isempty(varargin{i}) && ~isscalar(varargin{i}) ...
                    && ~ischar(varargin{i}))
                error('seizmo:stack:badInput',...
                    ['Arguments passed to INTERPOLATE must ' ...
                    'be empty, scalar, or a string!']);
            end
        end
        data=interpolate(data,varargin{:});
    end
    
    % add records
    data=addrecords(data,'ref','ignore');
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
