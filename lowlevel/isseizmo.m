function [lgc]=isseizmo(data,varargin)
%ISSEIZMO    True for SEIZMO data structures
%
%    Usage:    logical=isseizmo(data)
%              logical=isseizmo(data,field1,...,fieldN)
%
%    Description:
%     ISSEIZMO(DATA) returns logical true if DATA is a SEIZMO data
%     structure and false otherwise.  See SEIZMOCHECK for minimum 
%     requirements.
%
%     ISSEIZMO(DATA,FIELD1,...,FIELDN) allows adding more fields to the 
%     SEIZMO data structure requirements - see SEIZMOCHECK for details.
%
%    Notes:
%     - ISSEIZMO is just a logical frontend for SEIZMOCHECK
%
%    Examples:
%     % To see if the function READHEADER returns a valid SEIZMO structure
%     % after reading in files from the current directory:
%     if(isseizmo(readheader('*'))); disp('Valid Structure'); end
%
%    See also: SEIZMOCHECK, SEIZMODEF

%     Version History:
%        Feb. 28, 2008 - initial version
%        Mar.  4, 2008 - doc update
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - doc update, input checks
%        Nov. 13, 2008 - renamed from ISSEIS to ISSEIZ
%        Nov. 15, 2008 - update for new name schema (now ISSEIZMO)
%        Apr. 23, 2009 - move usage up
%        June 12, 2009 - force SEIZMOCHECK state to on for the check
%        Jan. 29, 2010 - update for new state functions
%        Jan. 28, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2012 at 03:45 GMT

% todo:

% check input
if(nargin<1)
    error('seizmo:isseizmo:notEnoughInputs','Not enough input arguments.');
elseif(nargin>1)
    for i=1:nargin-1
        if(~ischar(varargin{i}))
            error('seizmo:isseizmo:badInput',...
                'Additional argument FIELD%d must be a string!',i);
        end
    end
end

% test output of seizmocheck call
oldseizmocheckstate=seizmocheck_state(true);
lgc=isempty(seizmocheck(data,varargin{:}));
seizmocheck_state(oldseizmocheckstate);

end
