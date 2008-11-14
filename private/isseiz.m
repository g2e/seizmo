function [lgc]=isseiz(data,varargin)
%ISSEIZ    True for SAClab data structures
%
%    Description: ISSEIZ(DATA) returns logical true if DATA is a SAClab 
%     data structure and false otherwise.  See SEIZCHK for minimum 
%     requirements.
%
%     ISSEIZ(DATA,FIELD1,...,FIELDN) allows adding more fields to the 
%     SAClab data structure requirements - see SEIZCHK for details.
%
%    Notes:
%     - ISSEIZ is just a logical frontend for SEIZCHK
%
%    Tested on: Matlab r2007b
%
%    Usage:    logical=isseiz(data)
%              logical=isseiz(data,field1,...,fieldN)
%
%    Examples:
%     To see if the function rh returns a valid SAClab structure after
%     reading in files from the current directory:
%
%       if(isseiz(rh('*'))); disp('Valid Structure'); end
%
%    See also: seizchk, seizdef

%     Version History:
%        Feb. 28, 2008 - initial version
%        Mar.  4, 2008 - doc update
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - doc update, input checks
%        Nov. 13, 2008 - renamed from ISSEIS to ISSEIZ
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 13, 2008 at 04:50 GMT

% todo:

% check input
if(nargin<1)
    error('SAClab:isseiz:notEnoughInputs','Not enough input arguments.');
elseif(nargin>1)
    for i=1:nargin-1
        if(~ischar(varargin{i}))
            error('SAClab:isseiz:badInput',...
                'Additional argument FIELD%d must be a string!',i);
        end
    end
end

% test output of seischk call
lgc=isempty(seizchk(data,varargin{:}));

end
