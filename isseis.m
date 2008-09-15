function [lgc]=isseis(data,varargin)
%ISSEIS    True for SAClab data structures
%
%    Description: ISSEIS(DATA) returns logical true if DATA is a SAClab 
%     data structure and false otherwise.  See SEISCHK for minimum 
%     requirements.
%
%     ISSEIS(DATA,FIELD1,...,FIELDN) allows adding more fields to the 
%     SAClab data structure requirements - see SEISCHK for details.
%
%    Notes:
%     - ISSEIS is just a logical frontend for SEISCHK
%
%    System requirements: Matlab 7
%
%    Input/Output requirements: first arg can be anything, additional args 
%     must be strings
%
%    Header changes: N/A
%
%    Usage:    logical=isseis(data)
%              logical=isseis(data,field1,...,fieldN)
%
%    Examples:
%     To see if the function rh returns a valid SAClab structure after
%     reading in files from the current directory:
%
%       if(isseis(rh('*'))); disp('Valid Structure'); end
%
%    See also: seischk, seisdef

%     Version History:
%        Feb. 28, 2008 - initial version
%        Mar.  4, 2008 - doc update
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - doc update, input checks
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 14, 2008 at 17:50 GMT

% todo:

% check input
if(nargin<1)
    error('SAClab:isseis:notEnoughInputs','Not enough input arguments.');
elseif(nargin>1)
    for i=1:nargin-1
        if(~ischar(varargin{i}))
            error('SAClab:isseis:badInput',...
                'Additional argument FIELD%d must be a string!',i);
        end
    end
end

% test output of seischk call
lgc=isempty(seischk(data,varargin{:}));

end
