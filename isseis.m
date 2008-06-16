function [lgc]=isseis(data,varargin)
%ISSEIS    True for SAClab data structures
%
%    Description: ISSEIS(DATA) returns logical true if DATA is a SAClab 
%     data structure and false otherwise.  Optional extra inputs are for
%     adding more fields to the seischk input - see seischk for more.
%
%    Usage: logical=isseis(data)
%
%    Examples:
%     To see if the function rh returns a valid SAClab structure after
%     reading in files from the current directory:
%
%       if(isseis(rh('*'))); disp('Valid Structure'); end
%
%    See also: seischk, seisdef

%     Version History:
%        ????????????? - Initial Version
%        June 12, 2008 - Documentation Update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2008 at 16:30 GMT

% test output of seischk call
lgc=isempty(seischk(data,varargin{:}));

end
