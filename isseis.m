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
%       files=dir;
%       if(isseis(rh(files.name))); disp('Valid Structure'); end
%
%    See also: seischk, seishi

% test output of seischk call
lgc=isempty(seischk(data,varargin{:}));

end
