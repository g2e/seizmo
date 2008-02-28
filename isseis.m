function [lgc]=isseis(data,varargin)
%ISSEIS    True for seislab data structures
%
%    Description: ISSEIS(DATA) returns true if DATA is a seislab data
%     structure and false otherwise.  Optional extra inputs are for adding
%     more fields to the seischk input - see seischk for more.
%
%    Examples:
%     To see if the function rh returns a valid seislab structure after
%     reading in files from the current directory:
%
%       files=dir;
%       isseis(rh(files.name))
%
%    See also: seischk, seishi

% test output of seischk call
lgc=isempty(seischk(data,varargin{:}));

end
