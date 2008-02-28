function [valid]=vvseis()
%VVSEIS    Returns valid seislab data structure versions
%
%    Description: Returns a numeric vector of valid versions of seismic
%     data files that seislab can work with.
%
%    Notes:
%     - version 6 corresponds to SAC's v6 formatted files
%     - versions 101,102,200,201,202 are modifications of v6.
%    
%    Usage:    valid_versions=vvseis()
%
%    See also:  seischk, isseis, seishi

% no args - return valid versions
valid=[6 101 102 200 201 202];

end