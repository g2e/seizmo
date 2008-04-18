function [valid]=vvseis()
%VVSEIS    Returns valid SAClab datafile versions
%
%    Description: Returns a vector of version numbers corresponding to
%     seismic datafiles that SAClab can work with.
%
%    Notes:
%     - version 6 corresponds to SAC's v6 binary format
%     - versions 101,200,201 are modifications of SAC's v6.
%    
%    Usage:    valid_versions=vvseis()
%
%    Examples:
%
%    See also:  seischk, isseis, seishi, gv

% no args - return valid versions
valid=[6 101 200 201];

end
