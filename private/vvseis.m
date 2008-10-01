function [valid]=vvseis()
%VVSEIS    Returns valid SAClab datafile versions
%
%    Description: VVSEIS() returns a vector of version numbers 
%     corresponding to seismic datafile formats that SAClab can work with.
%
%    Notes:
%     - version 6 corresponds to SAC's v6 binary format
%     - versions 101,200,201 are modifications of SAC's v6.
%     - see SEISDEF.M for details on the formats
%    
%    System requirements: Matlab
%
%    Input/Output requirements: NONE
%
%    Header changes: N/A
%    
%    Usage:    valid_versions=vvseis()
%
%    Examples:
%     Find out the number of different versions SAClab supports:
%      length(vvseis)
%
%    See also:  seischk, isseis, seisdef, gv

%     Version History:
%        Feb. 28, 2008 - initial version
%        Mar.  4, 2008 - doc update
%        Mar. 18, 2008 - removed 2 versions to match seisdef update
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 14, 2008 at 15:25 GMT

% todo:

% no args - return valid versions
valid=[6 101 200 201];

end
