function [valid]=validseizmo(filetype)
%VALIDSEIZMO    Returns valid SEIZMO datafile versions
%
%    Usage:    valid_versions=validseizmo(filetype)
%
%    Description: VALIDSEIZMO(FILETYPE) returns a vector of version numbers
%     of the specified filetype FILETYPE that SEIZMO can work with.  If the
%     filetype is not supported, VALIDSEIZMO will return an empty array.
%
%    Notes:
%     - currently only supports filetypes 'SAC Binary' and 'SEIZMO Binary'
%     - SEIZMO versions 101,200,201 are modifications of SAC version 6
%
%    Examples:
%     How many different SEIZMO binary file versions are supported:
%      length(validseizmo('SEIZMO Binary'))
%
%    See also:  seizmocheck, isseizmo, seizmodef, getfileversion

%     Version History:
%        Feb. 28, 2008 - initial version
%        Mar.  4, 2008 - doc update
%        Apr. 18, 2008 - drop a couple deprecated versions
%        June 12, 2008 - doc update, added history
%        Sep. 14, 2008 - minor doc update
%        Oct. 17, 2008 - filetype argument, history fix
%        Oct. 27, 2008 - minor doc update
%        Nov. 13, 2008 - renamed from VVSEIS to VALIDSEIZ
%        Nov. 15, 2008 - now VALIDSEIZMO and separated SAC from SEIZMO
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 27, 2009 - switch v101 from SEIZMO to SAC even though it is
%                        not supported by SAC -- this makes things a bit
%                        easier for multiple component support
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 27, 2009 at 20:50 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% get versions
if(strcmpi(filetype,'SAC Binary'))
    valid=[6 101];
elseif(strcmpi(filetype,'SEIZMO Binary'))
    valid=[200 201];
end

end
