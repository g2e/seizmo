function [valid]=validseizmo(filetype)
%VALIDSEIZMO    Returns valid SEIZMO datafile filetypes or versions
%
%    Usage:    valid_filetypes=validseizmo()
%              valid_versions=validseizmo(filetype)
%
%    Description:
%     VALID_FILETYPES=VALIDSEIZMO() returns a cellstr array of filetypes
%     that SEIZMO will work with.
%
%     VALID_VERSIONS=VALIDSEIZMO(FILETYPE) returns a vector of version
%     numbers of the specified filetype FILETYPE that SEIZMO can work with.
%     If the filetype is not supported, VALIDSEIZMO will return an empty
%     array.
%
%    Notes:
%     - currently only supports filetypes 'SAC Binary' and 'SEIZMO Binary'
%     - SEIZMO versions 101,200,201 are modifications of SAC version 6
%
%    Examples:
%     % What filetypes are supported?
%     validseizmo()
%
%     % How many different SEIZMO binary file versions are supported:
%     length(validseizmo('SEIZMO Binary'))
%
%    See also:  SEIZMOCHECK, ISSEIZMO, SEIZMODEF, GETFILEVERSION

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
%        Sep. 25, 2009 - undid hack mentioned above
%        Oct.  6, 2009 - new usage form to return valid filetypes
%        Oct. 16, 2009 - dropped redundant nargin check (2x faster)
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 02:25 GMT

% todo:

% check nargin
if(nargin)
    % handle empty
    if(isempty(filetype)); valid=[]; return; end
    
    % get versions
    switch lower(filetype)
        case 'sac binary'
            valid=6;
        case 'seizmo binary'
            valid=[101 200 201];
        otherwise
            valid=[];
    end
else
    valid={'SAC Binary'; 'SEIZMO Binary'};
end

end
