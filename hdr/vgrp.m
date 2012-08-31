function [varargout]=vgrp(varargin)
%VGRP    List/Return members of SEIZMO virtual group header fields
%
%    Usage:    vgrp
%              fields=vgrp('groupfield')
%              [fields1,...,fieldN]=vgrp('groupfield1',...,'groupfieldN')
%
%    Description:
%     VGRP prints a list of the fields for each group field.  The printing
%     is done such that it is a valid command to create a cellstr array.
%
%     FIELDS=VGRP('GROUPFIELD') returns a cellstr array of the fields
%     within the group field GROUPFIELD.
%
%     [FIELDS1,...,FIELDN]=VGRP('GROUPFIELD1',...,'GROUPFIELDN') allows
%     multiple group field requests in a single call.
%
%    Notes:
%
%    Examples:
%     % I forget the order of ST & EV all the time:
%     vgrp('st')
%
%    See also: CHANGEHEADER, GETHEADER, LISTHEADER, QUERYHEADER,
%              COMPAREHEADER, SEIZMODEF

%     Version History:
%        Aug. 30, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 30, 2012 at 15:15 GMT

% todo:

% using group fields under SEIZMO's definition of SACv6
def=seizmodef('SAC Binary',6);

% print everything if no input
if(nargin==0)
    % loop over groupfields
    gf=fieldnames(def.vgrp);
    for i=1:numel(gf)
        f=def.vgrp.(gf{i});
        fprintf(['%s={' sprintf('''%s'' ',f{:}) '};\n'],gf{i});
    end
else % return group fields
    % require strings
    if(~iscellstr(varargin))
        error('seizmo:vgrp:badInput',...
            'All group field inputs must be strings!');
    end
    
    % loop over inputs
    for i=1:nargin
        try
            varargout{i}=def.vgrp.(lower(varargin{i}));
        catch
            error('seizmo:vgrp:badInput',...
                'Unknown group field: %s !',varargin{i});
        end
    end
end

end
