function [s]=getsubfield(s,varargin)
%GETSUBFIELD    Get substructure field contents
%
%    Usage:    f=getsubfield(s,'field','subfield',...)
%
%    Description:
%     F=GETSUBFIELD(S,'FIELD','SUBFIELD',...) extracts the subfield
%     S.('FIELD').('SUBFIELD').  Substructures should be concatenateable -
%     meaning they must have the same fields and dimensions 1 and 3+ must
%     be equal sized.  Deeper subfields may be extracted by providing more
%     fields.  F by default is returned as a cell array, unless all cells
%     contain non-empty elements of class struct, logical, or double.
%
%    Notes:
%
%    Examples:
%     % Get subfield c, which contains both char and double types:
%     a(1).b(2).c='fff';
%     a(2).b(1).c=1;
%     a(3).b(5).c='b';
%     f=getsubfield(a,'b','c')
%
%    See also: GETFIELD

%     Version History:
%        Sep. 22, 2009 - initial version
%        Sep. 25, 2009 - unwraps logical array now
%        Feb. 11, 2011 - mass nargchk fix
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check inputs
if(~isstruct(s))
    error('seizmo:getsubfield:badInput','S must be a struct!');
elseif(~iscellstr(varargin))
    error('seizmo:getsubfield:badInput','FIELD must be a string!');
end

% loop over varargin, extracting subfields
for i=1:numel(varargin)
    % check if field
    if(~isfield(s,varargin{i}))
        error('seizmo:getsubfield:badField',...
            'Subfield %s not found!',varargin{i});
    end
    
    % extract and concatenate
    if(i<nargin-1)
        % extract substructure
        s=[s.(varargin{i})];
    else
        % extract final depth
        s={s.(varargin{end})};
        
        % alter output based on final values
        if(iscellstr(s) || any(cellfun('isempty',s)))
            return;
        elseif(all(cellfun('isclass',s,'struct')))
            s=cell2mat(s);
        elseif(all(cellfun('isclass',s,'double')))
            s=cell2mat(s);
        elseif(all(cellfun('isclass',s,'logical')))
            s=cell2mat(s);
        else
            return;
        end
    end
end

end
