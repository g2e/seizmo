function [varargout]=mt_v26(varargin)
%MT_V26    Converts moment tensor from Nx6 array to 6 Nx1 vectors
%
%    Usage:
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also:

%     Version History:
%        Mar.  8, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  8, 2010 at 13:50 GMT

% todo:

if(nargin==1)
    sz=size(varargin{1});
    if(~isreal(varargin{1}) || sz(2)~=6)
        error('seizmo:mt_v26:badInput',...
            'Input must be a real-valued Nx6 array!');
    end
    varargout=num2cell(varargin{1},1);
else
    error('seizmo:mt_v26:badNumInput',...
        'Incorrect number of inputs!');
end

end
