function [varargout]=expandscalars(varargin)
%EXPANDSCALARS    Expands scalars to match size of array inputs
%
%    Usage:    [a,b,c,...]=expandscalars(a,b,c,...)
%
%    Description:
%     [A,B,C,...]=EXPANDSCALARS(A,B,C,...) will expand any scalar inputs to
%     match the size of nonscalar inputs.  Note that all nonscalar inputs
%     must be the same size (checked via ISEQUALSIZEORSCALAR).
%
%    Notes:
%
%    Examples:
%     % Expand two scalars to 6x6 matrices:
%     a=magic(6);
%     b=999;
%     c=zeros(6);
%     [a,b,c,d]=expandscalars(a,b,c,3)
%
%    See also: ISEQUALSIZEORSCALAR, ISSCALAR, SIZE

%     Version History:
%        May  16, 2010 - initial version
%        Feb.  5, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  5, 2012 at 11:25 GMT

% todo:

% check inputs
[tf,sz,ns]=isequalsizeorscalar(varargin{:});
if(~tf)
    error('seizmo:expandscalars:badInput',...
        'Inputs must be equal sized or scalar!');
end

% now expand the scalars
for i=find(~ns); varargin{i}=varargin{i}(ones(sz)); end

% output
varargout=varargin;

end
