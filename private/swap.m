function [varargout]=swap(varargin)
%SWAP    Swap values
%
%    Description: SWAP(INPUT1,...,INPUTN) assigns input variables' values
%     directly to the output variables.  If there are more input variables
%     than output variables, the excess input variables are ignored.  If
%     there are more output veriables than input variables, the extra
%     output variables share the last input variable.  This is very useful
%     for intricate swapping of values without bothering with intermediate
%     variables.  It is also great for mass preallocating.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:  [a,d,b,c,...]=swap(b,a,d,d,...)
%
%    Examples:
%     Do a simple swap:
%      [b,a]=swap(a,b)
%
%     Preallocate some arrays with 4x4 nan array:
%      [a,b,c,d]=swap(nan(4))
%
%    See also: 

%     Version History:
%        Feb. 12, 2008 - initial version
%        Oct. 26, 2008 - added multiassign last input feature, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 26, 2008 at 08:05 GMT

% todo:

varargout=varargin([1:nargin nargin*ones(1,nargout-nargin)]);

end
