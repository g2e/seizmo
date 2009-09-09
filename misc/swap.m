function [varargout]=swap(varargin)
%SWAP    Swap values
%
%    Usage:  [a,d,b,c,...]=swap(b,a,d,d,...)
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
%     - wrote this before I knew about DEAL...probably should toss it
%
%    Examples:
%     Do a simple swap:
%      [b,a]=swap(a,b)
%
%     Preallocate some arrays with 4x4 nan array:
%      [a,b,c,d]=swap(nan(4))
%
%    See also: deal

%     Version History:
%        Feb. 12, 2008 - initial version
%        Oct. 26, 2008 - added multiassign last input feature, doc update
%        Apr. 23, 2009 - move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 21:15 GMT

% todo:

varargout=varargin([1:nargin nargin*ones(1,nargout-nargin)]);

end
