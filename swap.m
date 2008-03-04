function [varargout]=swap(varargin)
%SWAP    Swap values
%
%    Description: Assigns input values to output variables.  Useful for
%     trading values without intermediate variables.
%
%    Usage:  [a,d,b,c,...]=swap(b,a,d,c,...)
%
%    Examples:
%     [b,a]=swap(a,b); % trade values
%
%    See also: 

varargout=varargin;

end
