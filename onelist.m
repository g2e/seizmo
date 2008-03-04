function [varargin]=onelist(varargin)
%ONELIST    Compiles multiple char/cellstr arrays into a single column cellstr
%
%    Description:  Combines char/cellstr arrays into a single column
%     cellstr array.  Useful for taking file lists in various formats and
%     combining them into one easy to handle list.
%
%    Usage: list=onelist(list1,list2,...)
%
%    Examples:
%
%    See also: strnlen

% check, organize, and compile char/cellstr arrays
for i=1:nargin
    % check that input is char or cellstr array
    if(~ischar(varargin{i}) && ~iscellstr(varargin{i}))
        error('char and cellstr arrays only')
    end
    
    % char to cellstr then array to column vector
    varargin{i}=cellstr(varargin{i});
    varargin{i}=varargin{i}(:).';
end

% concatinate arguments
varargin=[varargin{:}].';

end
