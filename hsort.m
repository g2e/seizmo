function [data]=hsort(data,field,mode)
%HSORT   Sort SAClab records by a header field
%
%    Description: Sorts SAClab data records according to a header field. 
%     Mode is optional ('ascend' or 'descend' is allowed - 'ascend' is
%     the default).  Also will sort by filename if field is set to 'name'.
%     Uses the Matlab function sort.
%
%    Usage:  [data]=hsort(data,field,mode);
%
%    Example: 
%       sort for descending degree distance
%       [data]=hsort(data,'gcarc','descend');
%
%    by Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: gh

% check number of args
error(nargchk(1,3,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head'))
    error('data structure does not have proper fields')
end

% do nothing on no sorting field
if(nargin==1); return; end

% set mode if none
if(nargin==2 || isempty(mode)); mode='ascend'; end

% get field values (or filenames)
if(strcmpi(field,'name') && isfield(data,'name'))
    [values,indices]=sort({data.name});
else
    [values,indices]=sort(gh(data,field));
end

% flip if descend mode
if(strcmpi(mode,'descend'))
    indices=indices(end:-1:1);
end

% sort data
data=data(indices);

end
