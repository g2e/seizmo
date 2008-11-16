function [data]=hsort(data,field,mode)
%HSORT   Sort SEIZMO records by a header field
%
%    Description: Sorts SEIZMO data records according to a header field. 
%     Mode is optional ('ascend' or 'descend' is allowed - 'ascend' is
%     the default).  Also will sort by any field in data such as 'name',
%     'version' and 'endian'.  Data fields 'head','x' and 't' are not
%     allowed for sorting.  No header group fields.  Uses Matlab function
%     sort.
%
%    Usage:  [data]=hsort(data,field,mode)
%
%    Examples: 
%       sort as descending degree distance
%       [data]=hsort(data,'gcarc','descend')
%
%    See also: gh

% check number of args
error(nargchk(2,3,nargin))

% check data structure
error(seischk(data))

% set mode if none
if(nargin==2 || isempty(mode)); mode='ascend'; end

% get field values (or filenames/byte-orders/versions)
bad={'head' 'x' 't'};
if(~ischar(field) || any(strcmpi(field,bad)))
    error('seizmo:hsort:badField','Field not allowed')
elseif(isfield(data,field))
    if(isnumeric([data.(field)]))
        [values,indices]=sort([data.(field)]);
    else
        [values,indices]=sort({data.(field)});
    end
else
    [values,indices]=sort(gh(data,field));
end

% check indices size
if(numel(indices)~=numel(data))
    error('seizmo:hsort:tooManyIndices','Too many output indices')
end

% flip if descend mode
if(strcmpi(mode,'descend'))
    indices=indices(end:-1:1);
end

% sort data
data=data(indices);

end
