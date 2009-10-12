function [data]=sortbyfield(data,field,mode)
%SORTBYFIELD   Sort SEIZMO records by a header or SEIZMO struct field
%
%    Usage:    data=sortbyfield(data,field)
%              data=sortbyfield(data,field,mode)
%
%    Description: SORTBYFIELD(DATA,FIELD) sorts SEIZMO records
%     in DATA by the header field FIELD.  Also will sort by any field in
%     data such as 'name', 'version', 'byteorder', etc.  Data fields
%     'head', 'dep', and 'ind' are not supported.  Group header fields are
%     also not supported.
%
%     SORTBYFIELD(DATA,FIELD,MODE) sets the sorting order ('ascend' or
%     'descend' is allowed - 'ascend' is the default).
%
%    Notes:
%
%    Header changes: NONE
%
%    Examples: 
%     Sort by descending degree distance:
%      data=sortbyfield(data,'gcarc','descend')
%
%    See also: SORT, GETHEADER

%     Version History:
%        Oct. 31, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Jan. 28, 2008 - code cleaning
%        Feb. 28, 2008 - better checking, doc update, any DATA field
%        Mar.  4, 2008 - minor doc update
%        Nov. 23, 2008 - updated for new name schema (now SORTBYHEADER),
%                        history fix
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:45 GMT

% todo:

% check number of args
msg=nargchk(2,3,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% set mode if none
if(nargin==2 || isempty(mode)); mode='ascend'; end

% get field values (or filenames/byte-orders/versions)
bad={'head' 'dep' 'ind'};
if(~ischar(field) || any(strcmpi(field,bad)))
    error('seizmo:sortbyfield:badField','FIELD is bad!');
elseif(isfield(data,field))
    if(isnumeric([data.(field)]))
        [values,indices]=sort([data.(field)]);
    else
        [values,indices]=sort({data.(field)});
    end
else
    [values,indices]=sort(getheader(data,field));
end

% check indices size
if(numel(indices)~=numel(data))
    error('seizmo:sortbyfield:tooManyIndices',...
        'Too many output indices!')
end

% flip if descend mode
if(strcmpi(mode,'descend'))
    indices=indices(end:-1:1);
end

% sort data
data=data(indices);

end
