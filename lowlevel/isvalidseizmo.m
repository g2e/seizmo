function [lgc]=isvalidseizmo(filetype,version)
%ISVALIDSEIZMO    TRUE for valid filetype/version combinations
%
%    Usage:    lgc=isvalidseizmo(filetype,version)
%
%    Description:
%     LGC=ISVALIDSEIZMO(FILETYPE,VERSION) returns a logical array
%     indicating which elements in FILETYPE and VERSION are valid
%     filetype/version combinations in SEIZMO.  FILETYPE and VERSION should
%     have the same number of elements or be scalar (scalar expansion is
%     enabled).  LGC is a column vector with the same number of elements as
%     FILETYPE and VERSION, where LGC(i) indicates if FILETYPE(i) + 
%     VERSION(i) is valid.
%
%    Notes:
%
%    Examples:
%     % Are any SEIZMO versions valid SAC versions?:
%     isvalidseizmo('SAC Binary',validseizmo('SEIZMO Binary'))
%
%     % Quickly check that all filetype and
%     % version fields in a dataset are ok:
%     all(isvalidseizmo({data.filetype},{data.version});
%
%    See also: VALIDSEIZMO, SEIZMOCHECK, ISSEIZMO, SEIZMODEF,
%              GETFILEVERSION

%     Version History:
%        Oct.  6, 2009 - initial version
%        Aug. 21, 2010 - nargchk fix
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 18:10 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check filetype/version
if(ischar(filetype)); filetype=cellstr(filetype); end
if(~iscellstr(filetype))
    error('seizmo:isvalidseizmo:badInput',...
        'FILETYPE must be a char/cellstring array!');
end
if(iscell(version)); version=cell2mat(version); end
if(~isreal(version))
    error('seizmo:isvalidseizmo:badInput',...
        'VERSION must be a real array!');
end

% get sizes and expand
filetype=filetype(:);
version=version(:);
nft=numel(filetype);
nv=numel(version);
if(nft==1); filetype=filetype(ones(nv,1),1); end
if(nv==1); version=version(ones(nft,1),1); end
if(~isequal(size(filetype),size(version)))
    error('seizmo:isvalidseizmo:badInput',...
        ['FILETYPE & VERSION must be scalar or ' ...
        'have an equal number of elements!']);
end
n=max(nft,nv);

% loop over filetype/versions
lgc=false(n,1);
for i=1:n; lgc(i)=any(version(i)==validseizmo(filetype{i})); end

end
