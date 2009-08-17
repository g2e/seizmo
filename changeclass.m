function [data]=changeclass(data,class)
%CHANGECLASS    Change SEIZMO data storage (in memory)
%
%    Usage:    data=changeclass(data,class)
%
%    Description: CHANGECLASS(DATA,CLASS) changes the class type of SEIZMO
%     records in DATA to CLASS.  CLASS must be a string or cellstr of valid
%     class(es) ('double', 'single', etc).  This does not change the
%     storage type of data written to disk (requires a version change).
%
%    Notes:
%
%    Header changes: NONE
%
%    Examples:
%     Double precision records and fix the delta intervals
%      data=fixdelta(changeclass(data,'double'))
%
%    See also: fixdelta

%     Version History:
%        Feb. 21, 2008 - initial version
%        Feb. 23, 2008 - uses GLGC now
%        Feb. 28, 2008 - uses SEISCHK now
%        Mar.  4, 2008 - doc update, fixed LEVEN bug, uses LGCCHK now
%        June  6, 2008 - minor code cleaning
%        Oct.  8, 2008 - changed name from DOUBLEIT to CLASSIT, doc update,
%                        dropped reclassing header (keep it double!), allow
%                        data to be changed to any Matlab supported class,
%                        drop LGCCHK
%        Nov. 22, 2008 - renamed from CLASSIT to CHANGECLASS
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:45 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% number of records
nrecs=numel(data);

% check class and convert to function handle
if(isa(class,'function_handle'))
    if(isscalar(class))
        class(1:nrecs,1)={class};
    elseif(numel(class)~=nrecs)
        error('seizmo:changeclass:badInput',...
            'CLASS must be scalar or match number of records in DATA!');
    else
        class=mat2cell(class(:),ones(1,nrecs),1);
    end
elseif(ischar(class) || iscellstr(class))
    class=cellstr(class);
    for i=1:numel(class)
        class(i)={str2func(class{i})};
    end
else
    error('seizmo:changeclass:badInput',...
        'CLASS must be a string, cellstr, or function handle!');
end
if(isscalar(class))
    class(1:nrecs,1)=class;
elseif(numel(class)~=nrecs)
    error('seizmo:changeclass:badInput',...
        'CLASS must be scalar or match number of records in DATA!');
end

% reclass dependent data
for i=1:nrecs
    data(i).dep=class{i}(data(i).dep);
end

% reclass independent data
if(isfield(data,'ind'))
    for i=1:nrecs
        data(i).ind=class{i}(data(i).ind);
    end
end

end
