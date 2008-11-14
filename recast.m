function [data]=classit(data,class)
%CLASSIT    Change in memory SAClab data storage
%
%    Description: CLASSIT(DATA,CLASS) changes the Matlab data storage of 
%     records in DATA to CLASS.  CLASS must be a string or cellstr of valid
%     Matlab class(es) ('double', 'single', etc).  This does not change the
%     storage type of data written to disk (requires a version change).
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage:    data=classit(data,class)
%
%    Examples:
%     Double precision records and fix the delta intervals
%      data=fixdelta(classit(data,'double'))
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  8, 2008 at 06:50 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin))

% check data structure
error(seischk(data,'dep'))

% number of records
nrecs=numel(data);

% check class
if(isa(class,'function_handle'))
    if(isscalar(class))
        class(1:nrecs,1)={class};
    elseif(numel(class)~=nrecs)
        error('SAClab:classit:badInput',...
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
    error('SAClab:classit:badInput',...
        'CLASS must be a string, cellstr, or function handle!');
end
if(isscalar(class))
    class(1:nrecs,1)=class;
elseif(numel(class)~=nrecs)
    error('SAClab:classit:badInput',...
        'CLASS must be scalar or match number of records in DATA!');
end

% retreive header info
leven=glgc(data,'leven');

% loop through each record
for i=1:nrecs
    % reclass dependent data
    data(i).dep=class{i}(data(i).dep);
    
    % reclass independent data
    if(strcmp(leven(i),'false'))
        data(i).ind=class{i}(data(i).ind);
    end
end

end
