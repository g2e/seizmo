function [data]=addcorrection(data,file,field)
%ADDCORRECTION    Adds travel time corrections to SEIZMO headers

%     Version History:
%        Dec. 10, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec. 10, 2009 at 02:45 GMT

% todo:

% check nargin
msg=nargchk(3,3,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% get record filenames
names={data.name}';

% read in file
txt=readtxt(file);

% separate filename and corrections
words=reshape(getwords(txt),2,[])';

% now we need to match a record to a correction
% - returns last matching entry
[tf,idx]=ismember(names,words(:,1));

% error for records without an entry
if(any(~tf))
    error('seizmo:addcorrection:missingCorrection',...
        ['Record(s):\n' sprintf('%d ',find(~tf)) ...
        '\nNo matching entry found!']);
end

% check field
if(~ischar(field) || size(field,1)~=1 || ndims(field)>2)
    error('seizmo:addcorrection:badField',...
        'FIELD must be a string!');
end

% add to header
data=changeheader(data,field,str2double(words(idx,2)));

end
