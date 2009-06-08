function [data]=unwhiten(data)
%UNWHITEN    Undo spectral whitening of SEIZMO data records
%
%    Description: UNWHITEN(DATA) uses an inverse filter to unwhiten the
%     records in DATA.  That inverse filter is found in the 'pef' field of
%     DATA.  Records are required to have the 'whitened' field set to true.
%     See function WHITEN and the suggested reading in Notes for more info.
%
%    Notes:
%     - Suggested Reading:
%       - Vaidyanathan, P. P., The Theory of Linear Prediction, Synthesis
%         Lectures on Signal Processing #3, Morgan and Claypool Publishers,
%         183 pp.
%
%    Header changes: DEPMIN, DEPMAX, DEPMEN
%
%    Examples:
%     Try whitening and unwhitening first.  Then try comparing some
%     operation without whiten/unwhiten with one including it to get a
%     feel for how important/detrimental it is.  
%      plot1(unwhiten(whiten(data(1),order)))
%      plot1(subtractrecords(data(1),unwhiten(whiten(data(1),order))))
%
%    See also: whiten, levinson, filter

%     Version History:
%        June  8, 2009 - initial version
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  8, 2009 at 21:30 GMT

% todo:

% check number of inputs
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% get some header fields
leven=getlgc(data,'leven');
iftype=getenumdesc(data,'iftype');
ncmp=getncmp(data);

% require evenly-spaced time series, general x vs y
if(any(~strcmpi(leven,'true')))
    error('seizmo:unwhiten:illegalOperation',...
        'Illegal operation on unevenly spaced record!')
elseif(any(~strcmpi(iftype,'Time Series File')...
        & ~strcmpi(iftype,'General X vs Y file')))
    error('seizmo:unwhiten:illegalOperation',...
        'Illegal operation on spectral/xyz record!')
end

% error for nonwhitened
nrecs=numel(data);
if(isfield(data,'whitened') && ...
        isfield(data,'pef') && ...
        islogical([data.whitened]) && ...
        numel([data.whitened])==nrecs)
    idx=[data.whitened];
else
    % prepare a decent list for the error msg
    try
        % list records that are unset or false
        data(cellfun('isempty',{data.whitened})).whitened=false;
        idx=[data.whitened];
    catch
        % list them all
        idx=false(nrecs,1);
    end
end
if(any(~idx))
    i=find(~idx);
    error('seizmo:unwhiten:recordsNotWhitened',...
        ['Records: ' sprintf('%d ',i) '\n'...
        'Cannot unwhiten non-whitened records!']);
end

% loop through whitened records
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=find(idx)
    % get class
    oclass=str2func(class(data(i).dep));
    
    % check pef matches ncmp
    if(size(data(i).pef,1)~=ncmp(i))
        error('seizmo:unwhiten:ncmpInconsistent',...
            'Record: %d\nNCMP has changed since WHITEN operation!',i);
    end
    
    % unwhiten filter
    for j=1:ncmp(i)
        data(i).dep(:,j)=oclass(...
            filter(1,data(i).pef(j,:),double(data(i).dep(:,j))));
    end
    
    % reset whitened, clear pef
    data(i).whitened=false;
    data(i).pef=[];
    
    % update dep*
    if(isempty(data(i).dep)); continue; end
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=changeheader(data,'depmen',depmen,'depmin',depmin,'depmax',depmax);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
