function [data]=unprewhiten(data)
%UNPREWHITEN    Undo prewhitening of SEIZMO data records
%
%    Usage:    data=unprewhiten(data)
%
%    Description: UNPREWHITEN(DATA) restores the predictable portion of
%     records in DATA by inversely applying the prediction error filter
%     stored in the .misc.pef field in DATA.  This effectively undoes
%     PREWHITEN.  Records are required to have the .misc.prewhitened' field
%     set to TRUE.  See function PREWHITEN and the suggested reading in
%     Notes for more detailed information.  The returned DATA will have the
%     .misc.pef field set to an empty array and the .misc.prewhitened field
%     set to FALSE.
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
%     Try prewhitening and unprewhitening first.  Then try comparing some
%     operation without prewhiten/unprewhiten with one including it to get
%     a feel for how important/detrimental it is.  Plotting the difference:
%      plot1(subtractrecords(data,unprewhiten(prewhiten(data,order))))
%
%    See also: PREWHITEN, LEVINSON, FILTER, WHITEN

%     Version History:
%        June  8, 2009 - initial version
%        June  9, 2009 - renamed from UNWHITEN to UNPREWHITEN, doc fixes
%        Sep. 22, 2009 - pushed .pef & .prewhitened to .misc.pef &
%                        .misc.prewhitened (avoids struct concatination
%                        errors)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 22, 2009 at 00:20 GMT

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
    error('seizmo:unprewhiten:illegalOperation',...
        'Illegal operation on unevenly spaced record!')
elseif(any(~strcmpi(iftype,'Time Series File')...
        & ~strcmpi(iftype,'General X vs Y file')))
    error('seizmo:unprewhiten:illegalOperation',...
        'Illegal operation on spectral/xyz record!')
end

% pull .misc field out
misc=[data.misc];

% error for nonprewhitened
nrecs=numel(data);
if(isfield(misc,'prewhitened') && ...
        isfield(misc,'pef') && ...
        islogical([misc.prewhitened]) && ...
        numel([misc.prewhitened])==nrecs)
    idx=[misc.prewhitened];
else
    % prepare a decent list for the error msg
    try
        % list records that are unset or false
        [misc(cellfun('isempty',...
            {misc.prewhitened})).prewhitened]=deal(false);
        idx=[misc.prewhitened];
    catch
        % list them all
        idx=false(nrecs,1);
    end
end
if(any(~idx))
    i=find(~idx);
    error('seizmo:unprewhiten:recordsNotPrewhitened',...
        ['Records: ' sprintf('%d ',i) '\n'...
        'Cannot unprewhiten non-prewhitened records!']);
end

% loop through whitened records
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=find(idx)
    % get class
    oclass=str2func(class(data(i).dep));
    
    % check pef matches ncmp
    if(size(data(i).misc.pef,1)~=ncmp(i))
        error('seizmo:unprewhiten:ncmpInconsistent',...
            'Record: %d\nNCMP has changed since WHITEN operation!',i);
    end
    
    % unwhiten filter
    for j=1:ncmp(i)
        data(i).dep(:,j)=oclass(...
            filter(1,data(i).misc.pef(j,:),double(data(i).dep(:,j))));
    end
    
    % unset prewhitened, clear pef
    data(i).misc.prewhitened=false;
    data(i).misc.pef=[];
    
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
