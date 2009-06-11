function [data]=prewhiten(data,order)
%PREWHITEN    Prewhiten SEIZMO data records for spectral operations
%
%    Usage:    data=prewhiten(data)
%              data=prewhiten(data,order)
%
%    Description: PREWHITEN(DATA) returns the difference between the 
%     records in DATA with and without a prediction error filter of order 6
%     applied.  The returned record is thus the unpredictable (random)
%     portion of the data (which has a significantly whiter spectrum).  The
%     predictable portion of the data is stored in the 'pef' field of the
%     data structure as the prediction error filter.  The original record
%     may be restored (as much as possible) by applying the prediction
%     error filter with an inverse filter to the whitened record (see
%     UNPREWHITEN).  The whitened record has several advantages but the
%     main one of interest in seismology is the improved stability of
%     spectral operations.  In particular, this operation will produce a
%     better conditioned matrix in a deconvolution.  For more info please
%     consider looking through the suggested reading in the Notes section.
%     Prewhitened records will have the 'prewhitened' field set to TRUE.
%
%     PREWHITEN(DATA,ORDER) allows specifying the order (number of samples
%     or poles) in the prediction error filter.  The higher the order, the
%     better the prediction error filter can capture the predictable
%     portion of the data.  The returned record will better represent a
%     random signal (with a flatter spectrum).  This does tend to converge
%     to some maximum entropy with increasing ORDER.  Restoration of the
%     original record may also degrade at higher ORDER.  Some tuning will
%     likely be required to find the best ORDER for your operation.  ORDER
%     must be a positive integer <NPTS, where NPTS is the number of points
%     in a record.  Currently ORDER must be scalar (for efficiency).
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
%     a feel for how important/detrimental it is.  The difference can be
%     found by plotting the difference:
%      plot1(subtractrecords(data,unprewhiten(prewhiten(data,order))))
%
%    See also: unprewhiten, levinson, filter, whiten

%     Version History:
%        June  8, 2009 - initial version
%        June  9, 2009 - renamed from WHITEN to PREWHITEN, doc fixes
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
%     Last Updated June  9, 2009 at 16:50 GMT

% todo:

% check number of inputs
msg=nargchk(1,2,nargin);
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

% require evenly-spaced time series, general x vs y
if(any(~strcmpi(leven,'true')))
    error('seizmo:prewhiten:illegalOperation',...
        'Illegal operation on unevenly spaced record!')
elseif(any(~strcmpi(iftype,'Time Series File')...
        & ~strcmpi(iftype,'General X vs Y file')))
    error('seizmo:prewhiten:illegalOperation',...
        'Illegal operation on spectral/xyz record!')
end

% error for already prewhitened
nrecs=numel(data);
if(isfield(data,'prewhitened') && ...
        isfield(data,'pef') && ...
        islogical([data.prewhitened]) && ...
        numel([data.prewhitened])==nrecs)
    idx=[data.prewhitened];
else
    % try to list those that are prewhitened with the correct index
    try
        % list records that are unset or false
        data(cellfun('isempty',{data.prewhitened})).prewhitened=false;
        idx=[data.prewhitened];
    catch
        % list them all
        idx=false(nrecs,1);
    end
end
if(any(idx))
    i=find(idx);
    error('seizmo:prewhiten:recordsNotWhitened',...
        ['Records: ' sprintf('%d ',i) '\n'...
        'PREWHITEN will not prewhiten prewhitened records!']);
end

% combine records for numerical efficiency
[recs,idx,ind,idx2,store,npts]=combinerecords(data);

% default order
if(nargin==1); order=6; end

% check order
if(~isnumeric(order) || any(fix(order)~=order) || ~isscalar(order))
    error('seizmo:prewhiten:badOrder','ORDER must be a scalar integer!');
end
if(order<1 || any(order>=npts))
    error('seizmo:prewhiten:badOrder','ORDER must be >0 but <NPTS!');
end

% get autocorr
nrows=2^nextpow2(2*size(recs,1)-1);
x=fft(recs,nrows);
x=ifft(abs(x).^2);
%x=ifft(abs(x).^2).*npts(idx1,ones(nrows,1)).'; % biased autocorr

% get predition error filter
a=levinson(x(1:order+1,:),order);

% implement filters separately
for i=1:size(recs,2)
    recs(:,i)=recs(:,i)-filter([0 -a(i,2:end)],1,recs(:,i));
end

% store the whitening filter at the struct level
[data.prewhitened]=deal(true);
for i=1:numel(data)
    data(i).pef=a(idx==i,:);
end

% distribute records back
data=distributerecords(data,recs,idx,[],[],store,npts);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
