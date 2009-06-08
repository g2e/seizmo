function [data]=whiten(data,order)
%WHITEN    Spectral whitening of time-series SEIZMO data records
%
%    Description: WHITEN(DATA) returns the difference between the records
%     in DATA with and without a prediction error filter of order 6
%     applied.  This signal has several advantages but the main one of
%     interest in seismology is spectral whitening.  Spectral whitening is
%     useful for deconvolutions and other inversions because it will
%     produce a better conditioned matrix for those inversions which
%     stabilizes the operation. An inverse filter to unwhiten can then be
%     applied after the operations using UNWHITEN.  For more info please
%     consider looking through the suggested book(s) in the Notes section.
%
%     WHITEN(DATA,ORDER) allows specifying the order (number of poles) in
%     the prediction error filter.  The higher the order, the flatter the
%     resulting record's spectra will be but possibly at the cost of losing
%     important information in the original record.  Some tuning will
%     likely be required.  ORDER must be >0 and <NPTS, where NPTS is the
%     number of points in a record.  ORDER above 12 is allowed!  Currently
%     ORDER must be a scalar (for operational efficiency).
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
%    See also: unwhiten, levinson, filter

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
    error('seizmo:whiten:illegalOperation',...
        'Illegal operation on unevenly spaced record!')
elseif(any(~strcmpi(iftype,'Time Series File')...
        & ~strcmpi(iftype,'General X vs Y file')))
    error('seizmo:whiten:illegalOperation',...
        'Illegal operation on spectral/xyz record!')
end

% error for already whitened
nrecs=numel(data);
if(isfield(data,'whitened') && ...
        isfield(data,'pef') && ...
        islogical([data.whitened]) && ...
        numel([data.whitened])==nrecs)
    idx=[data.whitened];
else
    % try to list those that are whitened with the correct index
    try
        % list records that are unset or false
        data(cellfun('isempty',{data.whitened})).whitened=false;
        idx=[data.whitened];
    catch
        % list them all
        idx=false(nrecs,1);
    end
end
if(any(idx))
    i=find(idx);
    error('seizmo:whiten:recordsNotWhitened',...
        ['Records: ' sprintf('%d ',i) '\n'...
        'WHITEN will not whiten whitened records!']);
end

% combine records for numerical efficiency
[recs,idx,ind,idx2,store,npts]=combinerecords(data);

% default order
if(nargin==1); order=6; end

% check order
if(~isnumeric(order) || any(fix(order)~=order) || ~isscalar(order))
    error('seizmo:whiten:badOrder','ORDER must be a scalar integer!');
end
if(order<1 || any(order>=npts))
    error('seizmo:whiten:badOrder','ORDER must be >0 but <NPTS!');
end

% get biased autocorr
nrows=2^nextpow2(2*size(recs,1)-1);
x=fft(recs,nrows);
x=ifft(abs(x).^2);
%x=ifft(abs(x).^2).*npts(idx1,ones(nrows,1)).';

% get predition error filter
a=levinson(x(1:order+1,:),order);

% implement filters separately
for i=1:size(recs,2)
    recs(:,i)=recs(:,i)-filter([0 -a(i,2:end)],1,recs(:,i));
end

% store the whitening filter at the struct level
[data.whitened]=deal(true);
for i=1:numel(data)
    data(i).pef=a(idx==i,:);
end

% distribute records back
data=distributerecords(data,recs,idx,[],[],store,npts);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
