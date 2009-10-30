function [data]=prewhiten(data,order)
%PREWHITEN    Prewhiten SEIZMO data records for spectral operations
%
%    Usage:    data=prewhiten(data)
%              data=prewhiten(data,order)
%
%    Description: PREWHITEN(DATA) returns the difference between the 
%     records in DATA with and without a prediction error filter of order 6
%     applied.  The returned records are thus the unpredictable (noise)
%     portion of the data (which has a significantly whiter spectrum).  The
%     predictable portion of the data is stored as a prediction error
%     filter under the struct field .misc.pef in DATA.  The original record
%     may be restored (as much as possible) by applying the prediction
%     error filter with an inverse filter to the whitened record (see
%     UNPREWHITEN).  The whitened record has several advantages but the
%     main one of interest in seismology is the improved stability of
%     spectral operations.  In particular, this operation will produce a
%     better conditioned matrix in a deconvolution.  For more info please
%     consider looking through the suggested reading in the Notes section.
%     Prewhitened records will also have the struct field .misc.prewhitened
%     field set to TRUE.
%
%     PREWHITEN(DATA,ORDER) allows specifying the order (number of samples
%     or poles) in the prediction error filter.  The higher the order, the
%     better the prediction error filter can capture the predictable
%     portion of the data.  The returned record will better represent a
%     random signal (with a flatter spectrum).  This does tend to converge
%     to some maximum entropy with increasing ORDER.  Restoration of the
%     original record may also degrade at higher ORDER.  Some tuning will
%     likely be required to find the best ORDER for your operation.  ORDER
%     must be a positive integer or integers (one per record) <NPTS, where
%     NPTS is the number of points in a record.  Default ORDER is 6.
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
%      plot1(subtractrecords(data,unprewhiten(prewhiten(data))))
%
%    See also: UNPREWHITEN, LEVINSON, FILTER, WHITEN

%     Version History:
%        June  8, 2009 - initial version
%        June  9, 2009 - renamed from WHITEN to PREWHITEN, doc fixes
%        June 25, 2009 - update for RECORD2MAT/MAT2RECORD, process records
%                        individually rather than all together (faster)
%        Sep. 22, 2009 - pushed .pef & .prewhitened to .misc.pef &
%                        .misc.prewhitened (avoids struct cat errors)
%        Oct. 13, 2009 - minor doc update
%        Oct. 19, 2009 - global access to order
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 19, 2009 at 06:40 GMT

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
npts=getheader(data,'npts');
ncmp=getncmp(data);
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

% pull .misc field out
misc=[data.misc];

% error for already prewhitened
nrecs=numel(data);
if(isfield(misc,'prewhitened') && ...
        isfield(misc,'pef') && ...
        islogical([misc.prewhitened]) && ...
        numel([misc.prewhitened])==nrecs)
    idx=[misc.prewhitened];
else
    % try to list those that are prewhitened with the correct index
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
if(any(idx))
    i=find(idx);
    error('seizmo:prewhiten:recordsNotWhitened',...
        ['Records: ' sprintf('%d ',i) '\n'...
        'PREWHITEN will not prewhiten prewhitened records!']);
end

% default/global/check order
global SEIZMO
if(nargin==1 || isempty(order))
    order=6;
    try
        if(~isempty(SEIZMO.PREWHITEN.ORDER))
            order=SEIZMO.PREWHITEN.ORDER;
        end
    catch
    end
end
if(~isnumeric(order) || any(fix(order)~=order) ...
        || ~any(numel(order)==[1 nrecs]))
    error('seizmo:prewhiten:badOrder',...
        'ORDER must be a scalar or an array with 1 integer per record!');
end
order=order(:);
if(any(order<1 | order>=npts))
    error('seizmo:prewhiten:badOrder','ORDER must be >0 but <NPTS!');
end

% expand scalar order
if(isscalar(order))
    order(1:nrecs,1)=order;
end

% loop over records
npts2=2.^nextpow2n(2*npts);
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % get autocorr
    x=ifft(abs(fft(data(i).dep,npts2(i))).^2);
    
    % get predition error filter
    a=levinson(x(1:order(i)+1,:),order(i));
    
    % implement filter
    for j=1:ncmp(i)
        data(i).dep(:,j)=...
            data(i).dep(:,j)-filter([0 -a(j,2:end)],1,data(i).dep(:,j));
    end
    
    % store the prewhitening filter at the struct level
    data(i).misc.prewhitened=true;
    data(i).misc.pef=a;
    
    % change class back
    data(i).dep=oclass(data(i).dep);
    
    % get dep*
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
