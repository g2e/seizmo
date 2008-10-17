function [data]=slidingmean(data,nsamples,varargin)
%SLIDINGMEAN    Returns sliding-window mean of SAClab data records
%
%    Description: SLIDINGMEAN(DATA,N) applies a centered sliding-window
%     mean of 2N+1 samples to the dependent component(s) of SAClab data
%     records in DATA.  N can be a scalar (each record has the same window
%     size) or a vector (define each record's window size separately).
%     Sliding windows extending outside the record are truncated (look at
%     'EDGE' option to change this).
%
%     SLIDINGMEAN(...,'POSITION','CENTER'|'TRAIL'|'LEAD') sets the position
%     of the sliding window relative to the reference data point.  CENTER 
%     positions the window such that the reference point is at its center.
%     TRAIL positions the window to trail the reference point such that the
%     reference point has the highest index in the window.  LEAD sets the
%     window to lead the reference point such that the reference point has 
%     the lowest index in the window.  Note that the window size for a
%     CENTER positioning is 2N+1, while for TRAIL or LEAD the window size
%     is N.  Default position is CENTER.
%     
%     SLIDINGMEAN(...,'OFFSET',OFFSET) sets the offset of the sliding
%     window from to the reference data point in number of samples.  For a
%     centered window (see option 'POSITION') this introduces a gap of
%     2*OFFSET-1 in the window.  For example an OFFSET of 1 will exclude
%     the reference data point from the sliding window.  Negative OFFSETS
%     are allowed for a centered window, but they are complicated due to
%     overlap.  For example an OFFSET of -1 will make the window include
%     the reference data point twice and an OFFSET of -2 will cause the 3
%     centermost points to be included twice and so on.  Default OFFSET=0.
%     OFFSET may be a vector of offsets specifying each record's offset.
%     
%     SLIDINGMEAN(...,'EDGE','TRUNCATE'|'PAD') sets how to handle the edge
%     cases.  TRUNCATE eliminates points in the sliding-window that do not
%     reference a datapoint (for instance if the window extends before or
%     after the data, that portion of the window will be truncated).  PAD
%     adds zeros to the data so that all the points of the sliding-window
%     always reference some value.  Default setting is TRUNCATE.
%
%    Notes:
%     - Centered windows are of length 2N+1, while the others are just N
%     - SLIDINGMEAN is faster than SLIDEFUN because it uses Matlab's FILTER
%
%    System requirements: Matlab 7
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Usage:    data=slidingmean(data,n)
%              data=slidingmean(...,'position','center'|'trail'|'lead')
%              data=slidingmean(...,'offset',offset)
%              data=slidingmean(...,'edge','truncate'|'pad')
%
%    Examples:
%      Compare an envelope and a 21-sample sliding-window mean:
%       p2([envelope(data(1)) slidingmean(data(1),10)])
%
%    See also: slidingrms, slidingam, slidefun

%     Version History:
%        Oct.  7, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  7, 2008 at 00:15 GMT

% todo:

% check nargin
error(nargchk(2,8,nargin))

% check data structure
error(seischk(data,'dep'))

% number of records
nrecs=numel(data);

% fix/check nsamples
if(isscalar(nsamples)); nsamples=nsamples(ones(1,nrecs),1); end
if(~isnumeric(nsamples) || any(nsamples<1) || numel(nsamples)~=nrecs ...
        || any(fix(nsamples)~=nsamples))
    error('SAClab:slidingmean:badInput',...
        'N must be positive integer(s)')
end

% option defaults
option.POSITION='center';
option.OFFSET=0;
option.EDGE='truncate';

% get options set by SACLAB global
global SACLAB; fields=fieldnames(option).';
if(isfield(SACLAB,'SLIDINGMEAN'))
    for i=fields
        if(isfield(SACLAB.SLIDINGMEAN,i{:})); 
            option.(i{:})=SACLAB.SLIDINGMEAN.(i{:});
        end
    end
end

% options must be field-value pairs
nargopt=length(varargin);
if(mod(nargopt,2))
    error('SAClab:slidingmean:badNumOptions','Unpaired option!')
end

% get inline options
for i=1:2:nargopt
    varargin{i}=upper(varargin{i});
    if(isfield(option,varargin{i}))
        option.(varargin{i})=varargin{i+1};
    else
        warning('SAClab:slidingmean:badInput',...
            'Unknown Option: %s !',varargin{i}); 
    end
end

% check options
if(~ischar(option.POSITION) || ...
        ~any(strcmpi(option.POSITION,{'center' 'trail' 'lead'})))
    error('SAClab:slidingmean:badOptionValue',...
        'POSITION option must be ''CENTER'', ''TRAIL'' or ''LEAD''!');
end
if(isscalar(option.OFFSET))
    option.OFFSET=option.OFFSET(ones(1,nrecs),1); 
end
if(~isnumeric(option.OFFSET) || numel(option.OFFSET)~=nrecs || ...
        any(option.OFFSET~=fix(option.OFFSET)))
    error('SAClab:slidingmean:badOptionValue',...
        'OFFSET option must be integer(s)!');
end
if(~ischar(option.EDGE) || ...
        ~any(strcmpi(option.EDGE,{'pad' 'truncate'})))
    error('SAClab:slidingmean:badOptionValue',...
        'EDGE option must be ''PAD'' or ''TRUNCATE''');
end

% loop through each record
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % get storage class of data
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % build window indices
    switch option.POSITION
        case 'center'
            hw=(1:nsamples(i))+option.OFFSET(i);
            window=[-fliplr(hw) zeros(~logical(option.OFFSET(i)),1) hw];
        case 'trail'
            window=(1-nsamples(i):0)+option.OFFSET(i);
        case 'lead'
            window=(0:nsamples(i)-1)+option.OFFSET(i);
    end
    
    % translate window to filter input
    nw=numel(window);
    minw=min(window);
    maxw=max(window);
    aminw=abs(minw);
    rangew=maxw-minw+1;
    lw=zeros(1,rangew);
    lw(window+aminw+1)=1;
    
    % pad/cut
    [sz1,sz2]=size(data(i).dep);
    minpt=1+minw;
    maxpt=sz1+maxw;
    gminpt=max([minpt 1]);
    gmaxpt=min([maxpt sz1]);
    temp=zeros(maxpt-minpt+1,sz2);
    tminpt=gminpt-minpt+1;
    trange=max([gmaxpt-gminpt+1 0]);
    temp(tminpt:tminpt-1+trange,:)=data(i).dep(gminpt:gmaxpt,:);
    
    % create divisors to average the sum
    divisors=ones(sz1,sz2)*nw;
    switch lower(option.EDGE)
        case 'truncate'
            % find where window is out of data range
            for j=1:nw
                gminw=min([max([-window(j) 0]) sz1]);
                gmaxw=max([sz1-max([window(j) 0])+1 1]);
                divisors(1:gminw,:)=divisors(1:gminw,:)-1;
                divisors(gmaxw:sz1,:)=divisors(gmaxw:sz1,:)-1;
            end
    end
    
    % filter to avoid loop (does a running selective sum)
    temp=filter(lw,1,temp);
    
    % trim off excess, mean
    temp=temp(end-sz1+1:end,:)./divisors;

    % assign back to data and change storage back
    data(i).dep=oclass(temp);
    
    % dep*
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
data=ch(data,'depmen',depmen,'depmin',depmin,'depmax',depmax);

end
