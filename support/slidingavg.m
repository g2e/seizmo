function [y]=slidingavg(x,nsamples,varargin)
%SLIDINGAVG    Returns sliding-window average of data
%
%    Usage:    y=slidingavg(x,n)
%              y=slidingavg(...,'position','center'|'trail'|'lead')
%              y=slidingavg(...,'offset',offset)
%              y=slidingavg(...,'edge','truncate'|'pad')
%              y=slidingavg(...,'dim',n)
%              y=slidingavg(...,'custom',window)
%
%    Description: SLIDINGAVG(X,N) applies a centered sliding-window average
%     of 2N+1 samples down the columns of numeric array X.  Sliding windows
%     extending outside the record are truncated (look at the 'EDGE' option
%     to change this).
%
%     SLIDINGAVG(...,'POSITION','CENTER'|'TRAIL'|'LEAD') sets the position
%     of the sliding window relative to the reference data point (the data
%     point that the average of the current window is assigned to).  CENTER 
%     positions the window such that the reference point is at its center.
%     TRAIL positions the window to trail the reference point such that the
%     reference point has the highest index in the window.  LEAD sets the
%     window to lead the reference point such that the reference point has 
%     the lowest index in the window.  Note that the window size for a
%     CENTER positioning is 2N+1 samples, while for TRAIL or LEAD the
%     window size is N samples.  Default position is CENTER.
%     
%     SLIDINGAVG(...,'OFFSET',OFFSET) sets the offset of the sliding
%     window from the reference data point in number of samples.  For a
%     centered window (see option 'POSITION') this introduces a gap of
%     2*OFFSET-1 samples in the center of the window.  For example an
%     OFFSET of 1 will exclude the reference data point from the sliding
%     window.  An OFFSET of 0 is treated specially and will not alter the
%     sliding window.  Negative OFFSETS are allowed for a centered window,
%     but they are complicated due to overlap.  For example an OFFSET of -1
%     will make the window include the reference data point twice and an
%     OFFSET of -2 will cause the 3 centermost points to be included twice
%     and so on.  Default OFFSET=0.
%     
%     SLIDINGAVG(...,'EDGE','TRUNCATE'|'PAD') sets how to handle the edge
%     cases.  TRUNCATE eliminates points in the sliding window that do not
%     reference a datapoint (for instance if the window extends before or
%     after the data, that portion of the window will be truncated).  PAD
%     adds zeros to the data so that all the points of the sliding-window
%     always reference some value.  Default setting is TRUNCATE.
%
%     SLIDINGAVG(...,'DIM',N) specifies an alternative dimension to slide
%     across.  Default seeting is N=1 (slides across rows).
%
%     SLIDINGAVG(...,'CUSTOM',WINDOW) allows a custom sliding window
%     average.  This might be useful for a Gaussian average or similar.
%     WINDOW must be formatted as [index; weight] where index is relative
%     to the reference data point and weight does not include the averaging
%     divisor (this will be automatically computed).  An example WINDOW:
%        [ -3  -2  -1    0   1  2    3;
%         0.1   1  10  100  10  1  0.1]
%     gives a severe weighting to the center point (reference data point).
%
%    Notes:
%     - Centered windows are of length 2N+1, while the others are just N
%     - SLIDINGAVG is much faster than SLIDINGFUN but less flexible
%
%    Tested on: Matlab r2007b
%
%    Examples:
%     Get a smoothed amplitude spectra:
%      y=slidingavg(abs(fft(x)))
%
%     An example of a custom call:
%      y=slidingavg(x,[],'custom',[-10:10; gausswin(21).'])
%
%    See also: slidingrms, slidingabsmean, slidingfun

%     Version History:
%        Oct.  7, 2008 - initial version
%        Nov. 13, 2008 - renamed from SLIDINGMEAN to SLIDINGAVG, made it a
%                        normal function, added options DIM and CUSTOM
%        Nov. 15, 2008 - update for new toolbox name
%        Dec.  3, 2008 - fixed offset bugs for center, fix offset bug for
%                        trail/lead, fixed bug for handling custom filters,
%                        fix bug in handling multiple entries for the same
%                        element, doc update
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 22:20 GMT

% todo:

% check nargin
msg=nargchk(2,12,nargin);
if(~isempty(msg)); error(msg); end

% quick return if empty
if(isempty(x)); y=x; return; end

% check input array
if(~isnumeric(x))
    error('seizmo:slidingavg:badInput','X must be numeric!');
end

% option defaults
option.POSITION='center';
option.OFFSET=0;
option.EDGE='truncate';
option.DIM=1;
option.CUSTOM=[];

% get options set by SEIZMO global
global SEIZMO; fields=fieldnames(option).';
if(isfield(SEIZMO,'SLIDINGAVG'))
    for i=fields
        if(isfield(SEIZMO.SLIDINGAVG,i{:})); 
            option.(i{:})=SEIZMO.SLIDINGAVG.(i{:});
        end
    end
end

% options must be field-value pairs
nargopt=length(varargin);
if(mod(nargopt,2))
    error('seizmo:slidingavg:badNumOptions','Unpaired option!')
end

% get inline options
for i=1:2:nargopt
    varargin{i}=upper(varargin{i});
    if(isfield(option,varargin{i}))
        option.(varargin{i})=varargin{i+1};
    else
        warning('seizmo:slidingavg:badInput',...
            'Unknown Option: %s !',varargin{i}); 
    end
end

% check options
if(~ischar(option.POSITION) || ...
        ~any(strcmpi(option.POSITION,{'center' 'trail' 'lead'})))
    error('seizmo:slidingavg:badOptionValue',...
        'POSITION option must be ''CENTER'', ''TRAIL'' or ''LEAD''!');
end
if(~isnumeric(option.OFFSET) || ~isscalar(option.OFFSET) || ...
        option.OFFSET~=fix(option.OFFSET))
    error('seizmo:slidingavg:badOptionValue',...
        'OFFSET option must be an integer!');
end
if(~ischar(option.EDGE) || ...
        ~any(strcmpi(option.EDGE,{'pad' 'truncate'})))
    error('seizmo:slidingavg:badOptionValue',...
        'EDGE option must be ''PAD'' or ''TRUNCATE''!');
end
if(~isnumeric(option.DIM) || ~isscalar(option.DIM) || option.DIM<1 || ...
        option.DIM~=fix(option.DIM))
    error('seizmo:slidingavg:badOptionValue',...
        'DIM option must be a positive integer!');
end
csz=size(option.CUSTOM);
if(~isempty(option.CUSTOM) && (~isnumeric(option.CUSTOM) || ...
        numel(csz)>2 || csz(1)~=2 || ...
        any(option.CUSTOM(1,:)~=fix(option.CUSTOM(1,:)))))
    error('seizmo:slidingavg:badOptionValue',...
        'CUSTOM option must be a 2xN numeric array!');
end

% build window if not custom
if(isempty(option.CUSTOM))
    % fix/check nsamples
    if(~isnumeric(nsamples) || ~isscalar(nsamples)...
            || nsamples<1 || fix(nsamples)~=nsamples)
        error('seizmo:slidingavg:badInput',...
            'N must be a positive integer!')
    end
    
    % build window indices
    switch option.POSITION
        case 'center'
            hw=(1:nsamples)+option.OFFSET-(option.OFFSET>0);
            window=[-fliplr(hw) zeros(~option.OFFSET,1) hw];
        case 'trail'
            window=(1-nsamples:0)+option.OFFSET;
        case 'lead'
            window=(0:nsamples-1)+option.OFFSET;
    end
    sw=numel(window); nw=sw;
    weights=ones(sw,1);
else
    window=option.CUSTOM(1,:);
    weights=option.CUSTOM(2,:);
    sw=sum(weights); nw=numel(window);
end

% translate window to filter input
minw=min(window);
maxw=max(window);
rangew=maxw-minw+1;
lw=zeros(1,rangew);
owidx=maxw-window+1;
for i=1:nw
    lw(owidx(i))=lw(owidx(i))+weights(i);
end

% permute x
sz=size(x);
perm=1:numel(sz);
perm(option.DIM)=[];
perm=[option.DIM perm];
x=permute(x,perm);
sz1=sz(option.DIM);
sz(option.DIM)=[];

% pad/cut
minpt=1+minw;
maxpt=sz1+maxw;
gminpt=max([minpt 1]);
gmaxpt=min([maxpt sz1]);
y=zeros([maxpt-minpt+1 sz]);
tminpt=gminpt-minpt+1;
trange=max([gmaxpt-gminpt+1 0]);
y(tminpt:tminpt-1+trange,:)=x(gminpt:gmaxpt,:);

% create divisors to average the sum
divisors=ones([sz1 sz])*sw;
switch lower(option.EDGE)
    case 'truncate'
        % find where window is out of data range
        for j=1:nw
            gminw=min([max([-window(j) 0]) sz1]);
            gmaxw=max([sz1-max([window(j) 0])+1 1]);
            divisors(1:gminw,:)=divisors(1:gminw,:)-weights(j);
            divisors(gmaxw:sz1,:)=divisors(gmaxw:sz1,:)-weights(j);
        end
end

% filter to avoid loop (does a running selective weighted sum)
y=filter(lw,1,y,[],option.DIM);

% trim off excess, avg, permute
y=ipermute(y(end-sz1+1:end,:)./divisors,perm);

end
