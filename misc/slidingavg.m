function [y]=slidingavg(x,nsamples,varargin)
%SLIDINGAVG    Returns sliding-window average of data
%
%    Usage:    y=slidingavg(x,n)
%              y=slidingavg(...,'position','center'|'trail'|'lead')
%              y=slidingavg(...,'offset',offset)
%              y=slidingavg(...,'edge','truncate'|'pad')
%              y=slidingavg(...,'nans','skipped'|'too'|'only')
%              y=slidingavg(...,'dim',n)
%              y=slidingavg(...,'custom',window)
%
%    Description:
%     Y=SLIDINGAVG(X,N) applies a centered sliding-window average of 2N+1
%     samples down the columns of numeric array X.  Sliding windows
%     extending outside the record are truncated (look at the 'EDGE' option
%     to change this).
%
%     Y=SLIDINGAVG(...,'POSITION','CENTER'|'TRAIL'|'LEAD') sets the layout
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
%     Y=SLIDINGAVG(...,'OFFSET',OFFSET) sets the offset of the sliding
%     window from the reference data point in number of samples.  For the
%     LEAD or TRAIL positioned window (see option 'POSITION') the offset is
%     as expected.  For a centered window, any nonzero offset introduces a
%     gap of 2*OFFSET-1 samples in the center of the window and removes one
%     sample (making the window now 2N samples).  For example, an OFFSET of
%     1 will exclude the reference data point from the sliding window.
%     Negative OFFSET for a centered window are complicated due to overlap.
%     For example an OFFSET of -1 will include the reference data point
%     twice and an OFFSET of -2 will cause the 3 centermost points to be
%     included twice and so on.  Default OFFSET=0.
%     
%     Y=SLIDINGAVG(...,'EDGE','TRUNCATE'|'PAD') sets how to handle the edge
%     cases.  TRUNCATE eliminates points in the sliding window that do not
%     reference a datapoint (for instance if the window extends before or
%     after the data, that portion of the window will be truncated).  PAD
%     adds zeros to the data so that all the points of the sliding-window
%     always reference some value.  Default setting is TRUNCATE.
%
%     Y=SLIDINGAVG(...,'NANS','SKIPPED'|'TOO'|'ONLY') sets the treatment of
%     nans and nonnans in X.  In all cases the nans are eliminated from the
%     sliding window (windows without any nonnans will return nan).  These
%     options just control if the nonnan/nan values are updated with the
%     sliding window average or restored back to their original value after
%     the sliding average.  SKIPPED will restore the nan elements while
%     updating the nonnans.  TOO will update both nonnan and nan values.
%     ONLY updates nan elements while restoring the nonnan elements.  The
%     default is 'TOO'.
%
%     Y=SLIDINGAVG(...,'DIM',N) specifies an alternative dimension to slide
%     across.  Default seeting is [] (slides across the 1st non-singleton
%     dimension or rows otherwise).
%
%     Y=SLIDINGAVG(...,'CUSTOM',WINDOW) allows a custom sliding window
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
%    Examples:
%     % Get a smoothed amplitude spectra:
%     y=slidingavg(abs(fft(x)),2)
%
%     % An example of a custom call:
%     y=slidingavg(x,[],'custom',[-10:10; gausswin(21).'])
%
%    See also: SLIDINGRMS, SLIDINGABSMEAN, SLIDINGFUN

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
%        May  10, 2009 - minor doc fix
%        June  9, 2009 - allow for unlimited varargin
%        Feb. 11, 2011 - mass nargchk fix
%        Jan. 28, 2012 - drop SEIZMO global, doc update
%        May  30, 2012 - allow N=0 for position=center
%        Aug.  1, 2013 - nans option for nan handling, dim option defaults
%                        to 1st non-singleton dimension, major code cleanup
%        Aug.  2, 2013 - fixed dimension bug, fixed other bugs introduced
%                        with recent code changes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  2, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% quick return if empty
if(isempty(x)); y=x; return; end

% check input array
if(~isnumeric(x))
    error('seizmo:slidingavg:badInput','X must be numeric!');
end

% fix/check nsamples
if(~isnumeric(nsamples) || ~isscalar(nsamples) ...
        || nsamples<0 || fix(nsamples)~=nsamples)
    error('seizmo:slidingavg:badInput','N must be a positive integer!');
end

% option defaults
option.POSITION='center';
option.OFFSET=0;
option.EDGE='truncate';
option.NANS='too';
option.DIM=[];
option.CUSTOM=[];

% options must be option-value pairs
nargopt=length(varargin);
if(mod(nargopt,2))
    error('seizmo:slidingavg:badNumOptions','Unpaired option!')
end

% get options
for i=1:2:nargopt
    switch lower(varargin{i})
        case {'position' 'pos' 'p'}
            if(isempty(varargin{i+1})); continue; end
            option.POSITION=varargin{i+1};
        case {'offset' 'off' 'o'}
            if(isempty(varargin{i+1})); continue; end
            option.OFFSET=varargin{i+1};
        case {'edge' 'e'}
            if(isempty(varargin{i+1})); continue; end
            option.EDGE=varargin{i+1};
        case {'nans' 'nan' 'n'}
            if(isempty(varargin{i+1})); continue; end
            option.NANS=varargin{i+1};
        case {'dim' 'dimension' 'd'}
            option.DIM=varargin{i+1};
        case {'custom' 'cust' 'c'}
            option.CUSTOM=varargin{i+1};
        otherwise
            error('seizmo:slidingavg:badInput',...
                'Unknown Option: %s !',varargin{i});
    end
end

% default empty dim option to 1st non-singleton dimension, otherwise 1
sx=size(x);
if(isempty(option.DIM))
    option.DIM=find(sx~=1,1);
    if(isempty(option.DIM)); option.DIM=1; end
end

% check options
csz=size(option.CUSTOM);
if(~ischar(option.POSITION) || size(option.POSITION,1)~=1 ...
        || ndims(option.POSITION)~=2 ...
        || ~any(strcmpi(option.POSITION,{'center' 'trail' 'lead'})))
    error('seizmo:slidingavg:badOptionValue',...
        'POSITION option must be ''CENTER'', ''TRAIL'' or ''LEAD''!');
elseif(~isnumeric(option.OFFSET) || ~isscalar(option.OFFSET) ...
        || option.OFFSET~=fix(option.OFFSET))
    error('seizmo:slidingavg:badOptionValue',...
        'OFFSET option must be an integer!');
elseif(~ischar(option.EDGE) || size(option.EDGE,1)~=1 ...
        || ndims(option.EDGE)~=2 ...
        || ~any(strcmpi(option.EDGE,{'pad' 'truncate'})))
    error('seizmo:slidingavg:badOptionValue',...
        'EDGE option must be ''PAD'' or ''TRUNCATE''!');
elseif(~ischar(option.NANS) || size(option.NANS,1)~=1 ...
        || ndims(option.NANS)~=2 ...
        || ~any(strcmpi(option.NANS,{'skipped' 'too' 'only'})))
    error('seizmo:slidingavg:badOptionValue',...
        'NANS option must be ''SKIPPED'', ''TOO'', or ''ONLY''!');
elseif(~isnumeric(option.DIM) || ~isscalar(option.DIM) || option.DIM<1 ...
        || option.DIM~=fix(option.DIM))
    error('seizmo:slidingavg:badOptionValue',...
        'DIM option must be a positive integer!');
elseif(~isempty(option.CUSTOM) && (~isnumeric(option.CUSTOM) ...
        || numel(csz)>2 || csz(1)~=2 ...
        || any(option.CUSTOM(1,:)~=fix(option.CUSTOM(1,:)))))
    error('seizmo:slidingavg:badOptionValue',...
        'CUSTOM option must be a 2xN numeric array!');
end

% build filter window
if(isempty(option.CUSTOM))
    % build preset window indices
    switch option.POSITION
        case 'center'
            % catch no window case
            if(nsamples==0 && option.OFFSET~=0)
                error('seizmo:slidingavg:badInput',...
                    'N must be >0 for nonzero OFFSET on a CENTER window!');
            end
            hw=(1:nsamples)+option.OFFSET-(option.OFFSET>0);
            window=[-fliplr(hw) zeros(option.OFFSET==0,1) hw];
        case 'trail'
            if(nsamples<1)
                error('seizmo:slidingavg:badInput',...
                    'N must be a positive integer!');
            end
            window=(1-nsamples:0)+option.OFFSET;
        case 'lead'
            if(nsamples<1)
                error('seizmo:slidingavg:badInput',...
                    'N must be a positive integer!');
            end
            window=(0:nsamples-1)+option.OFFSET;
    end
    weights=ones(1,numel(window));
else
    window=option.CUSTOM(1,:);
    weights=option.CUSTOM(2,:);
end

% translate the window to filter input
minw=min(window);
maxw=max(window);
owidx=maxw-window+1;         % note indices are reversed
lw=zeros(1,maxw-minw+1);     % default is points have no weight
%lw(owidx)=lw(owidx)+weights; % weight specified (no multi-entry support)
for i=1:numel(weights); lw(owidx(i))=lw(owidx(i))+weights(i); end
% NOTE: the above loop is necessary to support windows with repeated
%       indices (this only really applies to offset<0 & pos=center)

% x permutation setup
sx(end+1:option.DIM)=1; % handle dimension specified beyond last dimension
perm=1:numel(sx);
perm(option.DIM)=[];
perm=[option.DIM perm];
sx1=sx(option.DIM);
sx(option.DIM)=[];

% bring averaged dimension to front
x=permute(x,perm);

% keep track of nans and set them to 0 to stabilize filter
xnan=isnan(x);
x(xnan)=0;

% pad/cut setup
minpt=1+minw;   % min/max contributing in
maxpt=sx1+maxw; % edge=pad case
gminpt=max([minpt 1]);   % min/max contributing in
gmaxpt=min([maxpt sx1]); % edge=truncate case
tminpt=gminpt-minpt+1;           % first data point after padding
trange=max([gmaxpt-gminpt+1 0]); % number of data points

% padding with zeros and cutting out unused data
y=zeros([maxpt-minpt+1 sx]);
y(tminpt:tminpt-1+trange,:)=x(gminpt:gmaxpt,:);

% compute divisors using the filter on a logical array
% is padding included in divisor sum?
switch lower(option.EDGE)
    case 'pad'      % yes
        z=true([maxpt-minpt+1 sx]);
    case 'truncate' % no
        z=false([maxpt-minpt+1 sx]);
end
z(tminpt:tminpt-1+trange,:)=~xnan(gminpt:gmaxpt,:);
z=filter(lw,1,z,[],1);

% filter to avoid loop (does a running selective weighted sum)
y=filter(lw,1,y,[],1);

% trim off excess, divide to get avg
y=y(end-sx1+1:end,:)./z(end-sx1+1:end,:);

% restore values if desired
switch lower(option.NANS)
    case 'skipped'
        y(xnan)=nan;
    case 'only'
        y(~xnan)=x(~xnan);
end

% permute back
y=ipermute(y,perm);

end
