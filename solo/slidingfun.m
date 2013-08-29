function [data]=slidingfun(data,fun,nsamples,varargin)
%SLIDINGFUN    Apply a sliding window function to SEIZMO records
%
%    Usage:    data=slidingfun(data,fun,n)
%              data=slidingfun(...,'position','center'|'trail'|'lead')
%              data=slidingfun(...,'offset',offset)
%              data=slidingfun(...,'edge','truncate'|'pad')
%
%    Description:
%     SLIDINGFUN(DATA,FUN,N) applies the function defined by the function
%     handle FUN to a centered sliding window of 2N+1 samples to the
%     dependent component(s) of SEIZMO records in DATA.  FUN is expected to
%     handle a column vector (single component) or an array (multiple
%     component - components are distributed in columns).  Output of FUN is
%     assigned to the reference data point of the sliding window.  Vector
%     output will be distributed as multiple components.  Array output will
%     produce an error.  N can be a scalar (each record has the same window
%     size) or a vector (define each record's window size separately).
%
%     SLIDINGFUN(...,'POSITION','CENTER'|'TRAIL'|'LEAD') sets the position
%     of the sliding window relative to the reference data point.  CENTER 
%     positions the window such that the reference point is at its center.
%     TRAIL positions the window to trail the reference point such that the
%     reference point has the highest index in the window.  LEAD sets the
%     window to lead the reference point such that the reference point has 
%     the lowest index in the window.  Note that the window size for a
%     CENTER positioning is 2N+1, while for TRAIL or LEAD the window size
%     is N.  Default position is CENTER.
%     
%     SLIDINGFUN(...,'OFFSET',OFFSET) sets the offset of the sliding window
%     from to the reference data point in number of samples.  For a
%     centered window (see option 'POSITION') this introduces a gap of 
%     2*OFFSET-1 in the window.  For example an OFFSET of 1 will exclude
%     the reference data point from the sliding window.  Negative OFFSETS
%     are allowed for a centered window, but they are complicated due to
%     overlap.  For example an OFFSET of -1 will make the window include
%     the reference data point twice and an OFFSET of -2 will cause the 3
%     centermost points to be included twice and so on.  Default OFFSET=0.
%     
%     SLIDINGFUN(...,'EDGE','TRUNCATE'|'PAD') sets how to handle the
%     edge cases when applying the sliding window function.  TRUNCATE
%     eliminates points in the window that do not reference a datapoint
%     (for instance if a window extends before or after the data, that
%     portion of the window will be truncated).  PAD adds zeros to the data
%     so that all the points of the sliding window always reference some 
%     value.  Default setting is TRUNCATE.
%
%    Notes:
%     - The number of components in the output record need not match that
%       of the input record.
%     - Centered windows are of length 2N+1, while the others are just N.
%     - SLIDINGFUN is much _slooower_ than SLIDINGABSMEAN or SLIDINGRMS.
%       The main value is that it allows any function to be slid, not just
%       an average.  For example a median.
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, NCMP
%
%    Examples:
%     % Running absolute mean normalization like G. D. Bensen et al, 2007 
%     % where Tmin, Tmax bound the period band of the regional seismicity 
%     % to be downweighted (usually 5-150s):
%     f=iirfilter(data,'bandpass','butter',[1/Tmax 1/Tmin],4,2);
%     w=slidingfun(f,@(x)mean(abs(x)),...
%         ceil(Tmax./(4*getheader(data,'delta'))));
%     time_norm_data=dividerecords(data,w);
%
%     % Sliding Root-Mean-Square:
%     slidingfun(data,@(x)sqrt(mean(x.^2)),N)
%
%     % Sliding Root-Median-Square:
%     slidingfun(data,@(x)sqrt(median(x.^2)),N)
%
%    See also: SOLOFUN, SLIDINGRMS, SLIDINGABSMEAN, SLIDINGAVG

%     Version History:
%        Apr. 23, 2008 - initial version
%        May  12, 2008 - dep* fix
%        July 18, 2008 - doc update, dataless support, .dep
%                        rather than .x, added history, single ch call
%        Oct.  1, 2008 - big changes: now several options to control the
%                        window size and position as well as how edge cases
%                        are handled
%        Oct.  3, 2008 - fixed bug related to incrementing wrong variable
%        Oct.  5, 2008 - fixed bug that caused pad not to, fixed bug that
%                        did not change class to double for the function,
%                        allow OFFSET option to have value for each record
%        Nov. 22, 2008 - update for new name schema (now SLIDINGFUN)
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Oct.  6, 2009 - dropped use of LOGICAL function
%        Oct. 15, 2009 - increased nargin allowance
%        Feb.  3, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        versioninfo caching
%        Jan.  6, 2011 - drop versioninfo caching, nargchk fix,
%                        seizmofun/solofun rename
%        Jan. 28, 2012 - doc update, drop SEIZMO global
%        Aug.  2, 2013 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  2, 2013 at 09:45 GMT

% todo:

% check nargin
error(nargchk(3,inf,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% check input fun is a function
if(~isa(fun,'function_handle'))
    error('seizmo:slidingfun:badInput',...
        'FUN must be a function handle!')
end

% verbosity
verbose=seizmoverbose;

% number of records
nrecs=numel(data);

% fix/check nsamples
if(isscalar(nsamples)); nsamples=nsamples(ones(1,nrecs),1); end
if(~isnumeric(nsamples) || any(nsamples<1) || numel(nsamples)~=nrecs ...
        || any(fix(nsamples)~=nsamples))
    error('seizmo:slidingfun:badInput',...
        'N must be positive integer(s)')
end

% option defaults
option.POSITION='center';
option.OFFSET=0;
option.EDGE='truncate';

% options must be field-value pairs
nargopt=numel(varargin);
if(mod(nargopt,2))
    error('seizmo:slidingfun:badNumOptions','Unpaired option!')
end

% get inline options
for i=1:2:nargopt
    varargin{i}=upper(varargin{i});
    if(isfield(option,varargin{i}))
        option.(varargin{i})=varargin{i+1};
    else
        warning('seizmo:slidingfun:badInput',...
            'Unknown Option: %s !',varargin{i}); 
    end
end

% check options
if(~ischar(option.POSITION) || ...
        ~any(strcmpi(option.POSITION,{'center' 'trail' 'lead'})))
    error('seizmo:slidingfun:badOptionValue',...
        'POSITION option must be ''CENTER'', ''TRAIL'' or ''LEAD''!');
end
if(isscalar(option.OFFSET))
    option.OFFSET=option.OFFSET(ones(1,nrecs),1); 
end
if(~isnumeric(option.OFFSET) || numel(option.OFFSET)~=nrecs || ...
        any(option.OFFSET~=fix(option.OFFSET)))
    error('seizmo:slidingfun:badOptionValue',...
        'OFFSET option must be integer(s)!');
end
if(~ischar(option.EDGE) || ...
        ~any(strcmpi(option.EDGE,{'pad' 'truncate'})))
    error('seizmo:slidingfun:badOptionValue',...
        'EDGE option must be ''PAD'' or ''TRUNCATE''');
end

% detail message
if(verbose)
    disp('Applying Sliding Function to Record(s)');
    print_time_left(0,nrecs);
end

% loop through each record
ncmp=nan(nrecs,1); depmen=ncmp; depmin=ncmp; depmax=ncmp;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep))
        % detail message
        if(verbose); print_time_left(i,nrecs); end
        continue;
    end
    
    % get storage class of data
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % build window indices
    switch option.POSITION
        case 'center'
            hw=(1:nsamples(i))+option.OFFSET(i);
            window=[-fliplr(hw) zeros(option.OFFSET(i)==0,1) hw];
        case 'trail'
            window=(1-nsamples(i):0)+option.OFFSET(i);
        case 'lead'
            window=(0:nsamples(i)-1)+option.OFFSET(i);
    end
    
    % testbed (get the size of output)
    sz=size(data(i).dep);
    temp=fun(data(i).dep(1:numel(window),:));
    ncols=numel(temp);
    temp=nan(sz(1),ncols);
    
    % change method based on how we deal with edge cases
    switch option.EDGE
        case 'pad'
            % pad data and adjust window
            pwindow=window+abs(min(window));
            data(i).dep=[zeros(abs(min(window)),sz(2)); ...
                        data(i).dep; ...
                        zeros(abs(max(window)),sz(2))];
            % slide function across
            for j=1:sz(1)
                temp(j,:)=fun(data(i).dep(j+pwindow,:));
            end
        case 'truncate'
            % slide function across
            for j=1:sz(1)
                % allow only valid indices
                jwindow=j+window;
                tjwindow=jwindow(jwindow>=1 & jwindow<=sz(1));
                temp(j,:)=fun(data(i).dep(tjwindow,:));
            end
    end
    
    % assign back to data and change storage back
    data(i).dep=oclass(temp);
    
    % dep*
    depmen(i)=nanmean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
    ncmp(i)=size(data(i).dep,2);

    % detail message
    if(verbose); print_time_left(i,nrecs); end
end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt convolution
try
    % update header
    data=changeheader(data,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax,'ncmp',ncmp);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
