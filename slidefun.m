function [data]=slidefun(data,fun,nsamples,varargin)
%SLIDEFUN    Apply a sliding window function to SAClab data records
%
%    Description: SLIDEFUN(DATA,FUN,N) applies the function defined by the
%     function handle FUN to a centered sliding window of 2N+1 samples to 
%     the dependent component(s) of SAClab data records in DATA.  FUN is 
%     expected to handle a column vector (single component) or an array 
%     (multiple component - components are distributed in columns).  Output
%     of FUN is assigned to the reference data point of the sliding window.  
%     Vector output will be distributed as multiple components.  Array 
%     output will produce an error.  N can be a scalar (each record has the
%     same window size) or a vector (define each record's window size 
%     separately).
%
%     SLIDEFUN(...,'POSITION','CENTER'|'TRAIL'|'LEAD') sets the position of
%     the sliding window relative to the reference data point.  CENTER 
%     positions the window such that the reference point is at its center.
%     TRAIL positions the window to trail the reference point such that the
%     reference point has the highest index in the window.  LEAD sets the
%     window to lead the reference point such that the reference point has 
%     the lowest index in the window.  Note that the window size for a
%     CENTER positioning is 2N+1, while for TRAIL or LEAD the window size
%     is N.  Default position is CENTER.
%     
%     SLIDEFUN(...,'OFFSET',OFFSET) sets the offset of the sliding window
%     from to the reference data point in number of samples.  For a
%     centered window (see option 'POSITION') this introduces a gap of 
%     2*OFFSET-1 in the window.  For example an OFFSET of 1 will exclude
%     the reference data point from the sliding window.  Negative OFFSETS
%     are allowed for a centered window, but they are complicated due to
%     overlap.  For example an OFFSET of -1 will make the window include
%     the reference data point twice and an OFFSET of -2 will cause the 3
%     centermost points to be included twice and so on.  Default OFFSET=0.
%     
%     SLIDEFUN(...,'EDGE','TRUNCATE'|'PAD') sets how SLIDEFUN handles the
%     edge cases when applying the sliding window function.  TRUNCATE
%     eliminates points in the window that do not reference a datapoint
%     (for instance if a window extends before or after the data, that
%     portion of the window will be truncated).  PAD adds zeros to the data
%     so that all the points of the sliding window always reference some 
%     value.  Default setting is TRUNCATE.
%
%    Notes:
%     - The number of components in the output record need not match that
%       of the input record
%     - Centered windows are of length 2N+1, while the others are just N
%
%    System requirements: Matlab 7
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, NCMP
%
%    Usage:    data=slidefun(data,fun,nsamples)
%              data=slidefun(...,'position','center'|'trail'|'lead')
%              data=slidefun(...,'offset',offset)
%              data=slidefun(...,'edge','truncate'|'pad')
%
%    Examples:
%     Running absolute mean normalization like G. D. Bensen et al, 2007 
%     where Tmin, Tmax bound the period band of the regional seismicity 
%     to be downweighted (usually 5-150s):
%      f=iirfilter(data,'bandpass','butter',[1/Tmax 1/Tmin],4,2);
%      w=slidefun(f,@(x)mean(abs(x)),ceil(Tmax./(4*gh(data,'delta'))));
%      time_norm_data=divf(data,w);
%
%     Root-Mean-Square:
%      slidefun(data,@(x)sqrt(mean(x.^2)),N)
%
%     Root-Median-Square:
%      slidefun(data,@(x)sqrt(median(x.^2)),N)
%
%    See also: seisfun

%     Version History:
%        Apr. 23, 2008 - initial version
%        May  12, 2008 - dep* fix
%        July 18, 2008 - documentation update, dataless support, .dep
%                        rather than .x, added history, single ch call
%        Oct.  1, 2008 - big changes: now several options to control the
%                        window size and position as well as how edge cases
%                        are handled
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  1, 2008 at 21:55 GMT

% todo:

% check nargin
error(nargchk(3,9,nargin))

% check data structure
error(seischk(data,'dep'))

% check input fun is a function
if(~isa(fun,'function_handle'))
    error('SAClab:slidefun:badInput','FUN must be a function handle!')
end

% number of records
nrecs=numel(data);

% fix/check nsamples
if(isscalar(nsamples)); nsamples=nsamples(ones(1,nrecs),1); end
if(~isnumeric(nsamples) || any(nsamples<1) || numel(nsamples)~=nrecs ...
        || any(fix(nsamples)~=nsamples))
    error('SAClab:slidefun:badInput',...
        'NSAMPLES must be a scalar/vector of positive integer(s)')
end

% option defaults
option.POSITION='center';
option.OFFSET=0;
option.EDGE='truncate';

% get options set by SACLAB global
global SACLAB; fields=fieldnames(option).';
if(isfield(SACLAB,'SLIDEFUN'))
    for i=fields
        if(isfield(SACLAB.SLIDEFUN,i{:})); 
            option.(i{:})=SACLAB.SLIDEFUN.(i{:}); 
        end
    end
end

% options must be field-value pairs
nargopt=length(varargin);
if(mod(nargopt,2))
    error('SAClab:slidefun:badNumOptions','Unpaired option!')
end

% get inline options
for i=1:2:nargopt
    varargin{i}=upper(varargin{i});
    if(isfield(option,varargin{i}))
        option.(varargin{i})=varargin{i+1};
    else
        warning('SAClab:slidefun:badInput',...
            'Unknown Option: %s !',varargin{i}); 
    end
end

% check options
if(~ischar(option.POSITION) || ...
        ~any(strcmpi(option.POSITION,{'center' 'trail' 'lead'})))
    error('SAClab:slidefun:badOptionValue',...
        'POSITION option must be ''CENTER'', ''TRAIL'' or ''LEAD''!');
end
if(~isscalar(option.OFFSET) || ...
        option.OFFSET~=fix(option.OFFSET))
    error('SAClab:slidefun:badOptionValue',...
        'OFFSET option must be an integer!');
end
if(~ischar(option.EDGE) || ...
        ~any(strcmpi(option.EDGE,{'pad' 'truncate'})))
    error('SAClab:slidefun:badOptionValue',...
        'EDGE option must be ''PAD'' or ''TRUNCATE''');
end

% loop through each record
ncmp=nan(nrecs,1); depmen=ncmp; depmin=ncmp; depmax=ncmp;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % get storage class of data
    oclass=str2func(class(data(i).dep));
    
    % build window indices
    switch option.POSITION
        case 'center'
            hw=(1:nsamples(i))+option.OFFSET;
            window=[-fliplr(hw) zeros(~logical(option.OFFSET),1) hw];
        case 'trail'
            window=(1-nsamples(i):0)+option.OFFSET;
        case 'lead'
            window=(0:nsamples(i)-1)+option.OFFSET;
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
            window=window+abs(min(window));
            data(i).dep=[zeros(abs(min(window)),sz(2)); ...
                        double(data(i).dep); ...
                        zeros(abs(max(window)),sz(2))];
            % slide function across
            for j=1:sz(1)
                temp(j,:)=fun(data(i).dep(j+window,:));
            end
        case 'truncate'
            % slide function across
            for j=1:sz(1)
                % allow only valid indices
                window=j+window;
                window=window(window>=1 & window<=sz(1));
                temp(j,:)=fun(data(i).dep(window,:));
            end
    end
    
    % assign back to data and change storage back
    data(i).dep=oclass(temp);
    
    % dep*
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
    ncmp(i)=size(data(i).dep,2);
end

% update header
warning('off','SAClab:ch:fieldInvalid');
data=ch(data,'depmen',depmen,'depmin',depmin,'depmax',depmax,'ncmp',ncmp);
warning('on','SAClab:ch:fieldInvalid');

end
