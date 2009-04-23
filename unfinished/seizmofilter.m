function [data]=seizmofilter(data,fs,idx,passes,mirror)
%SEIZMOFILTER    Implements a filter on SEIZMO records
%


% check number of input arguments
msg=nargchk(2,5,nargin);
if(~isempty(msg)); error(msg); end;

% check data structure
error(seizmocheck(data,'dep'))
nrecs=numel(data);

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check filter structure
nfilters=numel(fs);
if(iscell(fs))
    for i=1:nfilters
        if(~strmatch(class(fs{i}),'dfilt.'))
            error('seizmo:seizmofilter:badFilter',...
                'FS must be a cell array of dfilt objects!');
        end
    end
elseif(~strmatch(class(fs),'dfilt.'))
    error('seizmo:seizmofilter:badFilter',...
        'FS must be a dfilt object!');
end

% check indices
if(nargin<3 || isempty(idx))
    idx=ones(nrecs,1);
elseif(~nrecs==numel(idx) || any(idx~=fix(idx)) ...
        || any(idx<1) || any(idx>nfilters))
    error('seizmo:seizmofilter:badIndices',...
        'IDX must be a list giving the filter index for each record!');
end

% check passes
if(nargin<4 || isempty(passes))
    passes=1;
elseif(~isscalar(passes) || ~isreal(passes) || passes~=fix(passes) ...
        || passes<1 || passes>4)
    error('seizmo:seizmofilter:badNumPasses',...
        'PASSES must be 1, 2, 3, or 4!');
end

% check mirrorflip
if(nargin<5 || isempty(mirror))
    mirror=true;
elseif(~isscalar(mirror) || ~islogical(mirror))
    error('seizmo:seizmofilter:badMirrorFlip',...
        'MIRRORFLIP must be a scalar logical!');
end

% loop over each filter
for i=1:nfilters
    % combine records
    [recs,idx1,ind,idx2,store,npts]=combinerecords(data(idx==i));

    % implement
    if(passes==1)
        % single pass - forward
        recs=impfilt(recs,fs{i},mirror,true);
    elseif(passes==2)
        % double pass - forward then backward
        recs=impfilt(impfilt(recs,fs{i},mirror,true),fs{i},mirror,false);
    elseif(passes==3)
        % single pass - backward
        recs=impfilt(recs,fs{i},mirror,false);
    elseif(passes==4)
        % double pass - backward then forward
        recs=impfilt(impfilt(recs,fs{i},mirror,false),fs{i},mirror,true);
    end
    
    % distribute records back
    data(idx==i)=distributerecords(data(idx==i),recs,idx1,[],[],store,npts);
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end

function [recs]=impfilt(recs,fs,mirror,forward)
%IMPFILT   Implements filter
%   Implements a filter design on SEIZMO data records and makes appropriate
%   header updates.  Takes a mirror option which does pseudo-IC to limit 
%   edge effects.  Works with multiple records.

% reverse
if(~forward); recs=recs(end:-1:1,:); end

% mirror logical
if(mirror)
    % prepend a mirror-flip of the series to limit edge effects
    recs=filter(fs,[2*recs(ones(end-1,1),:)-recs(end:-1:2,:); recs]);
    recs=recs(ceil(end/2):end,:);
else
    % straight forward filter
    recs=filter(fs,recs);
end

% reverse
if(~forward); recs=recs(end:-1:1,:); end

end
