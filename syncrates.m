function [data]=syncrates(data,sr)
%SYNCRATES    Resample SEIZMO records to a common sample rate
%
%    Usage:    data=syncrates(data,sr)
%
%    Description: SYNCRATES(DATA,SR) syncronizes the sample rates of SEIZMO 
%     records in DATA to the sample rate SR.  A fir filter is used to
%     avoid aliasing issues, but this can cause edge effects if the records
%     deviate from zero strongly at the start/end of the record.  Typically
%     using REMOVETREND and TAPER on records beforehand helps to limit the
%     edge effects.  Uses the Matlab function resample (Signal Processing
%     Toolbox).
%
%    Notes:
%     - requires evenly sampled data (use INTERPOLATE for uneven data)
%
%    Tested on: Matlab r2007b
%
%    Header Changes: DELTA, NPTS, DEPMEN, DEPMIN, DEPMAX, E
%
%    Examples:
%     Change all records to 5 samples per second:
%      data=syncrates(data,5)
%
%    See also: interpolate, iirfilter, squish, stretch

%     Version History:
%        Feb. 16, 2008 - initial version
%        Feb. 23, 2008 - minor fix
%        Feb. 26, 2008 - minor fix
%        Mar.  4, 2008 - doc update, better checks
%        June 16, 2008 - doc update, rejects spectral records
%                        fixed no header update bug, and changed the way
%                        records are resampled when Matlab's resample fails
%        Nov. 22, 2008 - update for new name schema (now SYNCRATES), allow
%                        resampling spectral records, disallow resampling
%                        xyz records
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 21:05 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% check rate
if(~isnumeric(sr) || ~isscalar(sr) || sr<=0)
    error('seizmo:syncrates:badInput',...
        'SR must be a positive numeric scalar!');
end

% require evenly sampled records only
if(any(~strcmpi(getlgc(data,'leven'),'true')))
    error('seizmo:syncrates:evenlySpacedOnly',...
        'Illegal operation on unevenly spaced data!');
end

% require non-spectral records
iftype=getenumdesc(data,'iftype');
if(any(strcmp(iftype,'General XYZ (3-D) file')))
    error('seizmo:syncrates:illegalOperation',...
        'Illegal operation on xyz files!');
end

% make delta groups
delta=getheader(data,'delta');
gdt=unique(delta); ng=length(gdt); gi=cell(ng,1);
for i=1:ng; gi{i}=find(delta==gdt(i)).'; end

% find fraction numerator/denominator of sampling
% rate ratio expressed as integers
[n,d]=rat(gdt*sr);

% work on each group
for i=1:ng
    % try resample
    try
        data(gi{i})=chsr(data(gi{i}),n(i),d(i),sr);
    catch
        % interpolate to a sample rate at least 4 times the current
        % AND the desired sample rates to preserve the data and allow
        % for decimation.  here we find the lowest multiple of the
        % desired rate that is 4x the current rate...hopefully that
        % integer is fairly factorable...if not the function crashes
        lowest_multiple=max([4 ceil(4/(gdt*sr))]);
        
        % interpolate temporarily to higher rate
        data(gi{i})=interpolate(data(gi{i}),lowest_multiple*sr,'spline');
        
        % now resample
        data(gi{i})=chsr(data(gi{i}),1,lowest_multiple,sr);
    end
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end

function [data]=chsr(data,n,d,sr)
%CHSR  Resample records to new rate

% combine records
[dep,idx1,ind,idx2,store,npts]=combinerecords(data);

% resample
dep=resample(dep,n,d);

% new length
nnpts=floor((npts-1)*n/d)+1;

% redistribute records
e=getheader(data,'b')+(nnpts-1)./sr;
data=changeheader(data,'delta',1/sr,'npts',nnpts,'e',e);
data=distributerecords(data,dep,idx1,[],[],store,nnpts);

end
