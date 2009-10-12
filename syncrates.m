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
%    Header Changes: DELTA, NPTS, DEPMEN, DEPMIN, DEPMAX, E
%
%    Examples:
%     Change all records to 5 samples per second:
%      data=syncrates(data,5)
%
%    See also: INTERPOLATE, IIRFILTER, SQUISH, STRETCH

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
%        June 25, 2009 - update for RECORD2MAT/MAT2RECORD, process records
%                        individually, minor bug fix for rare case when
%                        resample does not work
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:30 GMT

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

% require non-xyz records
iftype=getenumdesc(data,'iftype');
if(any(strcmp(iftype,'General XYZ (3-D) file')))
    error('seizmo:syncrates:illegalOperation',...
        'Illegal operation on xyz files!');
end

% get header info
[delta,b,npts]=getheader(data,'delta','b','npts');

% find fraction numerator/denominator of sampling
% rate ratio expressed as integers
[n,d]=rat(delta*sr);

% loop over every record
nrecs=numel(data);
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % save class
    oclass=str2func(class(data(i).dep));
    
    % try resample
    try
        data(i).dep=oclass(resample(double(data(i).dep),n(i),d(i)));
    catch
        % interpolate to a sample rate at least 4 times the current
        % AND the desired sample rates to preserve the data and allow
        % for decimation.  here we find the lowest multiple of the
        % desired rate that is 4x the current rate...hopefully that
        % integer is fairly factorable...if not the function crashes
        lm=max([4 ceil(4/(delta(i)*sr))]);
        
        % interpolate temporarily to higher rate
        data(i)=interpolate(data(i),lm*sr,'spline');
        
        % now resample
        data(i).dep=oclass(resample(double(data(i).dep),1,lm));
    end
    
    % get dep*
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header
npts=floor((npts-1).*n./d)+1; npts(npts<0)=0;
data=changeheader(data,'delta',1/sr,'npts',npts,'e',b+(npts-1)./sr,...
    'depmen',depmen,'depmin',depmin,'depmax',depmax);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
