function [data]=syncrates(data,sr,tol)
%SYNCRATES    Resample SEIZMO records to a common sample rate
%
%    Usage:    data=syncrates(data,sr)
%              data=syncrates(data,sr,tol)
%
%    Description: SYNCRATES(DATA,SR) syncronizes the sample rates of SEIZMO 
%     records in DATA to the sample rate SR.  A fir filter is used to
%     avoid aliasing issues, but this can cause edge effects if the records
%     deviate from zero strongly at the start/end of the record.  Typically
%     using REMOVETREND and TAPER on records beforehand helps to limit the
%     edge effects.  Uses the Matlab function RESAMPLE (Signal Processing
%     Toolbox) and RAT (see Notes below!).
%
%     SYNCRATES(DATA,SR,TOL) specifies the maximum tolerance TOL that the
%     fraction of 2 small integers must match the ratio of the old and new
%     sample rates of a record.  The integers specify the upsampling and
%     downsampling portions of the resampling operation.  See function RAT
%     for more details.  The default TOL is 1e-6:
%       abs(Up/Down - New/Old) / (New/Old) <= 1e-6
%
%    Notes:
%     - requires evenly sampled data (use INTERPOLATE for uneven data)
%     - Matlab r2007b function RAT has a bug in it - change line 116ish to:
%        if(x==0) || (abs((C(1,1)/C(2,1)-X(j))/X(j))<=max(tol,eps(X(j))))
%       This will force RAT to function as described.
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
%        Nov. 26, 2009 - document RAT bug, alter RAT call slightly to force
%                        better accuracy of the resampling operation, add
%                        TOL argument, fix NPTS handling
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 26, 2009 at 15:20 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
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

% check tol
if(nargin==2 || isempty(tol)); tol=1e-6; end
if(~isscalar(tol) || ~isreal(tol))
    error('seizmo:syncrates:badInput','TOL must be a scalar real!');
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
[delta,b]=getheader(data,'delta','b');

% find fraction numerator/denominator of sampling
% rate ratio expressed as integers
[n,d]=rat(delta*sr,tol);

% loop over every record
nrecs=numel(data);
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen; npts=depmen;
for i=1:nrecs
    % skip dataless
    if(isempty(data(i).dep)); continue; end
    
    % save class
    oclass=str2func(class(data(i).dep));
    
    % try resample
    try
        data(i).dep=oclass(resample(double(data(i).dep),n(i),d(i)));
    catch
        disp(['Had some trouble resampling record: ' num2str(i)])
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
    
    % get npts
    npts(i)=size(data(i).dep,1);
    
    % get dep*
    if(npts(i))
        depmen(i)=mean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
    end
end

% update header
data=changeheader(data,'delta',1/sr,'npts',npts,'e',b+(npts-1)./sr,...
    'depmen',depmen,'depmin',depmin,'depmax',depmax);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
