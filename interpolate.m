function [data]=interpolate(data,sr,method,new_b,new_e)
%INTERPOLATE    Interpolate SEIZMO records to a new samplerate
%
%    Description: INTERPOLATE(DATA,RATE) interpolates SEIZMO records in
%     DATA to a new sample rate RATE.  As this is interpolation (the
%     default method is spline), edge effects are not an issue as they are
%     for SYNCSR, DECI, and STRETCH.  RATE can be a vector of rates with 
%     one element per record in DATA to interpolate records to different 
%     rates.
%
%     INTERPOLATE(DATA,RATE,METHOD) allows selection of the interpolation
%     method from one of the following: 'nearest' (Nearest Neighbor),
%     'linear', 'spline' (Cubic Spline), or 'pchip' (Piecewise Cubic
%     Hermite).  Default is 'spline'.  Method can also be a list of methods
%     to specify a different interpolation method for each record.
%
%     INTERPOLATE(DATA,RATE,METHOD,TIME_START,TIME_END) specifies the time
%     window for interpolation.  The window can be a vector list to specify
%     a separate window for each record.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: DELTA, NPTS, LEVEN, B, E, DEPMEN, DEPMIN, DEPMAX
%
%    Usage:    data=interpolate(data,dt)
%              data=interpolate(data,dt,method)
%              data=interpolate(data,dt,method,new_b,new_e)
%
%    Examples:
%     interpolate records at 5 sps
%      data=interpolate(data,5);  
%
%     interpolate records at 1 sps from 300 seconds to e
%      data=interpolate(data,1,[],300)
%
%     interpolate at 5 sps from 900 to 950 seconds using linear interp
%      data_pdiff=interpolate(data,5,'linear',900,950)
%
%    See also: syncrates, subsample, oversample, iirfilter

%     Version History:
%        Oct. 31, 2007 - initial version
%        Feb. 16, 2008 - doc update, code cleaning
%        Feb. 29, 2008 - better checks
%        Mar.  4, 2008 - better checks, class support
%        Mar. 20, 2008 - change input order
%        May  12, 2008 - fix dep* formula
%        June 15, 2008 - doc update, name changed from INTRPOL8 to
%                        INTERPOLATE
%        Nov. 22, 2008 - better checks, .dep & .ind rather than .x & .t,
%                        doc update, history fix, one CHANGEHEADER call,
%                        extrapolation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 22, 2008 at 20:40 GMT

% check number of arguments
error(nargchk(2,5,nargin))

% check data structure
error(seizmocheck(data,'dep'))

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% get timing info
leven=getlgc(data,'leven');
[b,e,npts,delta]=getheader(data,'b','e','npts','delta');

% defaults
if(nargin<5 || isempty(new_e)); new_e=e; end
if(nargin<4 || isempty(new_b)); new_b=b; end
if(nargin<3 || isempty(method)); method{1}='spline'; end

% number of records
nrecs=numel(data);

% check and expand inputs
if(~isnumeric(new_b))
    error('seizmo:interpolate:badInput',...
        'NEW_B must be numeric!');
elseif(isscalar(new_b))
    new_b(1:nrecs,1)=new_b;
elseif(numel(new_b)~=nrecs)
    error('seizmo:interpolate:badInput',...
        'NEW_B must be scalar or have the same num of elements as DATA!');
end
if(~isnumeric(new_e))
    error('seizmo:interpolate:badInput',...
        'NEW_E must be numeric!');
elseif(isscalar(new_e))
    new_e(1:nrecs,1)=new_e;
elseif(numel(new_e)~=nrecs)
    error('seizmo:interpolate:badInput',...
        'NEW_E must be scalar or have the same num of elements as DATA!');
end
if(~isnumeric(sr))
    error('seizmo:interpolate:badInput',...
        'SR must be numeric!');
elseif(isscalar(sr))
    sr(1:nrecs,1)=sr;
elseif(numel(sr)~=nrecs)
    error('seizmo:interpolate:badInput',...
        'SR must be scalar or have the same num of elements as DATA!');
end
if(ischar(method)); method=cellstr(method); end
if(~iscellstr(method))
    error('seizmo:interpolate:badInput','METHOD must be char/cellstr!')
end
if(isscalar(method))
    method(1:nrecs,1)=method;
elseif(numel(method)~=nrecs)
    error('seizmo:interpolate:badInput',...
        'METHOD must be scalar or have the same num of elements as DATA!');
end

% sampling interval
dt=1./sr;

% looping for each file
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    
    % old timing of data
    if(strcmp(leven(i),'true')); ot=b(i)+(0:npts(i)-1).'*delta(i);
    else ot=data(i).ind; data(i).ind=[]; end
    
    % make new timing array
    nt=(new_b(i):dt(i):new_e(i)).';
    
    % interpolate and convert class back
    data(i).dep=oclass(interp1(...
        double(ot),double(data(i).dep),double(nt),method{i},'extrap'));
    
    % get values (handling dataless)
    npts(i)=numel(nt);
    if(npts(i)==0)
        b(i)=nan;
        e(i)=nan;
    else
        b(i)=nt(1);
        e(i)=nt(end);
        depmen(i)=mean(data(i).dep(:)); 
        depmin(i)=min(data(i).dep(:)); 
        depmax(i)=max(data(i).dep(:));
    end
end

% update header
data=changeheader(data,'delta',dt,'b',b,'e',e,'npts',npts,...
    'leven',true,'depmin',depmin,'depmax',depmax,'depmen',depmen);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
