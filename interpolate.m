function [data]=interpolate(data,sr,method,new_b,new_e)
%INTERPOLATE    Interpolates SAClab data records to a new samplerate
%
%    Description: INTERPOLATE(DATA,RATE) interpolates SAClab data records
%     in DATA to a new sample rate RATE.  As this is interpolation (the
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
%    Header changes: DELTA, NPTS, LEVEN, B, E, DEPMEN, DEPMIN, DEPMAX
%
%    Usage: data=interpolate(data,dt)
%           data=interpolate(data,dt,method)
%           data=interpolate(data,dt,method,new_b,new_e)
%
%    Examples:
%     interpolate records at 5 sps
%      data=intrpolate(data,5);  
%
%     interpolate records at 1 sps from 300 seconds to e
%      data=interpolate(data,1,[],300)
%
%     interpolate at 5 sps from 900 to 950 seconds using linear interp
%      data_pdiff=interpolate(data,5,'linear',900,950)
%
%    See also: syncsr, deci, stretch, iirfilter

%     Version History:
%        ????????????? - Initial Version
%        June 15, 2008 - Documentation update, name change
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 15, 2008 at 00:20 GMT

% check number of arguments
error(nargchk(2,5,nargin))

% check data structure
error(seischk(data,'x'))

% get timing info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
[b,e,npts,delta]=gh(data,'b','e','npts','delta');

% defaults
if(nargin<5 || isempty(new_e)); new_e=e; end
if(nargin<4 || isempty(new_b)); new_b=b; end
if(nargin<3 || isempty(method)); method{1}='spline'; end

% number of records
nrecs=length(data);

% check and expand inputs
if(isscalar(new_b)); new_b(1:nrecs,1)=new_b;
elseif(~isvector(new_b) || length(new_b)~=nrecs)
    error('SAClab:interpolate:badInput','dimensions of new_b are bad')
end
if(isscalar(new_e)); new_e(1:nrecs,1)=new_e;
elseif(~isvector(new_e) || length(new_e)~=nrecs)
    error('SAClab:interpolate:badInput','dimensions of new_e are bad') 
end
if(isscalar(sr)); sr(1:nrecs,1)=sr;
elseif(~isvector(sr) || length(sr)~=nrecs)
    error('SAClab:interpolate:badInput','dimensions of sr are bad') 
end
if(ischar(method)); method=cellstr(method); end
if(~iscellstr(method))
    error('SAClab:interpolate:badInput','method must be char/cellstr')
end
if(isscalar(method)); method=method(ones(nrecs,1),1);
elseif(~isvector(method) || length(method)~=nrecs)
    error('SAClab:interpolate:badInput','dimensions of method are bad')
end

% sampling interval
dt=1./sr;

% looping for each file
for i=1:nrecs
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    
    % make new timing array
    nt=(new_b(i):dt(i):new_e(i)).';
    
    % old timing of data
    if(strcmp(leven(i),'true')); ot=b(i)+(0:npts(i)-1).'*delta(i);
    else ot=data(i).t; data(i).t=[]; end
    
    % interpolate and convert class back
    data(i).x=oclass(interp1(double(ot),double(data(i).x),double(nt),method{i}));
    
    % update header
    data(i)=ch(data(i),'delta',dt(i),'b',nt(1),'e',nt(end),...
        'npts',length(nt),'depmin',min(data(i).x(:)),'leven','t',...
        'depmax',max(data(i).x(:)),'depmen',mean(data(i).x(:)));
end

end
