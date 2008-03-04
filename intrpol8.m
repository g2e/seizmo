function [data]=intrpol8(data,sr,new_b,new_e,method)
%INTRPOL8    Interpolates SAClab data records to a new sampling frequency
%
%    Description: Interpolates data records to new sample rate sr using 
%     Matlab's 'interp1' function.  As this is interpolation, edge effects 
%     are not a problem as they are for syncsr, deci, and stretch.  Other 
%     inputs (optional) are the new start and end times desired (default is 
%     to use the original record start and end times) as well as the method
%     of interpolation (default is piececubic cubic). Note that the new end
%     times of the records cannot be honored if it does not fall on a 
%     sample interval, in which case the last point before the end time 
%     will be the set as the new end time.
%
%    Usage: [data]=intrpol8(data,dt,new_b,new_e,method)
%
%    Examples:
%     interpolate at 5 sps using b and e header values
%      data=intrpol8(data,5);  
%
%     interpolate at 1 sps from 300 seconds to e
%      data=intrpol8(data,1,300)
%
%     interpolate at 5 sps from 900 to 950 seconds using linear interp
%      data_pdiff=intrpol8(data,5,900,950,'linear')
%
%    See also: syncsr, deci, stretch, iirfilter

% check number of arguments
error(nargchk(2,5,nargin))

% check data structure
error(seischk(data,'x'))

% get timing info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
[b,e,npts,delta]=gh(data,'b','e','npts','delta');

% defaults
if(nargin<5 || isempty(method)); method{1}='pchip'; end
if(nargin<4 || isempty(new_e)); new_e=e; end
if(nargin<3 || isempty(new_b)); new_b=b; end

% number of records
nrecs=length(data);

% check and expand inputs
if(isscalar(new_b)); new_b(1:nrecs,1)=new_b;
elseif(~isvector(new_b) || length(new_b)~=nrecs)
    error('SAClab:intrpol8:badInput','dimensions of new_b are bad')
end
if(isscalar(new_e)); new_e(1:nrecs,1)=new_e;
elseif(~isvector(new_e) || length(new_e)~=nrecs)
    error('SAClab:intrpol8:badInput','dimensions of new_e are bad') 
end
if(isscalar(sr)); sr(1:nrecs,1)=sr;
elseif(~isvector(sr) || length(sr)~=nrecs)
    error('SAClab:intrpol8:badInput','dimensions of sr are bad') 
end
if(ischar(method)); method=cellstr(method); end
if(~iscellstr(method))
    error('SAClab:intrpol8:badInput','method must be char/cellstr array')
end
if(isscalar(method)); method=method(ones(nrecs,1),1);
elseif(~isvector(method) || length(method)~=nrecs)
    error('SAClab:intrpol8:badInput','dimensions of method are bad')
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
        'npts',length(nt),'depmin',-norm(min(data(i).x)),'leven','t',...
        'depmax',norm(max(data(i).x)),'depmen',norm(mean(data(i).x)));
end

end
