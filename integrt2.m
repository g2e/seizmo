function [data]=integrt2(data)
%INTEGRT2    Integrates SAClab data records using discrete additions
%
%    Description: Calculates and returns the integral of each record using
%     a cumulative summation.  Assumes the given record is the discrete
%     derivative at the midpoints of the to-be-found integrated record.  
%     Basically it undoes the discrete differences operation of SAClab's 
%     dif function. 
%
%    Notes:
%     - Time shifts to midpoints.
%     - Increases npts by 1 (first point is 0 and is one sample interval
%       before the first midpoint).
%     - integrt2(dif(data))~=data
%     - rmean(integrt2(dif(data)))==rmean(data)
%     - inexact on unevenly sampled records
%
%    Warning:
%      If an unevenly spaced record is unable to satisfy the midpoint
%      assumption, a warning is issued and the record is integrated
%      by calling integrt.  See integrt for more info.
%
%    Usage: [data]=integrt2(data)
%
%    Examples:
%
%    See also: dif, integrt

% check nargin
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'x'))

% retreive header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
[b,e,npts,delta]=gh(data,'b','e','npts','delta');

% integrate and update header
for i=1:length(data)
    % save class and convert to double precision
    oclass=str2func(class(data(i).x));
    data(i).x=double(data(i).x);
    
    % number of components
    ncmp=size(data(i).x,2);
    
    % evenly spaced
    if(strcmp(leven(i),'true'))
        data(i).x=[zeros(1,ncmp); delta(i)*cumsum(data(i).x)];
        b(i)=b(i)-delta(i)/2; e(i)=e(i)+delta(i)/2; npts(i)=npts(i)+1;
    % unevenly spaced
    else
        % Assume each time point represents a midpoint between
        % the original time points - we want original times.
        % 
        % 1 more unknown than knowns (n midpoints, n+1 orig points)
        %
        % Solve for available solution range under condition that all time
        % intervals are positive.  Midpoint of range gives a good solution.
        %
        % Range is just smallest time interval.  This is the lack of 
        % knowledge in this case.
        data(i).t=double(data(i).t);
        dt=diff(data(i).t);
        [range,pos]=min(dt);
        
        % quick check that the midpoint assumption is plausible
        if(range<=0 || all((dt(2:end-1)-dt(1:end-2))>dt(3:end)) || ...
                all(dt(1:end-2)<(dt(2:end-1)-dt(3:end))));
            % nope - trapezoidal rule
            warning('SAClab:integrt2:failedAssumption',...
                ['Midpoint assumption failed for record %d. \n'...
                'Using the trapezoidal rule instead.'],i);
            data(i)=integrt(data(i));
            data(i).x=oclass(data(i).x);
            data(i).t=oclass(data(i).t);
            continue;
        end
        
        % find original times
        orig=zeros(npts(i)+1,1);
        orig(pos+1)=data(i).t(pos)+range/2;
        for j=pos+1:npts(i)
            orig(j+1)=2*data(i).t(j)-orig(j);
        end
        for j=pos:-1:1
            orig(j)=2*data(i).t(j)-orig(j+1);
        end
        
        % check monotonicity of times
        dto=diff(orig);
        if(min(dto)<=0)
            % non-monotonic - bail out, use trapezoidal rule
            warning('SAClab:integrt2:failedAssumption',...
                ['Midpoint assumption failed for record %d. \n'...
                'Using the trapezoidal rule instead.'],i);
            data(i)=integrt(data(i));
            data(i).x=oclass(data(i).x);
            data(i).t=oclass(data(i).t);
            continue;
        end
        
        % cumsum
        data(i).t=oclass(orig);
        data(i).x=[zeros(1,ncmp); cumsum(dto(:,ones(1,ncmp)).*data(i).x)];
        npts(i)=npts(i)+1; b(i)=data(i).t(1); e(i)=data(i).t(end);
    end
    
    % change class back
    data(i).x=oclass(data(i).x);
    
    % update header
    data(i)=ch(data(i),'depmen',norm(mean(data(i).x)),...
        'depmin',-norm(min(data(i).x)),'depmax',norm(max(data(i).x)),...
        'b',b(i),'e',e(i),'npts',npts(i));
end

end
