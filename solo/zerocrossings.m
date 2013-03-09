function [zx]=zerocrossings(data)
%ZEROCROSSINGS    Returns the zero-crossing times of SEIZMO records
%
%    Usage:    zx=zerocrossings(data)
%
%    Description:
%     ZX=ZEROCROSSINGS(DATA) returns zero-crossing times for each record
%     in SEIZMO struct DATA.  ZX is a Nx1 cell array of times where ZX{1}
%     gives the zero-crossings for the 1st record and so on for all N
%     records.  Zero-crossings are determined by points with opposite sign
%     and linear interpolation to a better estimate of the crossing time.
%     The times are in the relative time convention.
%
%    Notes:
%
%    Examples:
%     % The number of zero-crossings in a portion of time gives an estimate
%     % of the dominant period in that time section:
%     zx=zerocrossings(data(1));
%     n=sum(zx{1}>=0 && zx{1}<=100); % number of crossings from 0 to 100s
%     T=(2*100)/n; % 2 crossings per period
%
%    See also: GETVALUEFUN

%     Version History:
%        Jan. 18, 2011 - initial version
%        Mar. 13, 2012 - leven bugfix, better checkheader usage,
%                        seizmocheck fix, use getheader improvements
%        Feb. 14, 2013 - bugfix the leven bugfix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2013 at 12:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check seizmo struct
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt zero-crossings discovery
try
    % check headers
    data=checkheader(data);
    
    % number of records
    nrecs=numel(data);
    
    % get pertinent header info
    [b,delta,leven]=getheader(data,'b','delta','leven lgc');
    leven=~strcmpi(leven,'false');
    
    % loop over each record
    zx=cell(nrecs,1);
    for i=1:nrecs
        % index of lower point for each crossing
        idx=find(filter([1 1],1,sign(data(i).dep))==0)-1;
        
        % get times based on whether evenly spaced or not
        if(~leven(i))
            t1=data(i).ind(idx);
            t2=data(i).ind(idx+1);
        else
            t1=b(i)+delta(i)*(idx-1);
            t2=b(i)+delta(i)*idx;
        end
        
        % interpolate to best point
        zx{i}=t1+(t2-t1)...
            .*abs(data(i).dep(idx)./(data(i).dep(idx+1)-data(i).dep(idx)));
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
