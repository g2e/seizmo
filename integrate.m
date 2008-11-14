function [data]=integreat(data,option)
%INTEGREAT    Integrates SAClab data records
%
%    Description:
%
%    Notes:
%     - why does trapezoid not lose a point here while in sac it does?
%       - sac is dumb and drops the first point (zero) which does provide
%         info so should not be dropped
%     - wtf is sac doing for rectangular?
%       - they dont qc their shit
%
%    Tested on: Matlab r2007b
%
%    Usage:    data=integreat(data)
%              data=integreat(data,'rectangular'|'trapezoidal'|'midpoint')
%
%    Examples:
%
%    See also: dif, divomega, mulomega

%     Version History:
%        Nov. 12, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 12, 2008 at 07:10 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin))

% check data structure
error(seischk(data,'dep'))

% turn off struct checking
oldstate=get_seischk_state;
set_seischk_state(true);

% check data header
data=chkhdr(data);

% number of records
nrecs=numel(data);

% check option
if(nargin==1 || isempty(option))
    % default
    option='trapezoidal';
elseif(~ischar(option) && ~iscellstr(option))
    error('SAClab:integreat:badOption',...
        'OPTION must be a character or cell string array!');
end

% expand scalar option
option=cellstr(option);
nopt=numel(option);
if(nopt==1)
    option(1:nrecs,1)=option;
elseif(nopt~=nrecs)
    error('SAClab:integreat:badOption',...
        'OPTION must be a single option or one option per record!');
end

% get header info
leven=glgc(data,'leven');
idep=genum(data,'idep');
[b,e,delta,npts]=gh(data,'b','e','delta','npts');

% integrate
[depmen,depmin,depmax]=swap(nan(nrecs,1));
for i=1:nrecs
    % skip dataless
    if(~data(i).hasdata); continue; end
    if(~npts(i)); continue; end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    switch option{i}
        case 'rectangular'
            if (strcmpi(leven(i),'true'))
                data(i).dep=delta(i)*cumsum(data(i).dep);
            else
                
            end
        case 'rectangular-sac'
            if (strcmpi(leven(i),'true'))
                data(i).dep=delta(i)*cumsum(data(i).dep);
            else
                
            end
        case 'trapezoidal'
            if (strcmpi(leven(i),'true'))
                data(i).dep=delta(i)*cumtrapz(data(i).dep);
            else
                data(i).dep=cumtrapz(double(data(i).ind),data(i).dep);
            end
        case 'trapezoidal-sac'
            if (strcmpi(leven(i),'true'))
                data(i).dep=delta(i)*cumtrapz(data(i).dep);
                data(i).dep(1,:)=[];
                b(i)=b(i)+delta(i)/2; 
                e(i)=e(i)-delta(i)/2; 
                npts(i)=npts(i)-1;
            else
                data(i).dep=cumtrapz(double(data(i).ind),data(i).dep);
                data(i).dep(1,:)=[];
                b(i)=b(i)+(data(i).ind(2)-data(i).ind(1))/2;
                e(i)=e(i)-(data(i).ind(end)-data(i).ind(end-1))/2;
                npts(i)=npts(i)-1;
            end
        case 'midpoint'
            
        otherwise
            error('SAClab:integreat:badOption','Unknown OPTION!');
    end
    
    % change class back
    data(i).dep=oclass(data(i).dep);
    
    % get dependent component values
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
end

% change dependent component type
isacc=strcmpi(idep,'iacc');
isvel=strcmpi(idep,'ivel');
idep(isacc)={'ivel'};
idep(isvel)={'idisp'};
idep(~isacc & ~isvel)={'iunkn'};

% update headers
data=ch(data,'depmen',depmen,'depmin',depmin,'depmax',depmax,...
    'b',b,'e',e,'npts',npts,'idep',idep);

% toggle struct checking back
set_seischk_state(oldstate);

end
