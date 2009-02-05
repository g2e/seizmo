function [data]=differentiate(data,option)
%DIFFERENTIATE    Differentiate SEIZMO records
%
%    Description: DIFFERENTIATE(DATA) returns the derivative of each
%     record in the SEIZMO structure DATA using the differences between
%     points as an approximation of the derivative at the midpoint.  Works
%     with unevenly spaced data.
%
%     DIFFERENTIATE(DATA,OPTION) allows specifying the difference operator.
%     Supported OPTIONS are 'two', 'three', and 'five'. The operators
%     are as follows:
%      'two'   = (dep(j+1)-dep(j))/delta(i)
%      'three' = (dep(j+1)-dep(j-1))/(2*delta(i))
%      'five'  = (-dep(j+2)+8*dep(j+1)-8*dep(j-1)+dep(j-2))/(12*delta(i))
%
%     The default option is 'two'.  Option 'five' throws an error for
%      unevenly spaced records.
%    
%    Notes:
%     - for option 'two'
%       - timing is shifted to midpoints
%       - B increased by DELTA/2, E decreased by DELTA/2
%       - NPTS decreased by 1
%     - for option 'three'
%       - B increased by DELTA, E decreased by DELTA
%       - NPTS decreased by 2
%     - for option 'five'
%       - B increased by 2*DELTA, E decreased by 2*DELTA
%       - NPTS decreased by 4
%       - NOT SAC COMPATIBLE!
%
%    Tested on: Matlab r2007b
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, NPTS, B, E
%
%    Usage:    data=differentiate(data)
%
%    Examples:
%     These are equal:
%      removemean(data)
%      removemean(integrate(differentiate(data)))
%
%    See also: integrate, multiplyomega, divideomega

%    Version History:
%       Jan. 28, 2008 - initial version
%       Feb. 23, 2008 - glgc support
%       Feb. 28, 2008 - seischk support
%       Mar.  4, 2008 - lgcchk support
%       May  12, 2008 - dep* fix
%       June 20, 2008 - minor doc update
%       June 29, 2008 - doc update, .dep & .ind rather than .x &
%                       .t, dataless support, only calls ch once, strict
%                       filetype check
%
%    Written by Garrett Euler (ggeuler at wustl dot edu)
%    Last Updated June 29, 2008 at 07:25 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin))

% check data structure
error(seizmocheck(data,'dep'))

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% check for unsupported filetypes
iftype=getenumdesc(data,'iftype');
if(strcmpi(iftype,'General XYZ (3-D) file'))
    error('seizmo:differentiate:illegalFiletype',...
        'Illegal operation on xyz file!');
elseif(any(strcmpi(iftype,'Spectral File-Real/Imag') | ...
        strcmpi(iftype,'Spectral File-Ampl/Phase')))
    error('seizmo:differentiate:illegalFiletype',...
        'Illegal operation on spectral file!');
end

% retreive header info
even=strcmpi('true',getlgc(data,'leven'));
[delta,b,e,npts]=getheader(data,'delta','b','e','npts');
idep=getenumid(data,'idep');

% number of records
nrecs=numel(data);

% check option
if(nargin==1 || isempty(option))
    option='two';
elseif(~ischar(option) || ~iscellstr(option))
    error('seizmo:differentiate:badOption',...
        'OPTION must be char or cellstr!');
end
two=strcmpi(option,'two');
three=strcmpi(option,'three');
five=strcmpi(option,'five');
if(any(~(two | three | five)))
    error('seizmo:differentiate:badOption',...
        'OPTION must be ''two'' ''three'' or ''five''!');
end
if(~any(numel(two)==[1 nrecs]))
    error('seizmo:differentiate:badOption',...
        'OPTION must be one option per record or one option for all!');
end
two(1:nrecs,1)=two;
three(1:nrecs,1)=three;

% take derivative and update header
depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
for i=1:nrecs
    % dataless support
    if(~npts(i)); continue; end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % 2pt = (dep(j+1) -dep(j))/delta
    if(two(i))
        % evenly spaced
        if(even(i))
            data(i).dep=diff(data(i).dep)/delta(i);
            b(i)=b(i)+delta(i)/2; e(i)=e(i)-delta(i)/2; npts(i)=npts(i)-1;
        % unevenly spaced
        else
            data(i).ind=double(data(i).ind);
            ind=diff(data(i).ind);
            data(i).dep=diff(data(i).dep)...
                ./ind(:,ones(1,size(data(i).dep,2)));
            data(i).ind=oclass(data(i).ind(1:end-1)+ind/2);
            npts(i)=npts(i)-1; b(i)=data(i).ind(1); e(i)=data(i).ind(end);
        end
    % 3pt = (dep(j+1) - dep(j-1))/(2*delta)
    elseif(three(i))
        % evenly spaced
        if(even(i))
            data(i).dep=filter([1 0 -1],1,data(i).dep)/(2*delta(i));
            data(i).dep=data(i).dep(3:end,:);
            b(i)=b(i)+delta(i); e(i)=e(i)-delta(i); npts(i)=npts(i)-2;
        % unevenly spaced
        else
            data(i).ind=double(data(i).ind);
            ind=diff(data(i).ind);
            ind=ind(:,ones(1,size(data(i).dep,2)));
            ind2=filter([1 0 -1],1,data(i).ind);
            ind2=ind2(3:end,ones(1,size(data(i).dep,2)));
            data(i).dep=(ind(1:end-1,:).*data(i).dep(3:end,:)...
                +ind(2:end).*data(i).dep(1:end-2,:))./ind2;
            data(i).ind=oclass(data(i).ind(2:end-1));
            npts(i)=npts(i)-2; b(i)=data(i).ind(1); e(i)=data(i).ind(end);
        end
    % 5pt = (-dep(j+2)+8*dep(j+1)-8*dep(j-1)+dep(j-2))/(12*delta)
    else
        % evenly spaced
        if(even(i))
            data(i).dep=filter([-1 8 0 -8 1],1,data(i).dep)/(12*delta(i));
            data(i).dep=data(i).dep(5:end,:);
            b(i)=b(i)+2*delta(i); e(i)=e(i)-2*delta(i); npts(i)=npts(i)-4;
        % unevenly spaced
        else
            
            error('seizmo:differentiate:undone','5pt not for uneven!');
        end
    end
    
    % change class back
    data(i).dep=oclass(data(i).dep);
    
    % get dep* info
    depmen(i)=mean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
end

% change dependent component type
isdis=strcmpi(idep,'idisp');
isvel=strcmpi(idep,'ivel');
idep(isdis)={'ivel'};
idep(isvel)={'iacc'};
idep(~isdis & ~isvel)={'iunkn'};

% update header
data=changeheader(data,'depmen',depmen,'depmin',depmin,'depmax',depmax,...
    'idep',idep,'b',b,'e',e,'npts',npts);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
