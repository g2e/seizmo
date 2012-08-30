function [data]=differentiate(data,option)
%DIFFERENTIATE    Differentiate SEIZMO records
%
%    Usage:    data=differentiate(data)
%              data=differentiate(data,option)
%
%    Description:
%     DIFFERENTIATE(DATA) returns the derivative of each record in the
%     SEIZMO structure DATA using the differences between points as an
%     approximation of the derivative at the midpoint.  Works with unevenly
%     spaced data.
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
%       - EDGES ARE NOT SAC COMPATIBLE!
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX, NPTS, B, E
%
%    Examples:
%     % These are equal:
%     removemean(data)
%     removemean(integrate(differentiate(data)))
%
%    See also: INTEGRATE, MULTIPLYOMEGA, DIVIDEOMEGA

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 23, 2008 - glgc support
%        Feb. 28, 2008 - seischk support
%        Mar.  4, 2008 - lgcchk support
%        May  12, 2008 - dep* fix
%        June 20, 2008 - minor doc update
%        June 29, 2008 - doc update, .dep & .ind rather than .x &
%                        .t, dataless support, only calls ch once, strict
%                        filetype check
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        May   8, 2009 - uses expanded idep unit set
%        Jan. 29, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        improved messaging, forced dim on some functions
%        May   6, 2010 - slimmer code for units exchange
%        Feb. 11, 2011 - mass nargchk fix
%        Mar. 13, 2012 - doc update, seizmocheck fix, better checkheader
%                        usage, use getheader improvements
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2012 at 15:05 GMT

% todo:
% - add 'five' option for uneven
% - 'three' and 'five' are not really 3pt and 5pt
%   - 'three' is just 'two' basically
%   - 'five' is just a 4pt
% - lets look for a more stencils and a better naming scheme
% - can we undo other integrate options?

% check nargin
error(nargchk(1,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt differentiation
try
    % check headers
    data=checkheader(data,...
        'NONTIME_IFTYPE','ERROR');
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);

    % retreive header info
    [delta,b,e,npts,leven,idep]=getheader(data,...
        'delta','b','e','npts','leven lgc','idep id');
    leven=~strcmpi(leven,'false');

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
    
    % detail message
    if(verbose)
        disp('Differentiating Record(s)');
        print_time_left(0,nrecs);
    end

    % take derivative and update header
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % dataless support
        if(~npts(i))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);

        % 2pt = (dep(j+1) -dep(j))/delta
        if(two(i))
            % evenly spaced
            if(leven(i))
                data(i).dep=diff(data(i).dep,1,1)/delta(i);
                b(i)=b(i)+delta(i)/2;
                e(i)=e(i)-delta(i)/2;
                npts(i)=npts(i)-1;
            % unevenly spaced
            else
                data(i).ind=double(data(i).ind);
                ind=diff(data(i).ind,1,1);
                data(i).dep=diff(data(i).dep,1,1)...
                    ./ind(:,ones(1,size(data(i).dep,2)));
                data(i).ind=oclass(data(i).ind(1:end-1)+ind/2);
                npts(i)=npts(i)-1;
                b(i)=data(i).ind(1);
                e(i)=data(i).ind(end);
            end
        % 3pt = (dep(j+1) - dep(j-1))/(2*delta)
        elseif(three(i))
            % evenly spaced
            if(leven(i))
                data(i).dep=...
                    filter([1 0 -1],1,data(i).dep,[],1)/(2*delta(i));
                data(i).dep=data(i).dep(3:end,:);
                b(i)=b(i)+delta(i);
                e(i)=e(i)-delta(i);
                npts(i)=npts(i)-2;
            % unevenly spaced
            else
                data(i).ind=double(data(i).ind);
                ind=diff(data(i).ind,1,1);
                ind=ind(:,ones(1,size(data(i).dep,2)));
                ind2=filter([1 0 -1],1,data(i).ind,[],1);
                ind2=ind2(3:end,ones(1,size(data(i).dep,2)));
                data(i).dep=(ind(1:end-1,:).*data(i).dep(3:end,:)...
                    +ind(2:end).*data(i).dep(1:end-2,:))./ind2;
                data(i).ind=oclass(data(i).ind(2:end-1));
                npts(i)=npts(i)-2;
                b(i)=data(i).ind(1);
                e(i)=data(i).ind(end);
            end
        % 5pt = (-dep(j+2)+8*dep(j+1)-8*dep(j-1)+dep(j-2))/(12*delta)
        else
            % evenly spaced
            if(leven(i))
                data(i).dep=...
                    filter([-1 8 0 -8 1],1,data(i).dep,[],1)/(12*delta(i));
                data(i).dep=data(i).dep(5:end,:);
                b(i)=b(i)+2*delta(i);
                e(i)=e(i)-2*delta(i);
                npts(i)=npts(i)-4;
            % unevenly spaced
            else
                % need to derive this using Lagrange multipliers IIRC
                error('seizmo:differentiate:undone',...
                    ['Record: %d\n5pt differentiation is not ' ...
                    'available for unevenly sampled records!']);
            end
        end

        % change class back
        data(i).dep=oclass(data(i).dep);

        % get dep* info
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % change dependent component type
    newunit={'ipop' 'icrackle' 'isnap' 'ijerk' 'iacc' 'ivel' 'idisp' ...
        'iabsmnt' 'iabsity' 'iabseler' 'iabserk' 'iabsnap' 'iabsackle'};
    [tf,idx]=ismember(idep,[newunit(2:end) {'iabspop'}]);
    idep(tf)=newunit(idx(tf));
    idep(~tf)={'iunkn'};

    % update header
    data=changeheader(data,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax,...
        'idep',idep,'b',b,'e',e,'npts',npts);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
