function [data]=make_adjoint_source(data,win)
%MAKE_ADJOINT_SOURCE    Convert records to adjoint sources
%
%    Usage:    data=make_adjoint_source(data,win)
%
%    Description:
%     DATA=MAKE_ADJOINT_SOURCE(DATA,WIN) converts synthetic records created
%     through a forward simulation to adjoint source records for creation
%     of finite frequency kernels.  This basically windows the negative
%     derivitive of the records in SEIZMO struct DATA using the WIN
%     parameter (must be size NRECSx2 where NRECS is the number of records
%     in DATA).  Records are then scaled by their autocorrelation.  Records
%     must be in displacement and must be timeseries or xy data.
%
%    Notes:
%     - Really just for SPECFEM3D.
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     % Generally this goes something along the following:
%     data=readseizmo('*');
%     gcarc=getheader(data,'gcarc');
%     data(gcarc<100 | gcarc>150)=[];
%     pdiff=findpicks(arrivals2picks(data,'Pdiff'),'Pdiff',1);
%     data=make_adjoint_source(data,pdiff+[-20 80]);
%     writeseizmo(data,'append','.adj');
%
%    See also: MAKE_SOURCE_TIMEFUNCTION, CONVOLVE_SOURCE_TIMEFUNCTION,
%              DECONVOLVE_SOURCE_TIMEFUNCTION

%     Version History:
%        Dec.  8, 2009 - initial version
%        Mar.  1, 2010 - updated for newer checking methods
%        Feb.  1, 2011 - update for triangletf changes
%        Mar. 13, 2012 - doc update, use getheader improvements,
%                        seizmoverbose support, better checkheader usage
%        Mar. 15, 2012 - fix example
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 02:15 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt rest
try
    % check header
    data=checkheader(data,...
        'FALSE_LEVEN','ERROR',...
        'NONTIME_IFTYPE','ERROR');
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);
    
    % check window
    if(~isreal(win) || ~isequal(size(win),[nrecs 2]))
        error('seizmo:make_adjoint_source:badWin',...
            'WIN must be a real-valued array of size Nx2!');
    end
    
    % get window parameters
    cwin=sum(win,2)/2;
    hwin=diff(win,1,2)/2;

    % get header info
    [b,npts,delta,idep]=getheader(data,'b','npts','delta','idep id');
    
    % check units
    if(any(~strcmpi(idep,'idisp')))
        error('seizmo:make_adjoint_source:badUnits',...
            ['Record(s):\n' sprintf('%d ',find(~strcmpi(idep,'idisp'))) ...
            '\nDependent component units must be in displacement!']);
    end
    
    % detail message
    if(verbose)
        disp('Converting Synthetic(s) for Adjoint Source Computation');
        print_time_left(0,nrecs);
    end

    % loop over records
    depmin=nan(nrecs,1); depmen=depmin; depmax=depmin;
    for i=1:numel(data)
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end
        
        % convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        
        % get window function (a triangle window)
        % - note need to set amp to 1 to better preserve amplitudes
        tw=triangletf(b(i)+(0:delta(i):delta(i)*(npts(i)-1)),...
            cwin(i),hwin(i),1);

        % get derivative
        data(i).dep=gradient(data(i).dep,delta(i));

        % get adjoint source
        % (negative velocity trace of a triangle windowed phase
        %  normalized by its similarly windowed autocorrolation)
        data(i).dep=(-data(i).dep.*tw(:))...
            /(delta(i)*sum(tw(:).*data(i).dep.^2,1));
        
        % convert back
        data(i).dep=oclass(data(i).dep);
        
        % get dep*
        depmin=min(data(i).dep(:));
        depmen=nanmean(data(i).dep(:));
        depmax=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % update header
    data=changeheader(data,...
        'depmin',depmin,'depmen',depmen,'depmax',depmax);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
