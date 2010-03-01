function [data]=make_adjoint_source(data,win)
%MAKE_ADJOINT_SOURCE    Convert records to adjoint sources
%
%    Usage:    data=make_adjoint_source(data,win)
%
%    Description: DATA=MAKE_ADJOINT_SOURCE(DATA,WIN) converts synthetic
%     records created through a forward simulation to adjoint source
%     records for creation of finite frequency kernels.  This basically
%     windows the negative derivitive of the records in SEIZMO struct DATA
%     using the WIN parameter (must be size NRECSx2 where NRECS is the
%     number of records in DATA).  Records are then scaled by their
%     autocorrelation.  Records must be in displacement and must be
%     timeseries or xy data.
%
%    Notes:
%     - Really just for SPECFEM3D.
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     Generally this goes something along the following:
%      data=readseizmo('*');
%      gcarc=getheader(data,'gcarc');
%      data(gcarc<100 | gcarc>150)=[];
%      data=addarrivals(data,'phases','Pdiff');
%      pdiff=getarrival(data,'Pdiff');
%      data=make_adjoint_source(data,pdiff+[-20 80]);
%      writeseizmo(data,'append','.adj');
%
%    See also: MAKE_SOURCE_TIMEFUNCTION, CONVOLVE_SOURCE_TIMEFUNCTION,
%              DECONVOLVE_SOURCE_TIMEFUNCTION

%     Version History:
%        Dec.  8, 2009 - initial version
%        Mar.  1, 2010 - updated for newer checking methods
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2010 at 02:15 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% attempt rest
try
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
    [b,npts,delta]=getheader(data,'b','npts','delta');
    [idep,iftype]=getenumid(data,'idep','iftype');
    leven=getlgc(data,'leven');
    
    % check sample spacing, file type, units
    if(any(~strcmpi(leven,'true')))
        error('seizmo:make_adjoint_source:badSpacing',...
            ['Record(s):\n' sprintf('%d ',find(~strcmpi(leven,'true'))) ...
            'Illegal operation on unevenly spaced record(s)!'])
    elseif(any(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy')))
        error('seizmo:make_adjoint_source:badType',...
            ['Record(s):\n' sprintf('%d ',...
            find(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy'))) ...
            'Illegal operation on spectral/xyz record(s)!'])
    elseif(any(~strcmpi(idep,'idisp')))
        error('seizmo:make_adjoint_source:badUnits',...
            ['Record(s):\n' sprintf('%d ',find(~strcmpi(idep,'idisp'))) ...
            '\nDependent component units must be in displacement!']);
    end

    % loop over records
    depmin=nan(nrecs,1); depmen=depmin; depmax=depmin;
    for i=1:numel(data)
        % skip dataless
        if(isempty(data(i).dep)); continue; end
        
        % convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        
        % get window function (a triangle window)
        tw=triangletf(b(i)+(0:npts(i)-1)*delta(i),cwin(i),hwin(i));

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
        depmen=mean(data(i).dep(:));
        depmax=max(data(i).dep(:));
    end
    
    % update header
    data=changeheader(data,...
        'depmin',depmin,'depmen',depmen,'depmax',depmax);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
