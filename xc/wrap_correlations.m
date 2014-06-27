function [data]=wrap_correlations(data)
%WRAP_CORRELATIONS    Wraps correlations to convert to the frequency domain
%
%    Usage:    data=wrap_correlations(data)
%
%    Description:
%     DATA=WRAP_CORRELATIONS(DATA) wraps correlation data in preparation
%     for conversion to the frequency-domain from the time-domain.  This
%     involves interpolation to restore the discrete time sampling exactly,
%     zero-padding to get to the proper number of points, and wrapping the
%     zero-lag value to the beginning (which is not actually the zero-lag
%     value if the correlated records did not start at the same time - an
%     unfortunate fact that requires correction of the phase in the
%     frequency-domain using DEPHASE_CORRELATIONS).
%
%    Notes:
%
%    Header changes: B, E, NPTS, DEP*
%
%    Examples:
%     % Correlate records to get time-domain correlations:
%     xc=correlate(data,'mcxc');
%
%     % now convert those to the frequency domain:
%     xc=wrap_correlations(xc);            % wrap & pad for conversion
%     xc=divide(xc,getheader(xc,'delta')); % account for scaling factor
%     xc=dft(xc,'rlim');                   % td to fd
%     xc=dephase_correlations(xc);         % phase adjust
%
%    See also: UNWRAP_CORRELATIONS, DEPHASE_CORRELATIONS, CORRELATE, IDFT,
%              DFT, REVERSE_CORRELATIONS, ROTATE_CORRELATIONS,
%              SYMMETRIC_CORRELATIONS

%     Version History:
%        June 16, 2014 - initial version
%        June 23, 2014 - added checks & docs
%        June 24, 2014 - example update
%        June 25, 2014 - reltime now uses t3 field
%        June 26, 2014 - speed edit: interpolate only when necessary
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 26, 2014 at 11:15 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% turn off verbosity
verbose=seizmoverbose(false);

% number of records
nrecs=numel(data);

% detail message
if(verbose)
    disp('Wrapping Time-Domain Correlations');
    print_time_left(0,nrecs);
end

% attempt header check
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'NONTIME_IFTYPE','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % toggle verbosity back
    seizmoverbose(verbose);
    
    % rethrow error
    error(lasterror);
end

% attempt wrap
try
    % necessary header info
    [nspts,nx,ny,delta,b1,b2,b01,b02,isabs,kuser0,kuser1,npts,b]=...
        getheader(data,'nspts','nxsize','nysize','delta','a utc',...
        't0 utc','a','t3','t2','kuser0','kuser1','npts','b');
    
    % correlogram signature
    if(~all(strcmp('MASTER',kuser0) & strcmp('SLAVE',kuser1)))
        error('seizmo:wrap_correlations:notCorrelogram',...
            'Correlograms appear to have malformed headers!');
    end
    
    % lag range that corresponds to the original data
    b0=-delta.*(ny-1);
    if(isabs); b0=timediff(cell2mat(b1),cell2mat(b2),'utc')+b0;
    else b0=b02-b01+b0;
    end
    npts0=nx+ny-1;
    
    % who does not have original data
    no=b~=b0 | npts~=npts0;
    
    % interpolate those to original data
    if(any(no))
        data(no)=interpolate(data(no),1./delta(no),[],b0(no),...
            b0(no)+delta(no).*(npts0(no)-1),0);
    end
    
    % pad & wrap
    for i=1:numel(data)
        data(i).dep=[data(i).dep(end-nx+1:end);
            zeros(nspts(i)-npts0(i),1);
            data(i).dep(1:ny-1)];
    end
    
    % update headers
    data=changeheader(data,'npts',nspts,'e',b0+delta.*(nspts-1));
    
    % detail message
    if(verbose); print_time_left(nrecs,nrecs); end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % toggle verbosity back
    seizmoverbose(verbose);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % toggle verbosity back
    seizmoverbose(verbose);
    
    % rethrow error
    error(lasterror);
end

end
