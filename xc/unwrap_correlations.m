function [data]=unwrap_correlations(data)
%UNWRAP_CORRELATIONS    Unwraps correlations converted to the time domain
%
%    Usage:    data=unwrap_correlations(data)
%
%    Description:
%     DATA=UNWRAP_CORRELATIONS(DATA) unwraps correlations that were
%     converted from the frequency-domain to the time-domain using IDFT.
%     This is necessary to get the lag times correct.  This also trims off
%     the "zero-padding" around the potentially non-zero values.
%
%    Notes:
%     - The trimming operation does not look at the values it is removing.
%       For coherency the values are typically non-zero (coherency is
%       spectral deconvolution) making this a loss-of-data situation.  For
%       normal cross correlation the values are barely non-zero.
%       Regardless, this operation is not reversible so converting between
%       frequency-domain & time-domain involves data loss when using this
%       function.
%
%    Header changes: E, NPTS
%
%    Examples:
%     % Correlation of records entirely in the frequency-domain:
%     data=dft(data,[],1i*2^nextpow2(2*max(getheader(data,'npts'))-1));
%     xc=correlate(data,'mcxc','fdout');
%
%     % now convert those to the time domain:
%     xc=dephase_correlations(xc,false);     % phase adjust for conversion
%     xc=idft(xc);                           % fd to td
%     xc=multiply(xc,getheader(xc,'delta')); % account for scaling
%     xc=unwrap_correlations(xc);            % unwrap to proper lag times
%
%    See also: WRAP_CORRELATIONS, DEPHASE_CORRELATIONS, CORRELATE, IDFT,
%              DFT, REVERSE_CORRELATIONS, ROTATE_CORRELATIONS,
%              SYMMETRIC_CORRELATIONS

%     Version History:
%        June 16, 2014 - initial version
%        June 23, 2014 - added checks, minor doc update
%        June 24, 2014 - example update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 24, 2014 at 11:15 GMT

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
    disp('Unwrapping Time-Domain Correlations');
    print_time_left(0,nrecs);
end

% attempt unwrap
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'NONTIME_IFTYPE','ERROR');
    
    % necessary header info
    [nspts,nx,ny,b,delta,kuser0,kuser1]=getheader(data,...
        'nspts','nxsize','nysize','b','delta','kuser0','kuser1');
    
    % correlogram signature
    if(~all(strcmp('MASTER',kuser0) & strcmp('SLAVE',kuser1)))
        error('seizmo:unwrap_correlations:notCorrelogram',...
            'Correlograms appear to have malformed headers!');
    end
    
    % operate on each record separately
    for i=1:numel(data)
        % extract & rearrange nonzero points
        data(i).dep=data(i).dep([nspts(i)-ny(i)+2:nspts(i) 1:nx(i)],1);
    end
    
    % update headers
    data=changeheader(data,'npts',nx+ny-1,'e',b+delta.*(nx+ny-2));
    
    % detail message
    if(verbose); print_time_left(nrecs,nrecs); end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % toggle verbosity back
    seizmoverbose(verbose);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % toggle verbosity back
    seizmoverbose(verbose);
    
    % rethrow error
    error(lasterror);
end

end
