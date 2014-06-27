function [data]=dephase_correlations(data,dflag)
%DEPHASE_CORRELATIONS    Fixes the phase of frequency-domain correlations
%
%    Usage:    data=dephase_correlations(data)
%              data=dephase_correlations(data,false)
%
%    Description:
%     DATA=DEPHASE_CORRELATIONS(DATA) removes the phase introduced by the
%     difference in start times for the underlying correlated records.
%     This is only a potential issue for frequency-domain correlations.
%
%     DATA=DEPHASE_CORRELATIONS(DATA,FALSE) applies the phase introduced by
%     the difference in start times for the underlying correlated records.
%     This undoes a dephase operation without this flag.
%
%    Notes:
%     - The function CORRELATE calls this function for frequency-domain
%       output.  So you probably will not use this function unless you are
%       converting between time-domain and frequency-domain correlations
%       and want to retain every scrap of data.
%
%    Header changes: SB, DEP*
%
%    Examples:
%     % Correlate records, returning frequency-domain
%     % results and convert to the time-domain:
%     xc=correlate(data,'mcxc','normxc','fdout');
%     xc=dephase_correlations(xc,false);    % phase adjust for conversion
%     xc=idft(xc);                          % fd to td
%     xc=multiply(xc,getheader(xc,'delta'); % account for scaling
%     xc=unwrap_correlations(xc);           % unwrap to proper lag times
%
%    See also: CORRELATE, UNWRAP_CORRELATIONS, WRAP_CORRELATIONS,
%              OMEGASHIFT, REVERSE_CORRELATIONS, DEPHASE_CORRELATIONS,
%              SYMMETRIC_CORRELATIONS, ROTATE_CORRELATIONS

%     Version History:
%        June 23, 2014 - initial version
%        June 24, 2014 - direction option
%        June 25, 2014 - doc update, reltime now uses t3 field
%        June 26, 2014 - adjust sb field so conversion to the td no longer
%                        requires undoing this operation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 26, 2014 at 11:15 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'NONSPECTRAL_IFTYPE','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% attempt wrap
try
    % default/check dephase direction flag
    if(nargin<2 || isempty(dflag)); dflag=true; end
    if(~isscalar(dflag) || ~islogical(dflag))
        error('seizmo:dephase_correlations:badInput',...
            'DFLAG must be TRUE or FALSE!');
    end
    
    % convert logical to sign
    dflag=-1+2*dflag;
    
    % necessary header info
    [b1,b2,b01,b02,isabs,kuser0,kuser1,sb]=getheader(data,...
        'a utc','t0 utc','a','t3','t2','kuser0','kuser1','sb');
    
    % correlogram signature
    if(~all(strcmp('MASTER',kuser0) & strcmp('SLAVE',kuser1)))
        error('seizmo:dephase_correlations:notCorrelogram',...
            'Correlograms appear to have malformed headers!');
    end
    
    % time shift that corresponds to the original data
    if(isabs); b0=timediff(cell2mat(b1),cell2mat(b2),'utc');
    else b0=b02-b01;
    end
    
    % dephase
    data=omegashift(data,dflag*b0);
    
    % adjust sb field
    data=changeheader(data,'sb',sb+dflag*b0);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end
