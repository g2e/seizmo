function [data]=symmetric_correlations(data)
%SYMMETRIC_CORRELATIONS    Combines causal & acausal correlation components
%
%    Usage:    data=symmetric_correlations(data)
%
%    Description:
%     DATA=SYMMETRIC_CORRELATIONS(DATA) stacks the acausal and causal
%     components of time-domain noise correlation functions returning the
%     "symmetric" component.
%
%    Notes:
%     - Reference:
%        Bensen et al 2007, GJI, doi:10.1111/j.1365-246X.2007.03374.x
%
%    Header changes: B, E, NPTS, DEP*
%
%    Examples:
%     % Read in correlations and get symmetric component:
%     ncfs=readseizmo('my/ncf/dir/*');
%     sym_ncfs=symmetric_correlations(ncfs);
%
%    See also: CORRELATE, NOISE_PROCESS, NOISE_STACK, NOISE_C3,
%              ROTATE_CORRELATIONS, UNWRAP_CORRELATIONS, WRAP_CORRELATIONS,
%              REVERSE_CORRELATIONS, DEPHASE_CORRELATIONS

%     Version History:
%        June 16, 2014 - initial version
%        June 23, 2014 - added checks & docs, bugfix: forgot to average
%        June 25, 2014 - use addrecords simple replacement for speed
%        July 21, 2014 - fast path for symmetrically sampled records (5x),
%                        improved speed for asymmetric too (~2x) by using
%                        solofun instead
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 21, 2014 at 11:15 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% turn off verbosity
verbose=seizmoverbose(false);

% detail message
if(verbose); disp('EXTRACTING SYMMETRIC COMPONENT OF CORRELATIONS'); end

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
    [delta,b,e,kuser0,kuser1]=getheader(data,...
        'delta','b','e','kuser0','kuser1');
    
    % correlogram signature
    if(~all(strcmp('MASTER',kuser0) & strcmp('SLAVE',kuser1)))
        error('seizmo:symmetric_correlations:notCorrelogram',...
            'Correlograms appear to have malformed headers!');
    end
    
    % maximum lag time of symmetric component
    maxl=delta.*floor(max(abs(b),abs(e))./delta);
    
    % change class to double
    [data,oclass]=changeclass(data,'double');
    
    % who needs interpolation
    slow=b~=-e | mod(b,delta);
    
    % interpolate
    if(any(slow))
        data(slow)=interpolate(data(slow),1./delta,[],-maxl,maxl,0);
    end
    
    % fold & stack
    data=solofun(data,@(x)(x((end+1)/2:-1:1)+x((end+1)/2:end))/2);
    
    % fix header
    data=changeheader(data,'b',0,'e',maxl);
    
    % convert class back
    data=changeclass(data,oclass);
    
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
