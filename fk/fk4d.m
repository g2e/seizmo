function [s4d]=fk4d(data,width,overlap,varargin)
%FK4D    Returns a map of energy in frequency-wavenumber-time space
%
%    Usage:    s4d=fk4d(data,width,overlap,smax,spts,frng)
%              s4d=fk4d(data,width,overlap,smax,spts,frng,polar)
%              s4d=fk4d(data,width,overlap,smax,spts,frng,polar,center)
%
%    Description: S4D=FK4D(DATA,WIDTH,OVERLAP,SMAX,SPTS,FRNG) calculates
%     the energy passing through an array in frequency-wavenumber-time
%     space.  Actually, to allow for easier interpretation between
%     frequencies, the energy is mapped into frequency-slowness-time space.
%     The addition of the time dimension (compared to FKVOLUME) is made
%     possible by taking windows that are a fraction of the input records
%     in SEIZMO struct DATA.  The array info is derived from the metadata
%     stored in the SEIZMO struct DATA, so make sure station location and
%     timing fields are set!  WIDTH is the width of the time windows in
%     percent of the total time extent of the records & OVERLAP is the
%     percent of the overlap between each window.  The default WIDTH is
%     1% and the default OVERLAP is 0%.  The range of the slowness space is
%     specified by SMAX (in s/deg) and extends from -SMAX to SMAX for both
%     East/West and North/South directions.  SPTS controls the number of
%     slowness points for both directions (SPTSxSPTS grid).  FRNG gives the
%     frequency range as [FREQLOW FREQHIGH] in Hz.  S4D is a struct array
%     whose elements are essentially time frames of the frequency-slowness
%     space and contain relevant info as well as the frequency-slowness
%     data.  The struct layout is:
%          .response - frequency-slowness array response
%          .nsta     - number of stations utilized in making map
%          .stla     - station latitudes
%          .stlo     - station longitudes
%          .stel     - station elevations (surface)
%          .stdp     - station depths (from surface)
%          .butc     - UTC start time of data
%          .eutc     - UTC end time of data
%          .npts     - number of time points
%          .delta    - number of seconds between each time point
%          .x        - east/west slowness or azimuth values
%          .y        - north/south or radial slowness values
%          .z        - frequency values
%          .polar    - true if slowness is sampled in polar coordinates 
%          .center   - array center or method
%          .normdb   - what 0dB actually corresponds to
%          .volume   - true if frequency-slowness volume (false for FKMAP)
%
%     S4D=FK4D(DATA,WIDTH,OVERLAP,SMAX,SPTS,FRNG,POLAR) specifies if the
%     slowness space is sampled regularyly in cartesian or polar
%     coordinates.  Polar coords are useful for slicing the volume by
%     azimuth (pie slice) or slowness magnitude (rings).  Cartesian coords
%     (the default) samples the slowness space regularly in the East/West &
%     North/South directions and so exhibits less distortion of the
%     slowness space.
%
%     S4D=FK4D(DATA,WIDTH,OVERLAP,SMAX,SPTS,FRNG,POLAR,CENTER) defines the
%     array center.  CENTER may be [LAT LON], 'center', 'coarray', or
%     'full'.  The default is 'coarray'.  The 'center' option finds the
%     center position of the array by averaging the station positions
%     (using ARRAYCENTER).  Both 'coarray' and 'full' are essentially
%     centerless methods using the relative positioning between every
%     possible pairing of stations in the array.  The 'full' method
%     includes redundant and same station pairings (and will always give
%     poorer results compared to 'coarray').
%
%    Notes:
%     - Records in DATA must have equal number of points, equal sample
%       spacing, the same start time (in absolute time), and be evenly
%       spaced time series records.  Use functions SYNCHRONIZE, SYNCRATES,
%       & INTERPOLATE to get the timing/sampling the same.
%
%    Examples:
%     Get frequency-slowness-time data for an array at 20-50s periods:
%      s4d=fk4d(data,1,75,50,201,[1/50 1/20]);
%
%    See also: FKFREQSLIDE, FKTIMESLIDE, FKVOLUME, FKMAP, FKARF

%     Version History:
%        May   9, 2010 - initial version
%        May  10, 2010 - first working version
%        May  26, 2010 - minor doc touch
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  26, 2010 at 10:40 GMT

% todo:

% check nargin
msg=nargchk(6,8,nargin);
if(~isempty(msg)); error(msg); end

% check struct
versioninfo(data,'dep');

% defaults for width/overlap
if(isempty(width)); width=1; end
if(isempty(overlap)); overlap=0; end

% check width/overlap
if(~isreal(width) || ~isscalar(width) || width<0 || width>100)
    error('seizmo:fk4d:badInput',...
        'WIDTH must be a percent value between 0 & 100!');
elseif(~isreal(overlap) || ~isscalar(overlap) || overlap<0 || overlap>100)
    error('seizmo:fk4d:badInput',...
        'OVERLAP must be a percent value between 0 & 100!');
end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'NONINTEGER_REFTIME','ERROR',...
        'UNSET_REFTIME','ERROR',...
        'OUTOFRANGE_REFTIME','ERROR',...
        'UNSET_ST_LATLON','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
end

% do fk analysis
try
    % verbosity
    verbose=seizmoverbose;
    
    % synchronize to start of first record
    data=synchronize(data,'b','first');
    
    % get time range
    [b,e,delta]=getheader(data,'b','e','delta');
    b=min(b);
    e=max(e);
    delta=delta(1); % required all the same
    
    % get time frame setup
    width=width/100;
    overlap=overlap/100;
    nframes=fix((1-width*overlap)/(width-width*overlap));
    width=width*(e-b);
    tstep=width*(1-overlap);
    
    % loop over time windows
    skip=0;
    for i=1:nframes
        % detail message
        if(verbose)
            fprintf('Getting fk Response Window %d of %d\n',i,nframes);
        end
        
        % get data window
        start=b+(i-1)*tstep;
        finish=start+width;
        data0=interpolate(data,1/delta,[],start,finish,nan);
        
        % throw out incomplete records (those with nans)
        data0(getvaluefun(data0,@(x)any(isnan(x(:)))))=[];
        
        % skip if <2 records
        if(numel(data0)<2); skip=skip+1; continue; end
        
        % pass to fkvolume
        % - each column corresponds to a time range
        % - each row corresponds to a frequency range
        if((i-skip)==1)
            s4d=fkvolume(data0,varargin{:});
            s4d=s4d(:);
        else
            s4d(:,i-skip)=fkvolume(data0,varargin{:});
        end
    end
    
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
