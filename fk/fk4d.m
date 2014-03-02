function [s4d]=fk4d(data,width,overlap,varargin)
%FK4D    Returns beamformer volumes in frequency-wavenumber-time space
%
%    Usage:    s4d=fk4d(data,width,overlap,smax,spts,frng)
%              s4d=fk4d(data,width,overlap,smax,spts,frng,polar)
%              s4d=fk4d(data,width,overlap,smax,spts,frng,polar,method)
%
%    Description:
%     S4D=FK4D(DATA,WIDTH,OVERLAP,SMAX,SPTS,FRNG) beamforms the wave energy
%     passing through an array in frequency-slowness-time space.  The
%     addition of the time dimension (compared to FKVOLUME) is made
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
%          .beam     - frequency-slowness beamforming volume
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
%          .freq     - frequency values
%          .polar    - true if slowness is sampled in polar coordinates
%          .npairs   - number of pairs
%          .method   - beamforming method (center, coarray, full, user)
%          .center   - array center as [lat lon]
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
%     S4D=FK4D(DATA,WIDTH,OVERLAP,SMAX,SPTS,FRNG,POLAR,METHOD) defines the
%     beamforming method.  METHOD may be 'center', 'coarray', 'full', or
%     [LAT LON].  The default is 'coarray' which utilizes information from
%     all unique record pairings in the beamforming and is the default.
%     The 'full' method will utilize all possible pairings including
%     pairing records with themselves and pairing records as (1st, 2nd) &
%     (2nd, 1st) making this method quite redundant and slow.  The 'center'
%     option only pairs each record against the array center (found using
%     ARRAYCENTER) and is extremely fast for large arrays compared to the
%     'coarray' & 'full' methods.  Both 'center' and 'full' methods give
%     slightly degraded results compared to 'coarray'.  Using [LAT LON] for
%     method is essentially the same as the 'center' method but uses the
%     defined coordinates as the center for the array.
%
%    Notes:
%     - Records do NOT have to start at the same time (they are windowed
%       here).  They do need to have a common sample rate and their
%       reference time and station location must be set in the header.
%
%    Examples:
%     % Get frequency-slowness-time data for an array at 20-50s periods:
%     s4d=fk4d(data,1,75,50,201,[1/50 1/20]);
%
%    See also: FKFREQSLIDE, FKFRAMESLIDE, FKVOLUME, FKMAP, FKARF

%     Version History:
%        May   9, 2010 - initial version
%        May  10, 2010 - first working version
%        May  26, 2010 - minor doc touch
%        June 16, 2010 - fix nargchk, error for no windows, improved note
%                        section, removed meaningless line
%        July  6, 2010 - major update to struct, doc update
%        Apr.  3, 2012 - minor doc update, use seizmocheck
%        Mar.  1, 2014 - improved verbose output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 14:05 GMT

% todo:

% check nargin
error(nargchk(6,8,nargin));

% check struct
error(seizmocheck(data,'dep'));

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
    
    % rethrow error
    error(lasterror);
end

% do fk analysis
try
    % verbosity
    verbose=seizmoverbose(false);
    
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
    
    % detail message
    if(verbose)
        disp('Getting Sliding Window fk Response');
        print_time_left(0,nframes);
    end
    
    % loop over time windows
    skip=0;
    for i=1:nframes
        % get data window
        start=b+(i-1)*tstep;
        finish=start+width;
        data0=interpolate(data,1/delta,[],start,finish,nan);
        
        % throw out incomplete records (those with nans)
        data0(getvaluefun(data0,@(x)any(isnan(x(:)))))=[];
        
        % skip if <2 records
        if(numel(data0)<2)
            % detail message
            if(verbose); print_time_left(i,nframes); end
            skip=skip+1;
            continue;
        end
        
        % pass to fkvolume
        % - each column corresponds to a time range
        % - each row corresponds to a frequency range
        if((i-skip)==1)
            s4d=fkvolume(data0,varargin{:});
        else
            s4d(:,i-skip)=fkvolume(data0,varargin{:});
        end
        
        % detail message
        if(verbose); print_time_left(i,nframes); end
    end
    
    % throw error if no output
    if((i-skip)==0)
        error('seizmo:fk4d:badData',...
            'No windows contained 2 or more records!');
    end
    
    % toggle verbosity back
    seizmoverbose(verbose);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle verbosity back
    seizmoverbose(verbose);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end
