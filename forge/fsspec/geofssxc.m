function [s]=geofssxc(data,ll,slow,frng,whiten,w)
%GEOFSSXC    Estimate frequency-slowness-position spectrum
%
%    Usage:    s=geofssxc(xcdata,latlon,slow,frng)
%              s=geofssxc(xcdata,latlon,slow,frng,whiten)
%              s=geofssxc(xcdata,latlon,slow,frng,whiten,weights)
%
%    Description:
%     S=GEOFSSXC(XCDATA,LATLON,SLOW,FRNG) computes an estimate of the
%     frequency-slowness-position power spectra for an array by frequency
%     domain beamforming the time series data XCDATA.  The dataset XCDATA
%     is a SEIZMO struct containing array info and precomputed cross
%     correlograms.  This function differs from FSSXC in that the waves are
%     assumed to be surface waves expanding and contracting on a sphere
%     rather than plane waves traveling across a planar surface.  This is
%     essential for characterizing surface wave sources that are within or
%     near the array (a rule of thumb: a source within one array aperture
%     width has significant near-field terms).  XCDATA is expected to be
%     correlograms with header formatting as output from CORRELATE.
%     LATLON contains the latitude and longitude wavefield beamforming
%     positions and must be in units of degrees and formatted as an Nx2
%     array of [LAT LON].  SLOW is the magnitude of the horizontal slowness
%     of the waves in sec/deg and should be a vector.  FRNG gives the
%     frequency range as [FREQLOW FREQHIGH] in Hz (the individual
%     frequencies are predetermined by the fft).  The output S is a struct
%     containing relevant info and the frequency-slowness-position spectra
%     itself (with size NPOSxNSLOWxNFREQ).  The struct layout is:
%          .spectra  - frequency-slowness-position spectra
%          .nsta     - number of stations utilized in making map
%          .stla     - station latitudes
%          .stlo     - station longitudes
%          .stel     - station elevations (surface)
%          .stdp     - station depths (from surface)
%          .butc     - UTC start time of data
%          .eutc     - UTC end time of data
%          .npts     - number of time points
%          .delta    - number of seconds between each time point
%          .latlon   - latitude/longitude positions (deg)
%          .slow     - horizontal slowness (sec/deg)
%          .freq     - frequency values (Hz)
%          .npairs   - number of pairs (aka correlograms)
%          .method   - beamforming method (always 'coarray' in this case)
%          .center   - array center as [LAT LON]
%          .vector   - [tf_multifreq tf_multislow]
%          .weights  - weights used in beamforming
%
%     S=GEOFSSXC(XCDATA,LATLON,SLOW,FRNG,WHITEN) whitens the cross spectral
%     matrix elements before beamforming if WHITEN is TRUE.  The default is
%     TRUE.
%
%     S=GEOFSSXC(XCDATA,LATLON,SLOW,FRNG,WHITEN,WEIGHTS) specifies the
%     relative weights for each correlogram in XCDATA (must match size of
%     XCDATA).
%
%    Notes:
%     - The records must have the same lag range & sample
%       spacing.
%     - Best/quickest results are obtained when XCDATA is only one
%       "triangle" of the cross correlation matrix.  This corresponds to
%       the 'coarray' method.
%
%    Examples:
%     % Do you see the 26s microseism in your data?:
%     [lat,lon]=meshgrid(-10:.5:10,-10:.5:10);
%     slow=27:0.5:33;
%     frng=[1/27 1/26];
%     s=geofssxc(xcdata,[lat(:) lon(:)],slow,frng);
%     plotgeofss(geofssavg(s));
%
%    See also: FSSXC, GEOFSS, FSS, GEOFSSAVG, GEOFSSSUB, PLOTGEOFSS,
%              GEOFSSARF, GEOFSSINFO, GEOFSSCORR, GEOFSSFRAMESLIDE,
%              GEOFSSFREQSLIDE, GEOFSSSLOWSLIDE, PLOTGEOFSSARF

%     Version History:
%        June 22, 2010 - initial version
%        July  6, 2010 - major update to struct, doc update
%        July  7, 2010 - removed deg to km conversions
%        Apr.  3, 2012 - use seizmocheck
%        May  30, 2012 - pow2pad=0 by default
%        June  4, 2012 - altered from geofkxcvolume
%        June 10, 2012 - no more conj (multiply shift by -1), verified test
%                        results
%        June 12, 2012 - add whiten option
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2012 at 14:05 GMT

% todo:

% check nargin
error(nargchk(4,6,nargin));

% check struct
error(seizmocheck(data,'dep'));

% number of correlograms
ncorr=numel(data);

% defaults for optionals
if(nargin<5 || isempty(whiten)); whiten=true; end
if(nargin<6 || isempty(w)); w=ones(numel(data),1); end

% check inputs
sf=size(frng);
if(~isreal(ll) || ndims(ll)~=2 || size(ll,2)~=2)
    error('seizmo:geofssxc:badInput',...
        'LATLON must be a Nx2 real matrix of [LAT LON]!');
elseif(~isreal(slow) || ~isvector(slow) || any(slow<=0))
    error('seizmo:geofssxc:badInput',...
        'SLOW must be a positive real vector in sec/deg!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:geofssxc:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
elseif(~isscalar(whiten) || ~islogical(whiten))
    error('seizmo:geofssxc:badInput',...
        'WHITEN must be TRUE or FALSE!');
elseif(numel(w)~=ncorr || any(w(:)<0) || ~isreal(w))
    error('seizmo:geofssxc:badInput',...
        'WEIGHTS must be equal sized with XCDATA & be positive numbers!');
end
nrng=sf(1);

% column vector slownesses
slow=slow(:);
nslow=numel(slow);

% fix lat/lon
[ll(:,1),ll(:,2)]=fixlatlon(ll(:,1),ll(:,2));
nll=size(ll,1);

% convert weights to column vector
w=w(:);
w=w./sum(w);

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
        'MULTIPLE_NPTS','ERROR',...
        'MULTIPLE_B','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% estimate fs spectra
try
    % verbosity
    verbose=seizmoverbose;
    
    % get some useful info from headers
    [npts,delta,st,ev,a,f]=getheader(data,...
        'npts','delta','st','ev','a utc','f utc');
    
    % find unique station locations
    loc=unique([st; ev],'rows');
    nsta=size(loc,1);
    
    % start/end time range (handles old correlograms)
    a=cell2mat(a); a=a(all(~isnan(a),2),:);
    f=cell2mat(f); f=f(all(~isnan(f),2),:);
    if(~isempty(a))
        [ai,ai]=min(timediff(a,a(1,:),'utc'));
        a=a(ai,:);
    else
        a=[0 0 0 0 0];
    end
    if(~isempty(f))
        [fi,fi]=max(timediff(f,f(1,:),'utc'));
        f=f(fi,:);
    else
        f=[0 0 0 0 0];
    end
    
    % check nyquist
    fnyq=1/(2*delta(1));
    if(any(frng>=fnyq))
        error('seizmo:geofssxc:badFRNG',...
            ['FRNG frequencies exceeds nyquist frequency (' ...
            num2str(fnyq) ')!']);
    end
    
    % array center
    [clat,clon]=arraycenter([st(:,1); ev(:,1)],[st(:,2); ev(:,2)]);
    
    % setup output
    [s(1:nrng,1).nsta]=deal(nsta);
    [s(1:nrng,1).stla]=deal(loc(:,1));
    [s(1:nrng,1).stlo]=deal(loc(:,2));
    [s(1:nrng,1).stel]=deal(loc(:,3));
    [s(1:nrng,1).stdp]=deal(loc(:,4));
    [s(1:nrng,1).butc]=deal(a);
    [s(1:nrng,1).eutc]=deal(f);
    [s(1:nrng,1).delta]=deal(delta(1));
    [s(1:nrng,1).npts]=deal(npts(1));
    [s(1:nrng,1).vector]=deal([true nslow~=1]);
    [s(1:nrng,1).latlon]=deal(ll);
    [s(1:nrng,1).slow]=deal(slow);
    [s(1:nrng,1).freq]=deal([]);
    [s(1:nrng,1).method]=deal('coarray');
    [s(1:nrng,1).npairs]=deal(ncorr);
    [s(1:nrng,1).center]=deal([clat clon]);
    [s(1:nrng,1).weights]=deal(w);
    
    % get frequencies
    nspts=2^nextpow2(npts(1));
    f=(0:nspts/2)/(delta(1)*nspts);  % only +freq
    
    % extract data (silently)
    seizmoverbose(false);
    data=splitpad(data);
    data=records2mat(data);
    seizmoverbose(verbose);
    
    % get fft
    data=fft(data,nspts,1);
    
    % distance difference for the phasors that "steer" the array
    % dd is NLLxNCORR
    ev=ev.'; st=st.';
    distm=sphericalinv(ll(:,ones(ncorr,1)),ll(:,2*ones(ncorr,1)),...
        ev(ones(nll,1),:),ev(2*ones(nll,1),:));
    dists=sphericalinv(ll(:,ones(ncorr,1)),ll(:,2*ones(ncorr,1)),...
        st(ones(nll,1),:),st(2*ones(nll,1),:));
    dd=-2*pi*1i*(distm-dists); % includes extra terms for efficiency
    
    % loop over frequency ranges
    for a=1:nrng
        % get frequencies
        fidx=find(f>=frng(a,1) & f<=frng(a,2));
        s(a).freq=f(fidx);
        nfreq=numel(fidx);
        if(nfreq==1); s(a).vector(1)=false; end
        
        % preallocate spectra
        s(a).spectra=zeros(nll,nslow,nfreq,'single');
        
        % warning if no frequencies
        if(~nfreq)
            warning('seizmo:geofssxc:noFreqs',...
                'No frequencies within the range %g to %g Hz!',...
                frng(a,1),frng(a,2));
            continue;
        end
        
        % extract complex cross spectra array
        % cs is normalized by the auto spectra.
        %
        % cs is NCORRxNFREQ
        cs=data(fidx,:).';
        if(whiten); cs=cs./abs(cs); end
        
        % apply weighting
        cs=cs.*w(:,ones(1,nfreq));
        
        % detail message
        if(verbose)
            fprintf('Getting Freq-Slow Spectra %d for %g to %g Hz\n',...
                a,frng(a,1),frng(a,2));
            print_time_left(0,nfreq);
        end
        
        % loop over frequencies
        for b=1:nfreq
            % loop over slownesses
            for c=1:nslow
                % get response
                s(a).spectra(:,c,b)=exp(f(fidx(b))*slow(c)*dd)*cs(:,b);
            end
            
            % detail message
            if(verbose); print_time_left(b,nfreq); end
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
    error(lasterror);
end

end
