function [r,t]=geofsshorz(n,e,ll,slow,frng,method,whiten,w)
%GEOFSSHORZ    Estimate frequency-slowness-position spectrum of horizontals
%
%    Usage:   [r,t]=geofsshorz(n,e,latlon,slow,frng)
%             [r,t]=geofsshorz(n,e,latlon,slow,frng,method)
%             [r,t]=geofsshorz(n,e,latlon,slow,frng,method,whiten)
%             [r,t]=geofsshorz(n,e,latlon,slow,frng,method,whiten,weights)
%
%    Description:
%     [R,T]=GEOFSSHORZ(N,E,LATLON,SLOW,FRNG) computes an estimate of the
%     radial & transverse frequency-slowness-position power spectra
%     for an array by frequency domain beamforming horizontal components in
%     the time series dataset DATA.  The dataset DATA is a SEIZMO struct
%     containing array info and time series recordings. This function
%     differs from FSSHORZ in that the waves are assumed to be surface
%     waves expanding and contracting on a sphere rather than plane waves
%     traveling across a planar surface.  This is essential for
%     characterizing surface wave sources that are within or near the array
%     (a rule of thumb: a source within one array aperture width has
%     significant near-field terms).  LATLON contains the latitude and
%     longitude wavefield beamforming positions and must be in units of
%     degrees and formatted as an Nx2 array of [LAT LON].  SLOW is the
%     magnitude of the horizontal slowness of the waves in sec/deg and
%     should be a vector.  FRNG gives the frequency range as
%     [FREQLOW FREQHIGH] in Hz (the individual frequencies are
%     predetermined by the fft).  The outputs R & T are structs containing
%     relevant info and the frequency-slowness-position spectra for the
%     radial and transverse wavefields (with size NPOSxNSLOWxNFREQ).  The
%     struct layout is:
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
%          .method   - beamforming method ('user', 'center' or 'coarray')
%          .center   - array center as [LAT LON]
%          .vector   - [tf_multifreq tf_multislow]
%          .weights  - weights used in beamforming
%
%     [R,T]=GEOFSSHORZ(N,E,LATLON,SLOW,FRNG,METHOD) sets the beamforming
%     method.  METHOD may be 'center', 'coarray', 'full' or [LAT LON].  The
%     default is 'center' which is extremely fast for large arrays compared
%     to the 'coarray' method as it only pairs each record against the
%     array center (found using ARRAYCENTER) to compute the real-valued
%     spectrum.  The 'coarray' method utilizes information from all unique
%     pairings of records to compute the complex slowness spectrum while
%     the 'full' method uses every possible pairing to do the same.  The
%     'full' method is significantly slower and gives degraded results
%     compared to the 'coarray' method and so is not recommended except in
%     verification.  The 'center' method gives results that are the same as
%     the 'full' method but does it far faster.  Using [LAT LON] for method
%     is algorithmically the same as the 'center' method but uses the
%     defined coordinates as the center for the array.
%
%     [R,T]=GEOFSSHORZ(N,E,LATLON,SLOW,FRNG,METHOD,WHITEN) whitens the
%     cross spectral matrix elements before beamforming if WHITEN is TRUE.
%     The default is TRUE.
%
%     [R,T]=GEOFSSHORZ(N,E,LATLON,SLOW,FRNG,METHOD,WHITEN,WEIGHTS) sets
%     the relative weights for each record in DATA (must match the size of
%     DATA).
%
%    Notes:
%     - Records must have equal and regular sample spacing.
%
%    Examples:
%     % Horizontals need to be rotated to North/East prior to calling
%     % GEOFSSHORZ.  Here an example of how to do that:
%     data=rotate(data,'to',0,'kcmpnm1','N','kcmpnm2','E');
%     [lat,lon]=meshgrid(-10:.5:10,-10:.5:10);
%     slow=27:0.5:33;
%     frng=[1/27 1/26];
%     [r,t]=geofsshorz(data(1:2:end),data(2:2:end),...
%         [lat(:) lon(:)],slow,frng);
%     plotgeofss(geofssavg(r));
%     plotgeofss(geofssavg(t));
%
%    See also: FSSHORZ, GEOFSSHORZXC, FSSHORZXC, GEOFSSAVG, GEOFSSSUB,
%              PLOTGEOFSS, GEOFSSARF, GEOFSSDBINFO, GEOFSSCORR,
%              GEOFSSFRAMESLIDE, GEOFSSFREQSLIDE, GEOFSSSLOWSLIDE,
%              PLOTGEOFSSARF

%     Version History:
%        Aug. 30, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 30, 2012 at 15:15 GMT

% todo:

% check nargin/nargout
error(nargchk(5,8,nargin));

% check struct
error(seizmocheck(n,'dep'));
error(seizmocheck(e,'dep'));

% number of records
nrecs=numel(n);
if(nrecs~=numel(e))
    error('seizmo:geofsshorz:badInput',...
        'N & E datasets are not the same size!');
end

% need 2+ records
if(nrecs<2)
    error('seizmo:geofsshorz:arrayTooSmall',...
        'N & E datasets must have 2+ records!');
end

% defaults for optionals
if(nargin<5 || isempty(method)); method='center'; end
if(nargin<6 || isempty(whiten)); whiten=true; end
if(nargin<7 || isempty(w)); w=ones(numel(n),1); end

% valid method strings
valid.METHOD={'center' 'coarray' 'full' 'capon'};

% check inputs
sf=size(frng);
if(~isreal(ll) || ndims(ll)~=2 || size(ll,2)~=2)
    error('seizmo:geofsshorz:badInput',...
        'LATLON must be a Nx2 real matrix of [LAT LON]!');
elseif(~isreal(slow) || ~isvector(slow) || any(slow<=0))
    error('seizmo:geofsshorz:badInput',...
        'SLOW must be a positive real vector in sec/deg!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:geofsshorz:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
elseif((isnumeric(method) && (~isreal(method) || ~numel(method)==2)) ...
        || (ischar(method) && ~any(strcmpi(method,valid.METHOD))))
    error('seizmo:geofsshorz:badInput',...
        ['METHOD must be one of the following:\n' ...
        '''CENTER'', ''COARRAY'', ''FULL'', ''CAPON'' or [LAT LON]!']);
elseif(~isscalar(whiten) || ~islogical(whiten))
    error('seizmo:geofsshorz:badInput',...
        'WHITEN must be TRUE or FALSE!');
elseif(~isnumeric(w) || numel(w)~=nrecs || any(w(:)<0))
    error('seizmo:geofsshorz:badInput',...
        'WEIGHTS must be equal sized with DATA & be positive!');
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
    n=checkheader(n,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'NONINTEGER_REFTIME','ERROR',...
        'UNSET_REFTIME','ERROR',...
        'OUTOFRANGE_REFTIME','ERROR',...
        'UNSET_ST_LATLON','ERROR');
    e=checkheader(e,...
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

% estimate fs spectra
try
    % verbosity
    verbose=seizmoverbose;
    
    % grab necessary header info
    [butc,eutc,b,delta,npts,st,kname,cmp]=getheader(n,...
        'b utc','e utc','b','delta','npts','st','kname','cmp');
    [butc2,b2,delta2,npts2,st2,kname2,cmp2]=getheader(e,...
        'b utc','b','delta','npts','st','kname','cmp');
    butc=cell2mat(butc); eutc=cell2mat(eutc);
    
    % need n & e to match on several fields:
    % kname, st, b utc, delta, npts
    if(~isequal(butc,butc2) || ~isequal(b,b2) || ~isequal(delta,delta2) ...
            || ~isequal(npts,npts2) || ~isequal(st,st2) ...
            || ~isequal(kname(:,1:3),kname2(:,1:3)))
        error('seizmo:geofsshorz:badInput',...
            'N & E are not synced!');
    end
    clear butc2 b b2 delta2 npts2 st2 kname kname2;
    delta=delta(1);
    
    % require north & east orientation
    if(any(cmp(:,1)~=90) || any(cmp2(:,1)~=90)) % cmpinc
        error('seizmo:geofsshorz:badInput',...
            'N & E must have CMPINC==90!');
    elseif(any(cmp(:,2)~=0))
        error('seizmo:geofsshorz:badInput',...
            'N must have CMPAZ==0!');
    elseif(any(cmp2(:,2)~=90))
        error('seizmo:geofsshorz:badInput',...
            'E must have CMPAZ==90!');
    end
    clear cmp cmp2;
    
    % get utc static time shifts
    dt=timediff(butc,butc(1,:),'utc').'; % need row vector
    switch method
        case 'coarray'
            [row,col]=find(triu(true(nrecs),1));
            w=w(col).*w(row);     % expand weights
            w=w./sum(w);          % renormalize weights
            dt=dt(row)-dt(col);   % alter to pair shifts
        case {'full' 'capon'}
            [row,col]=find(true(nrecs));
            w=w(col).*w(row);     % expand weights
            w=w./sum(w);          % renormalize weights
            dt=dt(row)-dt(col);   % alter to pair shifts
    end
    dt=2*pi*1i*dt(ones(nll,1),:);
    
    % longest record 
    [maxnpts,imaxnpts]=max(npts);
    
    % check nyquist
    fnyq=1/(2*delta(1));
    if(any(frng>=fnyq))
        error('seizmo:geofsshorz:badFRNG',...
            ['FRNG frequencies exceeds nyquist frequency (' ...
            num2str(fnyq) ')!']);
    end
    
    % fix method/center/npairs
    if(ischar(method))
        method=lower(method);
        [clat,clon]=arraycenter(st(:,1),st(:,2));
        switch method
            case 'coarray'
                npairs=nrecs*(nrecs-1)/2;
            case {'full' 'capon'}
                npairs=nrecs*nrecs;
            case 'center'
                npairs=nrecs;
        end
    else
        clat=method(1);
        clon=method(2);
        method='user';
        npairs=nrecs;
    end
    
    % setup output
    [r(1:nrng,1).nsta]=deal(nrecs);
    [r(1:nrng,1).stla]=deal(st(:,1));
    [r(1:nrng,1).stlo]=deal(st(:,2));
    [r(1:nrng,1).stel]=deal(st(:,3));
    [r(1:nrng,1).stdp]=deal(st(:,4));
    [r(1:nrng,1).butc]=deal(butc(1,:));
    [r(1:nrng,1).eutc]=deal(eutc(imaxnpts,:));
    [r(1:nrng,1).delta]=deal(delta(1));
    [r(1:nrng,1).npts]=deal(maxnpts);
    [r(1:nrng,1).vector]=deal([true nslow~=1]);
    [r(1:nrng,1).latlon]=deal(ll);
    [r(1:nrng,1).slow]=deal(slow);
    [r(1:nrng,1).freq]=deal([]);
    [r(1:nrng,1).method]=deal(method);
    [r(1:nrng,1).npairs]=deal(npairs);
    [r(1:nrng,1).center]=deal([clat clon]);
    [r(1:nrng,1).weights]=deal(w);
    t=r;
    
    % get frequencies
    nspts=2^nextpow2(maxnpts); % half xcorr
    %nspts=2^nextpow2(2*maxnpts-1); % full xcorr for verification
    f=(0:nspts/2)/(delta(1)*nspts); % only +freq
    
    % extract data (silently)
    seizmoverbose(false);
    n=records2mat(n);
    e=records2mat(e);
    seizmoverbose(verbose);
    
    % get fft
    n=fft(n,nspts,1);
    e=fft(e,nspts,1);
    
    % distance difference for the phasors that "steer" the array
    % dd is NLLxNPAIRS
    % - also getting azimuth info for rotators
    switch method
        case {'coarray' 'full' 'capon'}
            % distance difference for each pair from source
            ev=st(col,:).'; st=st(row,:).';
            [distm,azm,bazm]=sphericalinv(...
                ll(:,ones(npairs,1)),ll(:,2*ones(npairs,1)),...
                ev(ones(nll,1),:),ev(2*ones(nll,1),:));
            [dists,azs,bazs]=sphericalinv(...
                ll(:,ones(npairs,1)),ll(:,2*ones(npairs,1)),...
                st(ones(nll,1),:),st(2*ones(nll,1),:));
            dd=2*pi*1i*(distm-dists); % includes extra terms for efficiency
        otherwise
            % distance difference for station and center from source
            st=st.';
            [distm,azm,bazm]=sphericalinv(...
                ll(:,ones(npairs,1)),ll(:,2*ones(npairs,1)),...
                clat(ones(nll,1),ones(npairs,1)),...
                clon(ones(nll,1),ones(npairs,1)));
            [dists,azs,bazs]=sphericalinv(...
                ll(:,ones(npairs,1)),ll(:,2*ones(npairs,1)),...
                st(ones(nll,1),:),st(2*ones(nll,1),:));
            dd=2*pi*1i*(distm-dists); % includes extra terms for efficiency
    end
    
    die
    % - i need to get the rotator formulas worked out for n & e inputs
    % - need az/baz between pairs (for r & t)
    
    % rotators
    % u = cos theta
    % v = sin theta
    az=az.'; az=az(ones(nll,1),:);
    baz=baz.'; baz=baz(ones(nll,1),:);
    thetam=d2r*(bazm-az+180);
    thetas=d2r*(bazs-baz);
    um=cos(thetam);
    vm=sin(thetam);
    us=cos(thetas);
    vs=sin(thetas);
    
    % loop over frequency ranges
    for a=1:nrng
        % get frequencies
        fidx=find(f>=frng(a,1) & f<=frng(a,2));
        rvol(a).freq=f(fidx);
        tvol(a).freq=f(fidx);
        nfreq=numel(fidx);
        
        % preallocate fk space
        rvol(a).beam=zeros(nll,nslow,nfreq,'single');
        tvol(a).beam=zeros(nll,nslow,nfreq,'single');
        
        % warning if no frequencies
        if(~nfreq)
            warning('seizmo:geofsshorz:noFreqs',...
                'No frequencies within the range %g to %g Hz!',...
                frng(a,1),frng(a,2));
            continue;
        end
        
        % detail message
        if(verbose)
            fprintf('Getting geofk Volume %d for %g to %g Hz\n',...
                a,frng(a,1),frng(a,2));
            print_time_left(0,nfreq);
        end
        
        % loop over frequencies
        for b=1:nfreq
            % current freq idx
            cf=fidx(b);
            
            % loop over slownesses
            for c=1:nslow
                % get beam
                
                % rotate data into direction of plane wave for all pairs
                data=um.*us.*rr(cf*ones(1,nll),:) ...
                    -um.*vs.*rt(cf*ones(1,nll),:) ...
                    -vm.*us.*tr(cf*ones(1,nll),:) ...
                    +vm.*vs.*tt(cf*ones(1,nll),:);
                
                % normalize by auto spectra
                data=data./abs(data);
                
                % apply weights
                data=data.*w(ones(1,nll),:);
                
                % now getting fk beam for radial (rayleigh)
                rvol(a).beam(:,c,b)=10*log10(abs(real(...
                    sum(data.*exp(f(cf)*s(c)*dd),2))));
                
                % rotate data ortho to plane wave direction for all pairs
                data=vm.*vs.*rr(cf*ones(1,nll),:) ...
                    +vm.*us.*rt(cf*ones(1,nll),:) ...
                    +um.*vs.*tr(cf*ones(1,nll),:) ...
                    +um.*us.*tt(cf*ones(1,nll),:);
                
                % normalize by auto spectra
                data=data./abs(data);
                
                % apply weights
                data=data.*w(ones(1,nll),:);
                
                % now getting fk beam for tangential (love)
                tvol(a).beam(:,c,b)=10*log10(abs(real(...
                    sum(data.*exp(f(cf)*s(c)*dd),2))));
            end
            
            % detail message
            if(verbose); print_time_left(b,nfreq); end
        end
        
        % normalize so max peak is at 0dB
        rvol(a).normdb=max(rvol(a).beam(:));
        rvol(a).beam=rvol(a).beam-rvol(a).normdb;
        tvol(a).normdb=max(tvol(a).beam(:));
        tvol(a).beam=tvol(a).beam-tvol(a).normdb;
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
