function [s]=geofss(data,ll,slow,frng,method,whiten,w)
%GEOFSS    Estimate frequency-slowness-position spectrum
%
%    Usage:    s=geofss(data,latlon,slow,frng)
%              s=geofss(data,latlon,slow,frng,method)
%              s=geofss(data,latlon,slow,frng,method,whiten)
%              s=geofss(data,latlon,slow,frng,method,whiten,weights)
%
%    Description:
%     S=GEOFSS(DATA,LATLON,SLOW,FRNG) computes an estimate of the
%     frequency-slowness-position power spectra for an array by frequency
%     domain beamforming the time series dataset DATA.  The dataset DATA
%     is a SEIZMO struct containing array info and time series recordings.
%     This function differs from FSS in that the waves are assumed to be
%     surface waves expanding and contracting on a sphere rather than plane
%     waves traveling across a planar surface.  This is essential for
%     characterizing surface wave sources that are within or near the array
%     (a rule of thumb: a source within one array aperture width has
%     significant near-field terms).  LATLON contains the latitude and
%     longitude wavefield beamforming positions and must be in units of
%     degrees and formatted as an Nx2 array of [LAT LON].  SLOW is the
%     magnitude of the horizontal slowness of the waves in sec/deg and
%     should be a vector.  FRNG gives the frequency range as
%     [FREQLOW FREQHIGH] in Hz (the individual frequencies are
%     predetermined by the fft).  The output S is a struct containing
%     relevant info and the frequency-slowness-position spectra
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
%          .method   - beamforming method ('user', 'center' or 'coarray')
%          .center   - array center as [LAT LON]
%          .vector   - [tf_multifreq tf_multislow]
%          .weights  - weights used in beamforming
%
%     S=GEOFSS(DATA,LATLON,SLOW,FRNG,METHOD) defines the beamforming
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
%     S=GEOFSS(DATA,LATLON,SLOW,FRNG,METHOD,WHITEN) whitens the cross
%     spectral matrix elements before beamforming if WHITEN is TRUE.  The
%     default is TRUE.
%
%     S=GEOFSS(DATA,LATLON,SLOW,FRNG,METHOD,WHITEN,WEIGHTS) specifies the
%     relative weights for each record in DATA (must match size of DATA).
%
%    Notes:
%     - Records in DATA must have equal and regular sample spacing.
%
%    Examples:
%     % Do you see the 26s microseism in your data?:
%     [lat,lon]=meshgrid(-10:.5:10,-10:.5:10);
%     slow=27:0.5:33;
%     frng=[1/27 1/26];
%     s=geofss(data,[lat(:) lon(:)],slow,frng);
%     plotgeofss(geofssavg(s));
%
%    See also: FSS, GEOFSSXC, FSSXC, GEOFSSAVG, GEOFSSSUB, PLOTGEOFSS,
%              GEOFSSARF, GEOFSSDBINFO, GEOFSSCORR, GEOFSSFRAMESLIDE,
%              GEOFSSFREQSLIDE, GEOFSSSLOWSLIDE, PLOTGEOFSSARF

%     Version History:
%        June 22, 2010 - initial version
%        July  6, 2010 - major update to struct, doc update
%        July  7, 2010 - removed deg to km conversions
%        Apr.  3, 2012 - use seizmocheck
%        May  30, 2012 - pow2pad=0 by default
%        June  4, 2012 - altered from geofssxc & fkvolume
%        June 10, 2012 - fixed weighting, allow b to vary, added full
%                        method, verified test results
%        June 11, 2012 - default to center method
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2012 at 14:05 GMT

% todo:

% check nargin
error(nargchk(4,7,nargin));

% check struct
error(seizmocheck(data,'dep'));

% number of records
nrecs=numel(data);

% need 2+ records
if(nrecs<2)
    error('seizmo:geofss:arrayTooSmall',...
        'DATA must have 2+ records!');
end

% defaults for optionals
if(nargin<5 || isempty(method)); method='center'; end
if(nargin<6 || isempty(whiten)); whiten=true; end
if(nargin<7 || isempty(w)); w=ones(numel(data),1); end

% valid method strings
valid.METHOD={'center' 'coarray' 'full'};

% check inputs
sf=size(frng);
if(~isreal(ll) || ndims(ll)~=2 || size(ll,2)~=2)
    error('seizmo:geofss:badInput',...
        'LATLON must be a Nx2 real matrix of [LAT LON]!');
elseif(~isreal(slow) || ~isvector(slow) || any(slow<=0))
    error('seizmo:geofss:badInput',...
        'SLOW must be a positive real vector in sec/deg!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:geofss:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
elseif((isnumeric(method) && (~isreal(method) || ~numel(method)==2)) ...
        || (ischar(method) && ~any(strcmpi(method,valid.METHOD))))
    error('seizmo:geofss:badInput',...
        'METHOD must be ''CENTER'', ''COARRAY'', or [LAT LON]!');
elseif(~isscalar(whiten) || ~islogical(whiten))
    error('seizmo:geofss:badInput',...
        'WHITEN must be TRUE or FALSE!');
elseif(~isnumeric(w) || numel(w)~=nrecs || any(w(:)<0))
    error('seizmo:geofss:badInput',...
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

% estimate fs spectra
try
    % verbosity
    verbose=seizmoverbose;
    
    % require all records to have equal delta, b utc, and be single cmp
    % - we could drop some requirements but that means more effort
    %   - this is useful for surface waves & large aperture arrays
    %   - requires static time shift term
    [npts,delta,butc,eutc,st]=getheader(data,...
        'npts','delta','b utc','e utc','st');
    butc=cell2mat(butc); eutc=cell2mat(eutc);
    
    % get utc static time shifts
    dt=timediff(butc,butc(1,:),'utc').'; % need row vector
    switch method
        case 'coarray'
            [row,col]=find(triu(true(nrecs),1));
            w=w(col).*w(row);     % expand weights
            w=w./sum(w);          % renormalize weights
            dt=dt(row)-dt(col);   % alter to pair shifts
        case 'full'
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
        error('seizmo:geofss:badFRNG',...
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
            case 'full'
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
    [s(1:nrng,1).nsta]=deal(nrecs);
    [s(1:nrng,1).stla]=deal(st(:,1));
    [s(1:nrng,1).stlo]=deal(st(:,2));
    [s(1:nrng,1).stel]=deal(st(:,3));
    [s(1:nrng,1).stdp]=deal(st(:,4));
    [s(1:nrng,1).butc]=deal(butc(1,:));
    [s(1:nrng,1).eutc]=deal(eutc(imaxnpts,:));
    [s(1:nrng,1).delta]=deal(delta(1));
    [s(1:nrng,1).npts]=deal(maxnpts);
    [s(1:nrng,1).vector]=deal([true nslow~=1]);
    [s(1:nrng,1).latlon]=deal(ll);
    [s(1:nrng,1).slow]=deal(slow);
    [s(1:nrng,1).freq]=deal([]);
    [s(1:nrng,1).method]=deal(method);
    [s(1:nrng,1).npairs]=deal(npairs);
    [s(1:nrng,1).center]=deal([clat clon]);
    [s(1:nrng,1).weights]=deal(w);
    
    % get frequencies
    nspts=2^nextpow2(maxnpts); % half xcorr
    %nspts=2^nextpow2(2*maxnpts-1); % full xcorr for verification
    f=(0:nspts/2)/(delta(1)*nspts); % only +freq
    
    % extract data (silently)
    seizmoverbose(false);
    data=records2mat(data);
    seizmoverbose(verbose);
    
    % get fft
    data=fft(data,nspts,1);
    
    % distance difference for the phasors that "steer" the array
    % dd is NLLxNPAIRS
    switch method
        case {'coarray' 'full'}
            % distance difference for each pair from source
            ev=st(col,:).'; st=st(row,:).';
            distm=sphericalinv(...
                ll(:,ones(npairs,1)),ll(:,2*ones(npairs,1)),...
                ev(ones(nll,1),:),ev(2*ones(nll,1),:));
            dists=sphericalinv(...
                ll(:,ones(npairs,1)),ll(:,2*ones(npairs,1)),...
                st(ones(nll,1),:),st(2*ones(nll,1),:));
            dd=2*pi*1i*(distm-dists); % includes extra terms for efficiency
        otherwise
            % distance difference for station and center from source
            st=st.';
            distm=sphericalinv(...
                ll(:,ones(npairs,1)),ll(:,2*ones(npairs,1)),...
                clat(ones(nll,1),ones(npairs,1)),...
                clon(ones(nll,1),ones(npairs,1)));
            dists=sphericalinv(...
                ll(:,ones(npairs,1)),ll(:,2*ones(npairs,1)),...
                st(ones(nll,1),:),st(2*ones(nll,1),:));
            dd=2*pi*1i*(distm-dists); % includes extra terms for efficiency
    end
    
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
            warning('seizmo:geofss:noFreqs',...
                'No frequencies within the range %g to %g Hz!',...
                frng(a,1),frng(a,2));
            continue;
        end
        
        % build complex cross spectra array, cs
        % cs is normalized by the auto spectra.
        %
        % note that center/user do not require
        % an array but do need normalization
        %
        % cs is NPAIRSxNFREQ
        cs=data(fidx,:).';
        if(whiten); cs=cs./abs(cs); end
        switch method
            case {'coarray' 'full'}
                cs=cs(row,:).*conj(cs(col,:));
        end
        
        % apply weighting
        cs=cs.*w(:,ones(1,nfreq));
        
        % detail message
        if(verbose)
            fprintf('Getting Freq-Slow Spectra %d for %g to %g Hz\n',...
                a,frng(a,1),frng(a,2));
            print_time_left(0,nfreq);
        end
        
        % beamforming method
        switch method
            case {'coarray' 'full'}
                for b=1:nfreq
                    for c=1:nslow
                        s(a).spectra(:,c,b)=...
                            exp(f(fidx(b))*(slow(c)*dd-dt))*cs(:,b);
                    end
                    if(verbose); print_time_left(b,nfreq); end
                end
            otherwise
                for b=1:nfreq
                    for c=1:nslow
                        s(a).spectra(:,c,b)=abs(...
                            exp(f(fidx(b))*(slow(c)*dd-dt))*cs(:,b)).^2;
                    end
                    if(verbose); print_time_left(b,nfreq); end
                end
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
