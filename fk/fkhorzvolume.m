function [varargout]=fkhorzvolume(edat,ndat,smax,spts,frng,polar,method,w)
%FKHORZVOLUME    Returns frequency-wavenumber space for horz. data
%
%    Usage:    [rvol,tvol]=fkhorzvolume(e,n,smax,spts,frng)
%              [rvol,tvol]=fkhorzvolume(e,n,smax,spts,frng,polar)
%              [rvol,tvol]=fkhorzvolume(e,n,smax,spts,frng,polar,method)
%              [rvol,tvol]=fkhorzvolume(...
%                   e,n,smax,spts,frng,polar,method,weights)
%
%    Description:
%     [RVOL,TVOL]=FKHORZVOLUME(E,N,SMAX,SPTS,FRNG) calculates the Rayleigh
%     & Love energy moving through an array in frequency-wavenumber space
%     utilizing the horizontal datasets E & N.  To allow for easier
%     interpretation between frequencies, the energy is mapped into
%     frequency-slowness space.  The array info is contained in the SEIZMO
%     structs E & N.  It is expected that E contains the East-oriented
%     horizontal data and N has the North-oriented data.  This differs from
%     FKVOLUME in that horizontals are utilized to retreive both the
%     radial (RVOL) & transverse (TVOL) energy distributions.  FKVOLUME
%     does not account for the directional sensitivity of horizontals - it
%     is better suited for vertical components that do not vary in
%     directional sensitivity to plane waves with different propagation
%     directions.  The range of the slowness space is given by SMAX (in
%     s/deg) and extends from -SMAX to SMAX for both East/West and
%     North/South directions.  SPTS controls the number of slowness points
%     for both directions (SPTSxSPTS grid).  FRNG gives the frequency range
%     as [FREQLOW FREQHIGH] in Hz.  RVOL & TVOL are structs containing
%     relevant info and the frequency-slowness volume itself.  The struct
%     layout is:
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
%          .z        - frequency values
%          .polar    - true if slowness is sampled in polar coordinates 
%          .center   - array center or method
%          .normdb   - what 0dB actually corresponds to
%          .volume   - true if frequency-slowness volume (false for FKMAP)
%          .weights  - weights used in beamforming
%
%     Calling FKXCHORZVOLUME with no outputs will automatically slide
%     through the frequency-slowness volumes using FKFREQSLIDE.
%
%     [RVOL,TVOL]=FKXCHORZVOLUME(E,N,SMAX,SPTS,FRNG,POLAR) sets if the
%     slowness space is sampled regularly in cartesian or polar
%     coordinates.  Polar coords are useful for slicing the volume by
%     azimuth (pie slice) or slowness magnitude (rings).  Cartesian coords
%     (the default) samples the slowness space regularly in the East/West
%     & North/South directions and so exhibits less distortion of the
%     slowness space.
%
%     [RVOL,TVOL]=FKXCHORZVOLUME(E,N,SMAX,SPTS,FRNG,POLAR,METHOD) defines
%     the beamform method.  METHOD may be 'center', 'coarray', 'full', or
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
%     [RVOL,TVOL]=FKXCHORZVOLUME(E,N,SMAX,SPTS,FRNG,POLAR,METHOD,WEIGHTS)
%     specifies weights for each seismogram in E/N (the size of E, N, &
%     WEIGHTS must be the same) for use in beamforming.  The weights are
%     normalized internally to sum to 1.
%
%    Notes:
%     - Records in E & N Tmust have the same sample spacing, number of
%       points, start time and stop time.  Furthermore, the datasets must
%       be equal sized and every record should correspond to the record
%       with the same index in the other dataset (ie E(3) & N(3) must
%       correspond to the same station pair).
%
%    Examples:
%     One way to perform horizontal fk analysis:
%      [rvol,tvol]=fkhorzvolume(e,n,50,101,[1/50 1/20]);
%      fkfreqslide(rvol,0);
%      fkfreqslide(tvol,0);
%
%    See also: FKXCVOLUME, FKFREQSLIDE, FKMAP, FK4D, FKVOL2MAP, FKSUBVOL,
%              FKVOLUME, CORRELATE, ROTATE_CORRELATIONS

%     Version History:
%        Nov. 18, 2010 - initial version
%        May  30, 2012 - pow2pad=0 by default
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  30, 2012 at 15:25 GMT

% todo:

% check nargin
error(nargchk(5,8,nargin));

% define some constants
d2r=pi/180;
d2km=6371*d2r;

% check struct
versioninfo(edat,'dep');
versioninfo(ndat,'dep');

% make sure rr/tt are the same size
nrecs=numel(edat);
nrecs1=numel(ndat);
if(~isequal(nrecs,nrecs1))
    error('seizmo:fkhorzvolume:unmatchedXCdata',...
        'Datasets do not match in size!');
end

% need 2+ records
if(nrecs<2)
    error('seizmo:fkhorzvolume:arrayTooSmall',...
        'DATA must have 2+ records!');
end

% defaults for optionals
if(nargin<6 || isempty(polar)); polar=false; end
if(nargin<7 || isempty(method)); method='coarray'; end
if(nargin<8); w=[]; end

% check inputs
sf=size(frng);
if(~isreal(smax) || ~isscalar(smax) || smax<=0)
    error('seizmo:fkhorzvolume:badInput',...
        'SMAX must be a positive real scalar in s/deg!');
elseif(~any(numel(spts)==[1 2]) || any(fix(spts)~=spts) || any(spts<=2))
    error('seizmo:fkhorzvolume:badInput',...
        'SPTS must be a positive scalar integer >2!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:fkhorzvolume:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
elseif(~isscalar(polar) || (~islogical(polar) && ~isnumeric(polar)))
    error('seizmo:fkhorzvolume:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif((isnumeric(method) && (~isreal(method) || ~numel(method)==2)) ...
        || (ischar(method) && ~any(strcmpi(method,valid.METHOD))))
    error('seizmo:fkhorzvolume:badInput',...
        'METHOD must be ''CENTER'', ''COARRAY'', ''FULL'', or [LAT LON]!');
elseif(~isempty(w) && (any(w(:)<0) || ~isreal(w) || sum(w(:))==0))
    error('seizmofkhorzvolume:badInput',...
        'WEIGHTS must be positive real numbers!');
end
nrng=sf(1);

% convert weights to row vector
w=w(:).';
w=w./sum(w);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    edat=checkheader(edat,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'MULTIPLE_NPTS','ERROR',...
        'NONINTEGER_REFTIME','ERROR',...
        'UNSET_REFTIME','ERROR',...
        'OUTOFRANGE_REFTIME','ERROR',...
        'UNSET_ST_LATLON','ERROR');
    ndat=checkheader(ndat,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'MULTIPLE_NPTS','ERROR',...
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
    verbose=seizmoverbose;
    
    % require east & north components
    [eaz,einc]=getheader(edat,'cmpaz','cmpinc');
    [naz,ninc]=getheader(ndat,'cmpaz','cmpinc');
    if(any(eaz~=90) || any(einc~=90))
        error('seizmo:fkhorzvolume:badE',...
            'E does not appear to be East horizontal data!');
    elseif(any(naz~=90) || any(ninc~=90))
        error('seizmo:fkhorzvolume:badN',...
            'N does not appear to be North horizontal data!');
    end
    
    % get station locations & require that they match
    est=getheader(edat,'st');
    nst=getheader(ndat,'st');
    stla=est(:,1); stlo=est(:,2);
    if(~isequal(est,nst))
        error('seizmo:fkhorzvolume:DataMismatch',...
            'Datasets have different station locations!');
    end
    
    % require all records to have equal b, npts, delta
    [butc,eutc,npts,delta]=getheader([edat(:); ndat(:)],...
        'b utc','e utc','npts','delta');
    butc=cell2mat(butc); eutc=cell2mat(eutc);
    if(size(unique(butc,'rows'),1)~=1)
        error('seizmo:fkhorzvolume:badData',...
            'Records in DATA must have equal B (UTC)!');
    end
    if(~isscalar(unique(npts)) || ~isscalar(unique(delta)))
        error('seizmo:fkhorzvolume:badData',...
            'Records must have equal NPTS & DELTA fields!');
    end
    npts=npts(1); delta=delta(1);
    
    % check nyquist
    fnyq=1/(2*delta);
    if(any(frng>=fnyq))
        error('seizmo:fkhorzvolume:badFRNG',...
            ['FRNG frequencies must be under the nyquist frequency (' ...
            num2str(fnyq) ')!']);
    end
    
    % fix method/center/npairs
    if(ischar(method))
        method=lower(method);
        [clat,clon]=arraycenter(stla,stlo);
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

    % check that number of weights is equal to number of pairs
    if(isempty(w)); w=ones(npairs,1)/npairs; end
    if(numel(w)~=npairs)
        error('seizmo:fkhorzvolume:badInput',...
            ['WEIGHTS must have ' num2str(npairs) ...
            ' elements for method: ' method '!']);
    end
    
    % setup output
    [rvol(1:nrng,1).nsta]=deal(nrecs);
    [rvol(1:nrng,1).stla]=deal(est(:,1));
    [rvol(1:nrng,1).stlo]=deal(est(:,2));
    [rvol(1:nrng,1).stel]=deal(est(:,3));
    [rvol(1:nrng,1).stdp]=deal(est(:,4));
    [rvol(1:nrng,1).butc]=deal(butc(1,:));
    [rvol(1:nrng,1).eutc]=deal(eutc(1,:));
    [rvol(1:nrng,1).delta]=deal(delta);
    [rvol(1:nrng,1).npts]=deal(npts);
    [rvol(1:nrng,1).polar]=deal(polar);
    [rvol(1:nrng,1).npairs]=deal(npairs);
    [rvol(1:nrng,1).method]=deal(method);
    [rvol(1:nrng,1).center]=deal([clat clon]);
    [rvol(1:nrng,1).volume]=deal(true);
    [rvol(1:nrng,1).weights]=deal(w);
    
    % get frequencies
    nspts=2^nextpow2(npts);
    f=(0:nspts/2)/(delta*nspts);  % only +freq
    
    % extract data (silently)
    seizmoverbose(false);
    edat=records2mat(edat);
    ndat=records2mat(ndat);
    seizmoverbose(verbose);
    
    % get fft
    edat=fft(edat,nspts,1);
    ndat=fft(ndat,nspts,1);
    
    % get relative positions for each pair
    % r=(x  ,y  )
    %     ij  ij
    %
    % position of j as seen from i
    % x is km east
    % y is km north
    %
    % r is a 2xNPAIRS matrix
    switch method
        case 'coarray'
            % [ r   r   ... r
            %    11  12      1N
            %   r   r   ... r
            %    21  22      2N
            %    .   .  .    .
            %    .   .   .   .
            %    .   .    .  .
            %   r   r   ... r   ]
            %    N1  N2      NN
            %
            % then we just use the
            % upper triangle of that
            [e,n]=geographic2enu(stla,stlo,0,clat,clon,0);
            e=e(:,ones(nrecs,1))'-e(:,ones(nrecs,1));
            n=n(:,ones(nrecs,1))'-n(:,ones(nrecs,1));
            idx=triu(true(nrecs),1);
            e=e(idx);
            n=n(idx);
        case 'full'
            % retain full coarray
            [e,n]=geographic2enu(stla,stlo,0,clat,clon,0);
            e=e(:,ones(nrecs,1))'-e(:,ones(nrecs,1));
            n=n(:,ones(nrecs,1))'-n(:,ones(nrecs,1));
        case 'center'
            % each record relative to array center
            [e,n]=geographic2enu(clat,clon,0,stla,stlo,0);
        otherwise % user
            % each record relative to define array center
            [e,n]=geographic2enu(clat,clon,0,stla,stlo,0);
    end
    r=[e(:) n(:)]';
    clear e n
    
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%% CODE NEEDS WORK BELOW HERE
    %%%%%%%%%%%%%%
    
    % get phasors corresponding to each wave speed+direction (slowness)
    % to slant stack in the frequency domain
    % p=2*pi*i*s*r
    %
    % where s is the slowness vector s=(s ,s ) and is NSLOWx2
    %                                    x  y
    %
    % x is sec/km east
    % y is sec/km north
    %
    % p is the projection of all slownesses onto all of the position
    % vectors (multiplied by 2*pi*i for efficiency reasons)
    %
    % Also get rotators for each slowness/pair set
    %
    % u=cos(theta)
    % v=sin(theta)
    %
    % where theta is the angle from the position vector to the slowness
    % vector with positive being in the counter-clockwise direction
    % 
    % ie. theta = atan2(sy,sx)-atan2(ry,rx)
    %
    % p,u,v are NSLOWxNPAIRS
    smax=smax/d2km;
    if(polar)
        if(numel(spts)==2)
            bazpts=spts(2);
            spts=spts(1);
        else
            bazpts=181;
        end
        smag=(0:spts-1)/(spts-1)*smax;
        [rvol(1:nrng,1).y]=deal(smag'*d2km);
        smag=smag(ones(bazpts,1),:)';
        baz=(0:bazpts-1)/(bazpts-1)*360*d2r;
        [rvol(1:nrng,1).x]=deal(baz/d2r);
        baz=baz(ones(spts,1),:);
        p=2*pi*1i*[smag(:).*sin(baz(:)) smag(:).*cos(baz(:))]*r;
        % u = cos theta
        % v = sin theta
        theta1=atan2(r(2,:),r(1,:));
        theta2=baz(:);
        theta=theta2(:,ones(ncorr,1))-theta1(ones(spts*bazpts,1),:);
        u=cos(theta);
        v=sin(theta);
        % fix for s==0 (accept all azimuths b/c no directional dependance
        % in beam on vertical traveling waves)
        zeroslow=smag(:)==0;
        u(zeroslow,:)=1;
        v(zeroslow,:)=1;
        clear r smag baz zeroslow theta theta1 theta2
    else % cartesian
        spts=spts(1); bazpts=spts;
        sx=-smax:2*smax/(spts-1):smax;
        [rvol(1:nrng,1).x]=deal(sx*d2km);
        [rvol(1:nrng,1).y]=deal(fliplr(sx*d2km)');
        sx=sx(ones(spts,1),:);
        sy=fliplr(sx)';
        p=2*pi*1i*[sx(:) sy(:)]*r;
        % u = cos theta
        % v = sin theta
        theta1=atan2(r(2,:),r(1,:));
        theta2=atan2(sy(:),sx(:));
        theta=theta2(:,ones(ncorr,1))-theta1(ones(spts^2,1),:);
        u=cos(theta);
        v=sin(theta);
        % fix for s==0 (accept all azimuths b/c no directional dependance
        % in beam on vertical traveling waves)
        zeroslow=sy(:)==0 & sx(:)==0;
        u(zeroslow,:)=1;
        v(zeroslow,:)=1;
        clear r sx sy zeroslow theta theta1 theta2
    end
    
    % copy rvol to tvol
    tvol=rvol;
    
    % loop over frequency ranges
    for a=1:nrng
        % get frequencies
        fidx=find(f>=frng(a,1) & f<=frng(a,2));
        rvol(a).freq=f(fidx);
        tvol(a).freq=f(fidx);
        nfreq=numel(fidx);
        
        % preallocate fk space
        rvol(a).beam=zeros(spts,bazpts,nfreq,'single');
        tvol(a).beam=zeros(spts,bazpts,nfreq,'single');
        
        % warning if no frequencies
        if(~nfreq)
            warning('seizmo:fkxchorzvolume:noFreqs',...
                'No frequencies within the range %g to %g Hz!',...
                frng(a,1),frng(a,2));
            continue;
        end
        
        % detail message
        if(verbose)
            fprintf('Getting fk Volume %d for %g to %g Hz\n',...
                a,frng(a,1),frng(a,2));
            print_time_left(0,nfreq);
        end
        
        % loop over frequencies
        for b=1:nfreq
            % current freq idx
            cf=fidx(b);
            
            % rotate data into direction of plane wave for every pairing
            data=u.*u.*rr(cf*ones(1,spts*bazpts),:) ...
                -u.*v.*rt(cf*ones(1,spts*bazpts),:) ...
                -v.*u.*tr(cf*ones(1,spts*bazpts),:) ...
                +v.*v.*tt(cf*ones(1,spts*bazpts),:);
            
            % normalize by auto spectra
            data=data./abs(data);
            
            % apply weights
            data=data.*w(ones(1,spts*bazpts),:);
            
            % now getting fk beam for radial (rayleigh)
            rvol(a).beam(:,:,b)=reshape(10*log10(abs(real(...
                sum(data.*exp(f(cf)*p),2)))),spts,bazpts);
            
            % rotate data perpendicular to plane wave direction for all
            data=v.*v.*rr(cf*ones(1,spts*bazpts),:) ...
                +v.*u.*rt(cf*ones(1,spts*bazpts),:) ...
                +u.*v.*tr(cf*ones(1,spts*bazpts),:) ...
                +u.*u.*tt(cf*ones(1,spts*bazpts),:);
            
            % normalize by auto spectra
            data=data./abs(data);
            
            % apply weights
            data=data.*w(ones(1,spts*bazpts),:);
            
            % now getting fk beam for tangential (love)
            tvol(a).beam(:,:,b)=reshape(10*log10(abs(real(...
                sum(data.*exp(f(cf)*p),2)))),spts,bazpts);
            
            % detail message
            if(verbose); print_time_left(b,nfreq); end
        end
        
        % normalize so max peak is at 0dB
        rvol(a).normdb=max(rvol(a).beam(:));
        rvol(a).beam=rvol(a).beam-rvol(a).normdb;
        tvol(a).normdb=max(tvol(a).beam(:));
        tvol(a).beam=tvol(a).beam-tvol(a).normdb;
        
        % plot if no output
        if(~nargout)
            fkfreqslide(rvol(a));
            fkfreqslide(tvol(a));
        end
    end
    
    % return struct
    if(nargout); varargout{1}=rvol; varargout{2}=tvol; end
    
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
