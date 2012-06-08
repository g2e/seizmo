function [rvol,tvol]=geofkxchorzvolume(rr,rt,tr,tt,ll,s,frng,w)
%GEOFKXCHORZVOLUME    Geographic FK beamforming of horizontals
%
%    Usage:    [rgeo,tgeo]=geofkxchorzvolume(rr,rt,tr,tt,...
%                                            latlon,horzslow,frng)
%              [rgeo,tgeo]=geofkxchorzvolume(rr,rt,tr,tt,...
%                                            latlon,horzslow,frng,weights)
%
%    Description:
%     [RGEO,TGEO]=GEOFKXCHORZVOLUME(RR,RT,TR,TT,LATLON,HORZSLOW,FRNG)
%     computes the spherical wave coherency (beam strength) through an
%     array as a function of frequency, horizontal slowness & geographic
%     location. The array info and correlograms between sets of horizontal
%     components are contained in the SEIZMO structs RR, RT, TR, & TT.
%     This differs from FKXCHORZVOLUME in that the waves are assumed to be
%     spherical rather than planar.  This is essential for handling sources
%     that are near the array (my rule of thumb is that if it is within 1
%     array aperture width of the array then a near-field based fk analysis
%     should be done - particularly with horizontals).  RR, RT, TR, & TT
%     must contain correlograms with header formatting as output from
%     CORRELATE.  LATLON contains the latitude and longitude wavefield
%     source positions and must be in units of degrees and formatted as an
%     Nx2 array of [LAT LON].  HORZSLOW is the magnitude of the horizontal
%     slowness of the waves in sec/deg and should be a vector.  FRNG gives
%     the frequency range as [FREQLOW FREQHIGH] in Hz.  The outputs RGEO &
%     TGEO correspond to radially and transversely polarized spherical
%     wavefronts respectively.  They are structs containing relevant info
%     and the frequency-slowness-position volume itself (with size
%     NPOSxNSLOWxNFREQ).  The struct layout is:
%          .beam     - frequency-slowness-position beamforming volume
%          .nsta     - number of stations utilized in making map
%          .stla     - station latitudes
%          .stlo     - station longitudes
%          .stel     - station elevations (surface)
%          .stdp     - station depths (from surface)
%          .butc     - UTC start time of data (set to all zeros)
%          .eutc     - UTC end time of data (set to all zeros)
%          .npts     - number of time points
%          .delta    - number of seconds between each time point
%          .latlon   - latitude/longitude positions (deg)
%          .horzslow - horizontal slowness (sec/deg)
%          .freq     - frequency values (Hz)
%          .npairs   - number of pairs (aka correlograms)
%          .method   - beamforming method (user, center, coarray, full)
%          .center   - array center as [LAT LON]
%          .normdb   - what 0dB actually corresponds to
%          .volume   - true if frequency-slowness-position volume
%          .weights  - weights used in beamforming
%
%     SGEO=GEOFKXCVOLUME(DATA,LATLON,HORZSLOW,FRNG,WEIGHTS) specifies
%     weights for each correlogram in XCDATA (must match size of XCDATA).
%     The weights are normalized internally to sum to 1.
%
%    Notes:
%     - Records in RR, RT, TR, & TT must be correlograms following the
%       formatting from CORRELATE.  The records must have the same lag
%       range & sample spacing.  Furthermore, the datasets must be equal
%       sized and every record should correspond to the records with the
%       same index in the other datasets (ie RR(3), RT(3), TR(3) & TT(3)
%       must correspond to the same station pair).
%     - Best/quickest results are obtained when RR, RT, TR, & TT is one
%       "triangle" of the cross correlation matrix.  This corresponds to
%       the 'coarray' method.
%
%    Examples:
%
%    See also: FKXCHORZVOLUME, GEOFKXCVOLUME, CHKGEOFKSTRUCT, PLOTGEOFKMAP,
%              GEOFKFREQSLIDE, GEOFKSLOWSLIDE, GEOFKSUBVOL, GEOFKVOL2MAP

%     Version History:
%        June 22, 2010 - initial version
%        July  1, 2010 - freq field bug fix, d2r bug fixed
%        July  6, 2010 - major update to struct, doc update
%        July  7, 2010 - removed deg to km conversions
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July  7, 2010 at 15:15 GMT

% todo:

% check nargin/nargout
error(nargchk(7,8,nargin));
error(nargchk(1,2,nargout));

% define some constants
d2r=pi/180;

% check struct
versioninfo(rr,'dep');
versioninfo(rt,'dep');
versioninfo(tr,'dep');
versioninfo(tt,'dep');

% make sure rr/tt are the same size
ncorr=numel(rr);
ncorr1=numel(tt);
ncorr2=numel(tt);
ncorr3=numel(tt);
if(~isequal(ncorr,ncorr1,ncorr2,ncorr3))
    error('seizmo:fkxchorzvolume:unmatchedXCdata',...
        'XC datasets do not match in size!');
end

% defaults for optionals
if(nargin<8 || isempty(w)); w=ones(ncorr,1); end
method='coarray';

% check inputs
sf=size(frng);
if(~isreal(ll) || ndims(ll)~=2 || size(ll,2)~=2)
    error('seizmo:geofkxchorzvolume:badInput',...
        'LATLON must be a Nx2 positive real matrix of [LAT LON]!');
elseif(~isreal(s) || ~isvector(s) || any(s<=0))
    error('seizmo:geofkxchorzvolume:badInput',...
        'HORZSLOW must be a positive real vector in sec/deg!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:geofkxchorzvolume:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
elseif(numel(w)~=ncorr || any(w(:)<0) || ~isreal(w) || sum(w(:))==0)
    error('seizmo:geofkxchorzvolume:badInput',...
        'WEIGHTS must be equal sized with XCDATA & be positive numbers!');
end
nrng=sf(1);

% column vector slownesses
s=s(:);
nslow=numel(s);

% fix lat/lon
[ll(:,1),ll(:,2)]=fixlatlon(ll(:,1),ll(:,2));
nll=size(ll,1);

% convert weights to row vector
w=w(:).';
w=w./sum(w);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    rr=checkheader(rr,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'MULTIPLE_NPTS','ERROR',...
        'MULTIPLE_B','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    rt=checkheader(rt,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'MULTIPLE_NPTS','ERROR',...
        'MULTIPLE_B','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    tr=checkheader(tr,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'MULTIPLE_NPTS','ERROR',...
        'MULTIPLE_B','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    tt=checkheader(tt,...
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

% do fk analysis
try
    % verbosity
    verbose=seizmoverbose;
    
    % require radial & transverse components
    [rrmcn,rrscn]=getheader(rr,'kt3','kcmpnm');
    [rtmcn,rtscn]=getheader(rt,'kt3','kcmpnm');
    [trmcn,trscn]=getheader(tr,'kt3','kcmpnm');
    [ttmcn,ttscn]=getheader(tt,'kt3','kcmpnm');
    rrmcn=char(rrmcn); rrmcn=rrmcn(:,3);
    rrscn=char(rrscn); rrscn=rrscn(:,3);
    rtmcn=char(rtmcn); rtmcn=rtmcn(:,3);
    rtscn=char(rtscn); rtscn=rtscn(:,3);
    trmcn=char(trmcn); trmcn=trmcn(:,3);
    trscn=char(trscn); trscn=trscn(:,3);
    ttmcn=char(ttmcn); ttmcn=ttmcn(:,3);
    ttscn=char(ttscn); ttscn=ttscn(:,3);
    if(~isequal(lower(unique(rrmcn)),'r') ...
            || ~isequal(lower(unique(rrscn)),'r'))
        error('seizmo:geofkxchorzvolume:badRR',...
            'RR does not appear to be Radial-Radial XC data!');
    elseif(~isequal(lower(unique(rtmcn)),'r') ...
            || ~isequal(lower(unique(rtscn)),'t'))
        error('seizmo:geofkxchorzvolume:badRT',...
            'RT does not appear to be Radial-Transverse XC data!');
    elseif(~isequal(lower(unique(trmcn)),'t') ...
            || ~isequal(lower(unique(trscn)),'r'))
        error('seizmo:geofkxchorzvolume:badTR',...
            'TR does not appear to be Transverse-Radial XC data!');
    elseif(~isequal(lower(unique(ttmcn)),'t') ...
            || ~isequal(lower(unique(ttscn)),'t'))
        error('seizmo:geofkxchorzvolume:badTT',...
            'TT does not appear to be Transverse-Transverse XC data!');
    end
    
    % get station locations & require that they match
    [st1,ev1,az,baz]=getheader(rr,'st','ev','az','baz');
    [st2,ev2]=getheader(rt,'st','ev');
    [st3,ev3]=getheader(tr,'st','ev');
    [st4,ev4]=getheader(tt,'st','ev');
    if(~isequal([st1; ev1],[st2; ev2],[st3; ev3],[st4; ev4]))
        error('seizmo:geofkxchorzvolume:xcDataMismatch',...
            'XC Datasets have different station locations!');
    end
    clear st2 st3 st4 ev2 ev3 ev4
    
    % find unique station locations
    loc=unique([st1; ev1],'rows');
    nsta=size(loc,1);
    
    % array center
    [clat,clon]=arraycenter([st1(:,1); ev1(:,1)],[st1(:,2); ev1(:,2)]);
    
    % require all records to have equal b, npts, delta
    [b,npts,delta]=getheader([rr(:); rt(:); tr(:); tt(:)],...
        'b','npts','delta');
    if(~isscalar(unique(b)) ...
            || ~isscalar(unique(npts)) ...
            || ~isscalar(unique(delta)))
        error('seizmo:geofkxchorzvolume:badData',...
            'XC records must have equal B, NPTS & DELTA fields!');
    end
    npts=npts(1); delta=delta(1);
    
    % check nyquist
    fnyq=1/(2*delta(1));
    if(any(frng>=fnyq))
        error('seizmo:geofkxchorzvolume:badFRNG',...
            ['FRNG frequencies must be under the nyquist frequency (' ...
            num2str(fnyq) ')!']);
    end
    
    % setup output
    [rvol(1:nrng,1).nsta]=deal(nsta);
    [rvol(1:nrng,1).stla]=deal(loc(:,1));
    [rvol(1:nrng,1).stlo]=deal(loc(:,2));
    [rvol(1:nrng,1).stel]=deal(loc(:,3));
    [rvol(1:nrng,1).stdp]=deal(loc(:,4));
    [rvol(1:nrng,1).butc]=deal([0 0 0 0 0]);
    [rvol(1:nrng,1).eutc]=deal([0 0 0 0 0]);
    [rvol(1:nrng,1).delta]=deal(delta);
    [rvol(1:nrng,1).npts]=deal(npts);
    [rvol(1:nrng,1).volume]=deal(true(1,2));
    [rvol(1:nrng,1).latlon]=deal(ll);
    [rvol(1:nrng,1).horzslow]=deal(s);
    [rvol(1:nrng,1).npairs]=deal(ncorr);
    [rvol(1:nrng,1).method]=deal(method);
    [rvol(1:nrng,1).center]=deal([clat clon]);
    [rvol(1:nrng,1).weights]=deal(w);
    
    % get frequencies (note no extra power for correlations)
    nspts=2^nextpow2(npts);
    f=(0:nspts/2)/(delta*nspts);  % only +freq
    
    % extract data (silently)
    seizmoverbose(false);
    rr=splitpad(rr);
    rr=records2mat(rr);
    rt=splitpad(rt);
    rt=records2mat(rt);
    tr=splitpad(tr);
    tr=records2mat(tr);
    tt=splitpad(tt);
    tt=records2mat(tt);
    seizmoverbose(verbose);
    
    % get fft (conjugate is b/c my xc is flipped?)
    % - this is the true cross spectra
    rr=conj(fft(rr,nspts,1));
    rt=conj(fft(rt,nspts,1));
    tr=conj(fft(tr,nspts,1));
    tt=conj(fft(tt,nspts,1));
    
    % distance difference for the phasors that steer the array
    % dd is NLLxNCORR
    % - also getting azimuth info for rotators
    ev1=ev1.'; st1=st1.';
    [distm,azm,bazm]=sphericalinv(ll(:,ones(ncorr,1)),...
        ll(:,2*ones(ncorr,1)),ev1(ones(nll,1),:),ev1(2*ones(nll,1),:));
    [dists,azs,bazs]=sphericalinv(ll(:,ones(ncorr,1)),...
        ll(:,2*ones(ncorr,1)),st1(ones(nll,1),:),st1(2*ones(nll,1),:));
    dd=2*pi*1i*(distm-dists);
    
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
        rvol(a).beam=zeros(nll,nslow,nfreq,'single');
        tvol(a).beam=zeros(nll,nslow,nfreq,'single');
        
        % warning if no frequencies
        if(~nfreq)
            warning('seizmo:geofkxchorzvolume:noFreqs',...
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
