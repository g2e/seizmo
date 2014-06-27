function [varargout]=fkxchorzvolume(rr,rt,tr,tt,smax,spts,frng,polar,w)
%FKXCHORZVOLUME    Returns frequency-wavenumber space for horz. xc data
%
%    Usage:    [rvol,tvol]=fkxchorzvolume(rr,rt,tr,tt,smax,spts,frng)
%              [rvol,tvol]=fkxchorzvolume(rr,rt,tr,tt,smax,spts,frng,polar)
%              [rvol,tvol]=fkxchorzvolume(rr,rt,tr,tt,smax,spts,frng,...
%                                         polar,weights)
%
%    Description:
%     [RVOL,TVOL]=FKXCHORZVOLUME(RR,RT,TR,TT,SMAX,SPTS,FRNG) calculates
%     the Rayleigh & Love energy moving through an array in
%     frequency-wavenumber space utilizing the horizontal cross correlation
%     datasets RR, RT, TR & TT.  To allow for easier interpretation between
%     frequencies, the energy is mapped into frequency-slowness space.  The
%     array info and correlograms are contained in the SEIZMO structs RR,
%     RT, TR, & TT.  RR is expected to contain the pairwise radial-radial
%     correlations and TT is expected to contain the pairwise transverse
%     -transverse correlations (RT/TR give the radial-tranverse & the
%     transverse-radial components - energy on these records is indicative
%     of an anisotropic source distribution).  This also differs from
%     FKXCVOLUME in that horizontals are utilized to retreive both the
%     radial (RVOL) & transverse (TVOL) energy distributions.  FKXCVOLUME
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
%     [RVOL,TVOL]=FKXCHORZVOLUME(RR,RT,TR,TT,SMAX,SPTS,FRNG,POLAR) sets
%     if the slowness space is sampled regularly in cartesian or polar
%     coordinates.  Polar coords are useful for slicing the volume by
%     azimuth (pie slice) or slowness magnitude (rings).  Cartesian coords
%     (the default) samples the slowness space regularly in the East/West
%     & North/South directions and so exhibits less distortion of the
%     slowness space.
%
%     [RVOL,TVOL]=FKXCHORZVOLUME(RR,RT,TR,TT,SMAX,SPTS,FRNG,POLAR,WEIGHTS)
%     specifies weights for each correlogram in RR/RT/TR/TT (must match
%     size of RR all) for use in beamforming.  The weights are normalized
%     internally to sum to 1.
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
%       the 'coarray' method from FKVOLUME.
%
%    Examples:
%     % One way to perform horizontal fk analysis:
%     data=rotate(data,'to',0,'kcmpnm1','N','kcmpnm2','E');
%     xc=correlate(data,'mcxc');
%     [axc,xc]=split_auto_correlations(xc);
%     [b,e,delta]=getheader(xc,'b','e','delta');
%     xc=interpolate(xc,delta(1),'spline',max(b),min(e));
%     xc=rotate_correlations(xc);
%     [in,set,cmp]=horz_correlations_sets(xc);
%     [rr,rt,tr,tt]=deal(xc(cmp==1),xc(cmp==2),xc(cmp==3),xc(cmp==4));
%     [r,t]=fkxchorzvolume(rr,rt,tr,tt,50,101,[1/50 1/20]);
%     fkfreqslide(rvol,0);
%     fkfreqslide(tvol,0);
%
%    See also: FKXCVOLUME, FKFREQSLIDE, FKMAP, FK4D, FKVOL2MAP, FKSUBVOL,
%              FKVOLUME, CORRELATE, ROTATE_CORRELATIONS

%     Version History:
%        June  9, 2010 - initial version
%        June 12, 2010 - math is sound now
%        June 13, 2010 - passes several eq verifications
%        June 16, 2010 - allow any form of cross correlation matrix (not
%                        just one triangle), doc update, add example
%        June 18, 2010 - add weights
%        June 22, 2010 - default weights
%        July  1, 2010 - high latitude fix
%        July  6, 2010 - major update to struct, doc update
%        Nov. 18, 2010 - .weights bugfix
%        Jan. 29, 2013 - example update for new correlate functions
%        Sep. 23, 2013 - example update for new correlate functions,
%                        adjust for rotate_correlations change
%        May  30, 2014 - adjusted example for split_auto_correlations fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  30, 2014 at 15:25 GMT

% todo:

% check nargin
error(nargchk(7,9,nargin));

% define some constants
d2r=pi/180;
d2km=6371*d2r;

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
if(nargin<8 || isempty(polar)); polar=false; end
if(nargin<9 || isempty(w)); w=ones(ncorr,1); end
method='coarray';

% check inputs
sf=size(frng);
if(~isreal(smax) || ~isscalar(smax) || smax<=0)
    error('seizmo:fkxchorzvolume:badInput',...
        'SMAX must be a positive real scalar in s/deg!');
elseif(~any(numel(spts)==[1 2]) || any(fix(spts)~=spts) || any(spts<=2))
    error('seizmo:fkxchorzvolume:badInput',...
        'SPTS must be a positive scalar integer >2!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:fkxchorzvolume:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
elseif(~isscalar(polar) || (~islogical(polar) && ~isnumeric(polar)))
    error('seizmo:fkxchorzvolume:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif(numel(w)~=ncorr || any(w(:)<0) || ~isreal(w) || sum(w(:))==0)
    error('seizmofkxchorzvolume:badInput',...
        'WEIGHTS must be equal sized with XCDATA & be positive numbers!');
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
        error('seizmo:fkxchorzvolume:badRR',...
            'RR does not appear to be Radial-Radial XC data!');
    elseif(~isequal(lower(unique(rtmcn)),'r') ...
            || ~isequal(lower(unique(rtscn)),'t'))
        error('seizmo:fkxchorzvolume:badRT',...
            'RT does not appear to be Radial-Transverse XC data!');
    elseif(~isequal(lower(unique(trmcn)),'t') ...
            || ~isequal(lower(unique(trscn)),'r'))
        error('seizmo:fkxchorzvolume:badTR',...
            'TR does not appear to be Transverse-Radial XC data!');
    elseif(~isequal(lower(unique(ttmcn)),'t') ...
            || ~isequal(lower(unique(ttscn)),'t'))
        error('seizmo:fkxchorzvolume:badTT',...
            'TT does not appear to be Transverse-Transverse XC data!');
    end
    
    % get station locations & require that they match
    [st1,ev1]=getheader(rr,'st','ev');
    [st2,ev2]=getheader(rt,'st','ev');
    [st3,ev3]=getheader(tr,'st','ev');
    [st4,ev4]=getheader(tt,'st','ev');
    if(~isequal([st1; ev1],[st2; ev2],[st3; ev3],[st4; ev4]))
        error('seizmo:fkxchorzvolume:xcDataMismatch',...
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
        error('seizmo:fkxchorzvolume:badData',...
            'XC records must have equal B, NPTS & DELTA fields!');
    end
    npts=npts(1); delta=delta(1);
    
    % check nyquist
    fnyq=1/(2*delta(1));
    if(any(frng>=fnyq))
        error('seizmo:fkxchorzvolume:badFRNG',...
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
    [rvol(1:nrng,1).polar]=deal(polar);
    [rvol(1:nrng,1).npairs]=deal(ncorr);
    [rvol(1:nrng,1).method]=deal(method);
    [rvol(1:nrng,1).center]=deal([clat clon]);
    [rvol(1:nrng,1).volume]=deal(true);
    [rvol(1:nrng,1).weights]=deal(w);
    
    % get frequencies (note no extra power for correlations)
    nspts=2^nextpow2(npts);
    f=(0:nspts/2)/(delta*nspts);  % only +freq
    
    % extract data (silently)
    seizmoverbose(false);
    rr=splitpad(rr,0);
    rr=records2mat(rr);
    rt=splitpad(rt,0);
    rt=records2mat(rt);
    tr=splitpad(tr,0);
    tr=records2mat(tr);
    tt=splitpad(tt,0);
    tt=records2mat(tt);
    seizmoverbose(verbose);
    
    % get fft (conjugate is b/c my xc is flipped?)
    % - this is the true cross spectra
    rr=conj(fft(rr,nspts,1));
    rt=conj(fft(rt,nspts,1));
    tr=conj(fft(tr,nspts,1));
    tt=conj(fft(tt,nspts,1));
    
    % get relative positions for each pair
    % r=(x  ,y  )
    %     ij  ij
    %
    % position of j as seen from i
    % x is km east
    % y is km north
    %
    % r is a 2xNCORR matrix
    [e_ev,n_ev]=geographic2enu(ev1(:,1),ev1(:,2),0,clat,clon,0);
    [e_st,n_st]=geographic2enu(st1(:,1),st1(:,2),0,clat,clon,0);
    r=[e_st-e_ev n_st-n_ev]';
    clear e_ev e_st n_ev n_st
    
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
    % p is NSLOWxNPAIRS
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
    % p,u,v are NSLOWxNCORR
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
