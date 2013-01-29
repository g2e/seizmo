function [varargout]=fkxcmap(data,smax,spts,frng,polar,w)
%FKXCMAP    Returns beam map in frequency-wavenumber space for xc data
%
%    Usage:    smap=fkxcmap(xcdata,smax,spts,frng)
%              smap=fkxcmap(xcdata,smax,spts,frng,polar)
%              smap=fkxcmap(xcdata,smax,spts,frng,polar,weights)
%
%    Description:
%     SMAP=FKXCMAP(XCDATA,SMAX,SPTS,FRNG) beamforms the wave
%     energy moving through an array in frequency-wavenumber space.
%     Actually, to allow for easier interpretation between frequencies,
%     the energy is mapped into frequency-slowness space.  The array info
%     and correlograms are contained in the SEIZMO struct XCDATA.  This
%     differs from FKMAP in that DATA is expected to be correlograms
%     resulting from a multi-channel cross correlation with CORRELATE
%     rather than actual records.  This has the advantage of allowing time
%     averaging (aka stacking) of the cross spectra to improve resolution
%     of persistent sources of seismic energy like microseismic energy from
%     storms over the oceans.  The range of the slowness space is given by
%     SMAX (in s/deg) and extends from -SMAX to SMAX for both East/West and
%     North/South directions.  SPTS controls the number of slowness points
%     for both directions (SPTSxSPTS grid).  FRNG gives the frequency range
%     as [FREQLOW FREQHIGH] in Hz.  SMAP is a struct containing relevant
%     info and the frequency-slowness volume itself.  The struct layout is:
%          .beam     - frequency-slowness beamforming map
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
%          .npairs   - number of pairs (aka correlograms)
%          .method   - beamforming method (center, coarray, full, user)
%          .center   - array center as [lat lon]
%          .normdb   - what 0dB actually corresponds to
%          .volume   - true if frequency-slowness volume (false for FKMAP)
%          .weights  - weights used in beamforming
%
%     Calling FKXCMAP with no outputs will automatically plot the slowness
%     map using PLOTFKMAP.
%
%     SMAP=FKXCMAP(XCDATA,SMAX,SPTS,FRNG,POLAR) specifies if the slowness
%     space is sampled regularly in cartesian or polar coordinates.  Polar
%     coords are useful for slicing the volume by azimuth (pie slice) or
%     slowness magnitude (rings).  Cartesian coords (the default) samples
%     the slowness space regularly in the East/West & North/South
%     directions and so exhibits less distortion of the slowness space.
%     
%     SMAP=FKXCMAP(XCDATA,SMAX,SPTS,FRNG,POLAR,WEIGHTS) specifies
%     weights for each correlogram in XCDATA (must match size of XCDATA)
%     for use in beamforming.  The weights are normalized internally to sum
%     to 1.
%
%    Notes:
%     - Records in XCDATA must be correlograms following the formatting
%       from CORRELATE.  The records must have the same lag range & sample
%       spacing.
%     - Best/quickest results are obtained when XCDATA is only one
%       "triangle" of the cross correlation matrix.  This corresponds to
%       the 'coarray' method from FKMAP.
%
%    Examples:
%     % Show frequency-slowness map for a dataset at 20-50s periods:
%     fkxcmap(correlate(data,'mcxc','noauto'),50,201,[1/50 1/20]);
%
%    See also: FKMAP, FKFREQSLIDE, FKMAP, FK4D, FKVOL2MAP, FKSUBVOL,
%              FKXCHORZMAP, CORRELATE

%     Version History:
%        Nov. 18, 2010 - initial version
%        Apr.  3, 2012 - use seizmocheck
%        Jan. 29, 2013 - example update for new correlate functions
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 29, 2013 at 15:25 GMT

% todo:

% check nargin
error(nargchk(4,6,nargin));

% define some constants
d2r=pi/180;
d2km=6371*d2r;

% check struct
error(seizmocheck(data,'dep'));

% number of correlograms
ncorr=numel(data);

% defaults for optionals
if(nargin<5 || isempty(polar)); polar=false; end
if(nargin<6 || isempty(w)); w=ones(numel(data),1); end
method='coarray';

% check inputs
sf=size(frng);
if(~isreal(smax) || ~isscalar(smax) || smax<=0)
    error('seizmo:fkxcmap:badInput',...
        'SMAX must be a positive real scalar in s/deg!');
elseif(~any(numel(spts)==[1 2]) || any(fix(spts)~=spts) || any(spts<=2))
    error('seizmo:fkxcmap:badInput',...
        'SPTS must be a positive scalar integer >2!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:fkxcmap:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
elseif(~isscalar(polar) || (~islogical(polar) && ~isnumeric(polar)))
    error('seizmo:fkxcmap:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif(numel(w)~=ncorr || any(w(:)<0) ||  ~isreal(w) || sum(w(:))==0)
    error('seizmo:fkxcmap:badInput',...
        'WEIGHTS must be equal sized with XCDATA & be positive numbers!');
end
nrng=sf(1);

% convert weights to column vector and normalize
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

% do fk analysis
try
    % verbosity
    verbose=seizmoverbose;
    
    % get some useful info from headers
    [npts,delta,st,ev]=getheader(data,'npts','delta','st','ev');
    
    % find unique station locations
    loc=unique([st; ev],'rows');
    nsta=size(loc,1);
    
    % array center
    [clat,clon]=arraycenter([st(:,1); ev(:,1)],[st(:,2); ev(:,2)]);
    
    % check nyquist
    fnyq=1/(2*delta(1));
    if(any(frng>=fnyq))
        error('seizmo:fkxcmap:badFRNG',...
            ['FRNG frequencies must be under the nyquist frequency (' ...
            num2str(fnyq) ')!']);
    end
    
    % setup output
    [smap(1:nrng,1).nsta]=deal(nsta);
    [smap(1:nrng,1).stla]=deal(loc(:,1));
    [smap(1:nrng,1).stlo]=deal(loc(:,2));
    [smap(1:nrng,1).stel]=deal(loc(:,3));
    [smap(1:nrng,1).stdp]=deal(loc(:,4));
    [smap(1:nrng,1).butc]=deal([0 0 0 0 0]);
    [smap(1:nrng,1).eutc]=deal([0 0 0 0 0]);
    [smap(1:nrng,1).delta]=deal(delta(1));
    [smap(1:nrng,1).npts]=deal(npts(1));
    [smap(1:nrng,1).polar]=deal(polar);
    [smap(1:nrng,1).npairs]=deal(ncorr);
    [smap(1:nrng,1).method]=deal(method);
    [smap(1:nrng,1).center]=deal([clat clon]);
    [smap(1:nrng,1).volume]=deal(false);
    [smap(1:nrng,1).weights]=deal(w);
    
    % get frequencies (note no extra power for correlations)
    nspts=2^nextpow2(npts(1));
    f=(0:nspts/2)/(delta(1)*nspts);  % only +freq
    
    % extract data (silently)
    seizmoverbose(false);
    data=splitpad(data,0);
    data=records2mat(data);
    seizmoverbose(verbose);
    
    % get fft (conjugate is b/c my xc is flipped?)
    % - this is the true cross spectra
    data=conj(fft(data,nspts,1));
    
    % get relative positions for each pair
    % r=(x  ,y  )
    %     ij  ij
    %
    % position of j as seen from i
    % x is km east
    % y is km north
    %
    % r is a 2xNCORR matrix
    [e_ev,n_ev]=geographic2enu(ev(:,1),ev(:,2),0,clat,clon,0);
    [e_st,n_st]=geographic2enu(st(:,1),st(:,2),0,clat,clon,0);
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
    smax=smax/d2km;
    if(polar)
        if(numel(spts)==2)
            bazpts=spts(2);
            spts=spts(1);
        else
            bazpts=181;
        end
        smag=(0:spts-1)/(spts-1)*smax;
        [smap(1:nrng,1).y]=deal(smag'*d2km);
        smag=smag(ones(bazpts,1),:)';
        baz=(0:bazpts-1)/(bazpts-1)*360*d2r;
        [smap(1:nrng,1).x]=deal(baz/d2r);
        baz=baz(ones(spts,1),:);
        p=2*pi*1i*[smag(:).*sin(baz(:)) smag(:).*cos(baz(:))]*r;
        clear r smag baz
        [smap(1:nrng,1).beam]=deal(zeros(spts,bazpts,'single'));
    else % cartesian
        spts=spts(1);
        sx=-smax:2*smax/(spts-1):smax;
        [smap(1:nrng,1).x]=deal(sx*d2km);
        [smap(1:nrng,1).y]=deal(fliplr(sx*d2km)');
        sx=sx(ones(spts,1),:);
        sy=fliplr(sx)';
        p=2*pi*1i*[sx(:) sy(:)]*r;
        clear r sx sy
        [smap(1:nrng,1).beam]=deal(zeros(spts,'single'));
    end
    
    % loop over frequency ranges
    for a=1:nrng
        % get frequencies
        fidx=find(f>=frng(a,1) & f<=frng(a,2));
        smap(a).freq=f(fidx);
        nfreq=numel(fidx);
        
        % warning if no frequencies
        if(~nfreq)
            warning('seizmo:fkxcmap:noFreqs',...
                'No frequencies within the range %g to %g Hz!',...
                frng(a,1),frng(a,2));
            continue;
        end
        
        % extract complex cross spectra array
        % cs is normalized by the auto spectra.
        %
        % cs is NCORRxNFREQ
        cs=data(fidx,:).';
        cs=cs./abs(cs);
        
        % apply weighting
        cs=cs.*w(:,ones(1,nfreq));
        
        % detail message
        if(verbose)
            fprintf('Getting fk Map %d for %g to %g Hz\n',...
                a,frng(a,1),frng(a,2));
            print_time_left(0,nfreq);
        end
        
        % loop over frequencies
        for b=1:nfreq
            % form beam
            smap(a).beam(:)=smap(a).beam(:)+exp(f(fidx(b))*p)*cs(:,b);
            
            % detail message
            if(verbose); print_time_left(b,nfreq); end
        end
        
        % using full and real here gives the exact plots of
        % Koper, Seats, and Benz 2010 in BSSA
        smap(a).beam=10*log10(abs(real(smap(a).beam))/nfreq);
        
        % normalize so max peak is at 0dB
        smap(a).normdb=max(smap(a).beam(:));
        smap(a).beam=smap(a).beam-smap(a).normdb;
        
        % plot if no output
        if(~nargout); plotfkmap(smap(a)); drawnow; end
    end
    
    % return struct
    if(nargout); varargout{1}=smap; end
    
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
