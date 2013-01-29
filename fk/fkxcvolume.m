function [varargout]=fkxcvolume(data,smax,spts,frng,polar,w)
%FKXCVOLUME    Returns beam map in frequency-wavenumber space for xc data
%
%    Usage:    svol=fkxcvolume(xcdata,smax,spts,frng)
%              svol=fkxcvolume(xcdata,smax,spts,frng,polar)
%              svol=fkxcvolume(xcdata,smax,spts,frng,polar,weights)
%
%    Description:
%     SVOL=FKXCVOLUME(XCDATA,SMAX,SPTS,FRNG) beamforms the wave
%     energy moving through an array in frequency-wavenumber space.
%     Actually, to allow for easier interpretation between frequencies,
%     the energy is mapped into frequency-slowness space.  The array info
%     and correlograms are contained in the SEIZMO struct XCDATA.  This
%     differs from FKVOLUME in that DATA is expected to be correlograms
%     resulting from a multi-channel cross correlation with CORRELATE
%     rather than actual records.  This has the advantage of allowing time
%     averaging (aka stacking) of the cross spectra to improve resolution
%     of persistent sources of seismic energy like microseismic energy from
%     storms over the oceans.  The range of the slowness space is given by
%     SMAX (in s/deg) and extends from -SMAX to SMAX for both East/West and
%     North/South directions.  SPTS controls the number of slowness points
%     for both directions (SPTSxSPTS grid).  FRNG gives the frequency range
%     as [FREQLOW FREQHIGH] in Hz.  SVOL is a struct containing relevant
%     info and the frequency-slowness volume itself.  The struct layout is:
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
%          .npairs   - number of pairs (aka correlograms)
%          .method   - beamforming method (center, coarray, full, user)
%          .center   - array center as [lat lon]
%          .normdb   - what 0dB actually corresponds to
%          .volume   - true if frequency-slowness volume (false for FKMAP)
%          .weights  - weights used in beamforming
%
%     Calling FKXCVOLUME with no outputs will automatically plot the
%     frequency-slowness volume using FKFREQSLIDE.
%
%     SVOL=FKXCVOLUME(XCDATA,SMAX,SPTS,FRNG,POLAR) sets if the slowness
%     space is sampled regularly in cartesian or polar coordinates.  Polar
%     coords are useful for slicing the volume by azimuth (pie slice) or
%     slowness magnitude (rings).  Cartesian coords (the default) samples
%     the slowness space regularly in the East/West & North/South
%     directions and so exhibits less distortion of the slowness space.
%     
%     SVOL=FKXCVOLUME(XCDATA,SMAX,SPTS,FRNG,POLAR,WEIGHTS) specifies
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
%       the 'coarray' method from FKVOLUME.
%
%    Examples:
%     % Show frequency-slowness volume for a dataset at 20-50s periods:
%     xcdata=correlate(data,'mcxc','noauto');
%     svol=fkxcvolume(xcdata,50,201,[1/50 1/20]);
%     mov=fkfreqslide(svol,0); % show once with no delay (make movie too)
%     h=figure; % open figure so movie consumes figure
%     movie(h,mov,10); % show 10 times
%
%    See also: FKVOLUME, FKFREQSLIDE, FKMAP, FK4D, FKVOL2MAP, FKSUBVOL,
%              FKXCHORZVOLUME, CORRELATE

%     Version History:
%        May  13, 2010 - initial version
%        June  8, 2010 - minor touches
%        June 16, 2010 - allows any form of cross correlation matrix (not
%                        just one triangle), fix see also section
%        June 18, 2010 - add weights
%        July  1, 2010 - high latitude fix
%        July  6, 2010 - major update to struct, doc update
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
    error('seizmo:fkxcvolume:badInput',...
        'SMAX must be a positive real scalar in s/deg!');
elseif(~any(numel(spts)==[1 2]) || any(fix(spts)~=spts) || any(spts<=2))
    error('seizmo:fkxcvolume:badInput',...
        'SPTS must be a positive scalar integer >2!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:fkxcvolume:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
elseif(~isscalar(polar) || (~islogical(polar) && ~isnumeric(polar)))
    error('seizmo:fkxcvolume:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif(numel(w)~=ncorr || any(w(:)<0) ||  ~isreal(w) || sum(w(:))==0)
    error('seizmo:fkxcvolume:badInput',...
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
        error('seizmo:fkxcvolume:badFRNG',...
            ['FRNG frequencies must be under the nyquist frequency (' ...
            num2str(fnyq) ')!']);
    end
    
    % setup output
    [svol(1:nrng,1).nsta]=deal(nsta);
    [svol(1:nrng,1).stla]=deal(loc(:,1));
    [svol(1:nrng,1).stlo]=deal(loc(:,2));
    [svol(1:nrng,1).stel]=deal(loc(:,3));
    [svol(1:nrng,1).stdp]=deal(loc(:,4));
    [svol(1:nrng,1).butc]=deal([0 0 0 0 0]);
    [svol(1:nrng,1).eutc]=deal([0 0 0 0 0]);
    [svol(1:nrng,1).delta]=deal(delta(1));
    [svol(1:nrng,1).npts]=deal(npts(1));
    [svol(1:nrng,1).polar]=deal(polar);
    [svol(1:nrng,1).npairs]=deal(ncorr);
    [svol(1:nrng,1).method]=deal(method);
    [svol(1:nrng,1).center]=deal([clat clon]);
    [svol(1:nrng,1).volume]=deal(true);
    [svol(1:nrng,1).weights]=deal(w);
    
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
        [svol(1:nrng,1).y]=deal(smag'*d2km);
        smag=smag(ones(bazpts,1),:)';
        baz=(0:bazpts-1)/(bazpts-1)*360*d2r;
        [svol(1:nrng,1).x]=deal(baz/d2r);
        baz=baz(ones(spts,1),:);
        p=2*pi*1i*[smag(:).*sin(baz(:)) smag(:).*cos(baz(:))]*r;
        clear r smag baz
    else % cartesian
        spts=spts(1); bazpts=spts;
        sx=-smax:2*smax/(spts-1):smax;
        [svol(1:nrng,1).x]=deal(sx*d2km);
        [svol(1:nrng,1).y]=deal(fliplr(sx*d2km)');
        sx=sx(ones(spts,1),:);
        sy=fliplr(sx)';
        p=2*pi*1i*[sx(:) sy(:)]*r;
        clear r sx sy
    end
    
    % loop over frequency ranges
    for a=1:nrng
        % get frequencies
        fidx=find(f>=frng(a,1) & f<=frng(a,2));
        svol(a).freq=f(fidx);
        nfreq=numel(fidx);
        
        % preallocate fk space
        svol(a).beam=zeros(spts,bazpts,nfreq,'single');
        
        % warning if no frequencies
        if(~nfreq)
            warning('seizmo:fkxcvolume:noFreqs',...
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
            fprintf('Getting fk Volume %d for %g to %g Hz\n',...
                a,frng(a,1),frng(a,2));
            print_time_left(0,nfreq);
        end
        
        % loop over frequencies
        for b=1:nfreq
            % form beam & convert to dB
            % - following Koper, Seats, and Benz 2010 in BSSA
            %   by only using the real component
            svol(a).beam(:,:,b)=reshape(...
                10*log10(abs(real(exp(f(fidx(b))*p)*cs(:,b)))),...
                spts,bazpts);
            
            % detail message
            if(verbose); print_time_left(b,nfreq); end
        end
        
        % normalize so max peak is at 0dB
        svol(a).normdb=max(svol(a).beam(:));
        svol(a).beam=svol(a).beam-svol(a).normdb;
        
        % plot if no output
        if(~nargout); fkfreqslide(svol(a)); drawnow; end
    end
    
    % return struct
    if(nargout); varargout{1}=svol; end
    
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
