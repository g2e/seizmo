function [s]=fss(data,smax,spts,frng,varargin)
%FSS    Estimate frequency-slowness spectrum
%
%    Usage:    s=fss(data,smax,spts,frng)
%              s=fss(data,smax,[nspts espts],frng)
%              s=fss(data,[nsmin nsmax nspts],[esmin esmax espts],frng)
%              s=fss(...,'polar',true|false,...)
%              s=fss(...,'method',string,...)
%              s=fss(...,'whiten',true|false,...)
%              s=fss(...,'weights',w,...)
%              s=fss(...,'damping',d,...)
%              s=fss(...,'ntiles',nt,...)
%              s=fss(...,'fhwidth',n,...)
%              s=fss(...,'prep',true|false,...)
%              s=fss(...,'avg',true|false,...)
%
%    Description:
%     S=FSS(DATA,SMAX,SPTS,FRNG) estimates the frequency-slowness
%     power spectra for an array by frequency-domain beamforming the time
%     series dataset DATA in a cartesian grid of size SPTSxSPTS extending
%     from -SMAX to SMAX (sec/deg) in both North & East horizontal slowness
%     space.  The dataset DATA is a SEIZMO struct containing array info and
%     time series recordings.  FRNG gives the frequency range as
%     [FREQLOW FREQHIGH] in Hz (the individual frequencies are determined
%     by a fft of the data).  FRNG may also be given as a single frequency
%     in which case the nearest discrete frequency from the fft operation
%     is returned.  The output S is a struct containing relevant info and
%     the frequency-slowness spectra (with size SPTSxSPTSxNFREQ).  The
%     struct layout is:
%          .nsta     - number of stations
%          .st       - station positions [lat lon elev depth]
%          .butc     - UTC start time of data
%          .eutc     - UTC end time of data
%          .npts     - number of time points
%          .delta    - number of seconds between each time point
%          .polar    - true if slowness is sampled in polar coordinates
%          .x        - east slowness (sec/deg) or azimuth values (deg)
%          .y        - north or radial slowness values (sec/deg)
%          .freq     - frequency values (Hz)
%          .npairs   - number of pairs
%          .method   - beamforming method ('center', 'coarray', etc)
%          .center   - array center as [LAT LON]
%          .whiten   - is spectra whitened? true/false
%          .weights  - weights used in beamforming
%          .ntiles   - number of timesection tiles
%          .fhwidth  - frequency smoothing halfwidth in samples
%          .damping  - damping for cross spectral matrix inversion
%          .spectra  - frequency-slowness spectra estimate
%
%     S=FSS(DATA,SMAX,[NSPTS ESPTS],FRNG) specifies the horizontal
%     slowness sample size separately.  For example [200 100] would give
%     twice the sampling in the North slowness space compared to that of
%     the East slowness space.
%
%     S=FSS(DATA,[NSMIN NSMAX NSPTS],[ESMIN ESMAX ESPTS],FRNG) gives
%     the range and sampling of the North & East slowness.  This allows
%     targeted sampling of slowness space to save computation and memory.
%     The 3rd element (the number of samples) may be omitted (will default
%     to 101).
%
%     S=FSS(...,'POLAR',TRUE|FALSE,...) specifies if the spectra is
%     sampled regularly in cartesian or polar coordinates.  Polar coords
%     are useful for slicing the spectra by azimuth (pie slice) or slowness
%     (rings).  Cartesian coords (the default) samples the slowness space
%     regularly in the East/West & North/South directions and so exhibits
%     less distortion in plots of the slowness space. If POLAR=TRUE, SPTS
%     may be given as [SPTS BAZPTS] to control the azimuthal resolution
%     (default is BAZPTS=181 points).  BAZPTS goes from 0 to 360 including
%     both end points such that the true resolution in azimuthal space is
%     BAZPTS-1.  Note that SPTS goes from 0 to SMAX in this case rather
%     than from -SMAX to SMAX for cartesian (so reduce SMAX by 2 if you
%     want to keep a similar spacing).  You may also use:
%      [SMIN SMAX SPTS],[BAZMIN BAZMAX BAZPTS]
%     to explicitly specify the sampled region in absolute slowness and
%     back-azimuth.  The 3rd element (the number of samples) may be omitted
%     (SPTS defaults to 101, BAZPTS defaults to 181).
%
%     S=FSS(...,'METHOD',STRING,...) defines the beam method.  STRING
%     may be 'center', 'coarray', 'capon', 'full', or [LAT LON].  The
%     default is 'center' which is extremely fast for large arrays compared
%     to the other methods as it phase delays the auto-correlation of each
%     record using the distance to the array center (using ARRAYCENTER).
%     The 'coarray' method utilizes all unique pairings of records to
%     compute the spectrum with a little less than half the cross power
%     spectral density matrix while the 'full' method uses the entire cross
%     power spectral density matrix.  The 'full' method is significantly
%     slower and gives degraded results compared to the 'coarray' method
%     and so is not recommended except in verification (as it gives the
%     exact same result of the 'center' method).  The 'coarray' method has
%     the best resolution out of those three methods but is inferior to the
%     'capon' method.  The 'capon' method is similar to the 'full' method
%     in computation time (the only additional operation is a dampened
%     inversion of the full cross power spectral density matrix) but the
%     resolution at every slowness is optimized in a least-squares sense
%     using the maximum likelyhood method.  Because of the reliance of the
%     'capon' method on the actual data, the signal to incoherent noise
%     levels must be high otherwise the 'coarray' method is the best
%     choice.  Using [LAT LON] as a method is algorithmically the same as
%     the 'center' method but uses the defined coordinates as the beam
%     location for the array.
%
%     S=FSS(...,'WHITEN',TRUE|FALSE,...) whitens the spectras before
%     beamforming if WHITEN is TRUE.  The default is TRUE.
%
%     S=FSS(...,'WEIGHTS',W,...) gives the relative weights for each
%     record in DATA (must match size of DATA) or pairing (only if METHOD
%     is 'coarray', 'capon', or 'full').
%
%     S=FSS(...,'DAMPING',D,...) alters the dampening parameter used
%     in the inversion of the cross power spectral density matrix for the
%     'capon' method.  This is done by "diagonal loading" which means the
%     dampening value is added the elements along the diagonal of the cross
%     power spectral density matrix.  The default value of 0.001 likely
%     needs to be adjusted for the best results.
%
%     S=FSS(...,'NTILES',NT,...) sets how many nonoverlapping
%     timesections the data should be split into.  This is mainly for the
%     'capon' method as it needs time/frequency averaging of the cross
%     power spectral density matrix for stability.  Setting NT>=N where N
%     is the number of records in DATA is recommended by Capon 1969.  The
%     default is 1 which provides the best frequency resolution.  Use this
%     in combination with FHWIDTH option to balance between time,
%     frequency, and slowness resolution.
%
%     S=FSS(...,'FHWIDTH',N,...) sets the sliding window halfwidth used to
%     average neighboring frequencies of the cross spectral matrix.  This
%     is mainly for stabilizing the inversion in the 'capon' method but has
%     utility in smoothing spectra.  The window is size 2N+1 so for example
%     N=2 averages each frequency with the closest 2 discrete frequencies
%     above and below (a 5 point sliding average).  The default N=0 which
%     does not average.
%
%     S=FSS(...,'PREP',TRUE|FALSE,...) prepares the tiled data by
%     detrending and tapering.  The default is FALSE.
%
%     S=FSS(...,'AVG',TRUE|FALSE,...) indicates if the slowness spectra is
%     averaged across frequency during computation (thus collapsing the
%     frequency dimension to size 1).  This is not related to the FHWIDTH
%     option (and therefore does not help to stabilize the 'capon' method)
%     but does save a significant amount of memory.  The default is false.
%
%    Notes:
%     - Records in DATA must have equal and regular sample spacing.
%     - Attenuation across the array is not included.
%     - dB values for unwhitened data are probably not scaled correctly but
%       this does not have an impact on relative dB values.
%     - References:
%        Capon 1969, High-Resolution Frequency-Wavenumber Spectrum
%         Analysis, Proc. IEEE, Vol. 57, No. 8, pp. 1408-1418
%        Husebye & Ruud 1989, Array Seismology - Past, Present and Future
%         Developments, in Observatory Seismology, edited by Litehiser, pp.
%         123-153, Univ of Calif Press, Berkeley
%        Rost & Thomas 2002, Array Seismology: Methods and Applications,
%         Rev of Geoph, Vol. 40, No. 3, doi:10.1029/2000RG000100
%
%    Examples:
%     % Show slowness spectra for an artificial dataset at 40s periods:
%     plotfss(fss(capon1970,50,201,1/40));
%
%     % Now compare that to the Capon Estimation:
%     plotfss(fss(capon1970,50,201,1/40,'m','capon'));
%
%    See also: FSSXC, ARF, SNYQUIST, PLOTFSS, KXY2SLOWBAZ, SLOWBAZ2KXY,
%              FSSAVG, FSSSUB, FSSHORZ, FSSHORZXC, ARFHORZ, FSSDBINFO,
%              FSSFREQSLIDE, FSSFRAMESLIDE, PLOTARF, FSSCORRCOEF

%     Version History:
%        May   3, 2010 - initial version
%        May   7, 2010 - only doing one triangle gives better beam and
%                        takes less than half the time
%        May   8, 2010 - array math version (another big speed jump), added
%                        a couple options for doing rose vs grid slowness
%                        plots, also allow changing between coarray & a
%                        specified array center (specified is hugely faster
%                        but suffers in resolution slightly)
%        May   9, 2010 - struct changes
%        May  10, 2010 - use checkheader more effectively
%        May  12, 2010 - fixed an annoying bug (2 wrongs giving the right
%                        answer), minor code touches while debugging
%                        fkxcvolume
%        May  24, 2010 - beam is now single precision
%        June 16, 2010 - fixed nargchk, better verbose msg, improved see
%                        also section, create beam as s.p. array
%        July  1, 2010 - high latitude fix, allocation fix
%        July  6, 2010 - major update to struct, doc update
%        Nov. 18, 2010 - added weighting
%        Apr.  3, 2012 - use seizmocheck
%        May  30, 2012 - pow2pad=0 by default
%        June 13, 2012 - add capon method (not working!)
%        Sep. 11, 2012 - altered from geofss & fkmap
%        Sep. 21, 2012 - allow 0-Fnyq range (for full spectrum)
%        Sep. 22, 2012 - allow station or pair weights
%        Sep. 27, 2012 - pv pair inputs, capon method works with damping,
%                        doc update, less memory usage, error for no freq
%        Sep. 28, 2012 - tiling works for all methods
%        Sep. 30, 2012 - avg option
%        Oct. 10, 2012 - doc update, lower mem usage (expand dt on the
%                        fly), whitening for capon method
%        Jan.  9, 2013 - bugfix for norm. when record is all zeros, allow
%                        options to be any case
%        Jan. 14, 2013 - update history
%        Jan. 15, 2014 - allow far more flexible slowness space sampling,
%                        fhwidth option, tiles are averaged
%        Mar. 24, 2014 - added prep option, improved slowness space
%                        sampling some more, fixed implementation of ntiles
%                        and fhwidth options
%        Mar. 25, 2014 - single frequency option
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2014 at 14:05 GMT

% todo:

% check nargin
error(nargchk(4,inf,nargin));

% check struct
error(seizmocheck(data,'dep'));

% require not a xc dataset
if(isxc(data))
    error('seizmo:fss:badInput',...
        'Use FSSXC for correlations!');
end

% number of records
nrecs=numel(data);

% need 2+ records
if(nrecs<2)
    error('seizmo:fss:arrayTooSmall',...
        'DATA must have 2+ records!');
end

% check inputs
if(isscalar(smax))
    % (smax,spts) or (smax,[ypts xpts])
    %       ep=np or (smax,[npts epts])
    %      bp=181 or (smax,[spts bazpts])
    if(~isnumeric(smax) || ~isreal(smax) || smax<=0)
        error('seizmo:fss:badInput',...
            'SMAX must be a positive real scalar in sec/deg!');
    elseif(~isnumeric(spts) || ~isreal(spts) || ~isvector(spts) ...
            || isempty(spts) || any(spts<2) || any(fix(spts)~=spts))
        error('seizmo:fss:badInput',...
            'SPTS must be positive integer(s) >1!');
    elseif(numel(spts)>2)
        error('seizmo:fss:badInput',...
            'SPTS must contain 1 or 2 elements!');
    end
else
    % ([ymin ymax],[xmin xmax]) or ([ymin ymax ypts],[xmin xmax xpts])
    % ([nmin nmax],[emin emax]) with npts=101
    % ([smin smax],[bazmin bazmax]) with spts=101, bazpts=181
    if(~isnumeric(smax) || ~isreal(smax))
        error('seizmo:fss:badInput',...
            '[YMIN YMAX YPTS] input must be real-valued!');
    elseif(~isnumeric(spts) || ~isreal(spts))
        error('seizmo:fss:badInput',...
            '[XMIN XMAX XPTS] input must be real-valued!');
    end
    if(isempty(smax) || numel(smax)>3)
        error('seizmo:fss:badInput',...
            '2nd argument must be [YMIN YMAX YPTS] or [YMIN YMAX]!');
    elseif(isempty(smax) || numel(smax)>3)
        error('seizmo:fss:badInput',...
            '3rd argument must be [XMIN XMAX XPTS] or [XMIN XMAX]!');
    end
    if(numel(smax)==3 && (smax(3)<1 || smax(3)~=fix(smax(3))))
        error('seizmo:fss:badInput',...
            'YPTS must be a positive integer!');
    elseif(numel(spts)==3 && (spts(3)<1 || spts(3)~=fix(spts(3))))
        error('seizmo:fss:badInput',...
            'XPTS must be a positive integer!');
    end
end
sf=size(frng);
if(~isreal(frng) || numel(sf)~=2 || any(frng(:)<0))
    error('seizmo:fss:badInput',...
        'FRNG must be 1 or more positive real values in Hz!');
end
if(sf(2)==1); frng=frng(:,[1 1]); sf(2)=2; end
if(sf(2)~=2 || any(frng(:,1)>frng(:,2)))
    error('seizmo:fss:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
end
nrng=sf(1);

% parse options
pv=parse_fss_pv_pairs(varargin{:});

% defaults for optionals
if(isempty(pv.w)); pv.w=ones(nrecs,1); end

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

% extract header info & data
try
    % verbosity
    verbose=seizmoverbose;
    
    % grap necessary header info
    [npts,delta,butc,eutc,st]=getheader(data,...
        'npts','delta','b utc','e utc','st');
    butc=cell2mat(butc); eutc=cell2mat(eutc);
    
    % extract data (silently)
    seizmoverbose(false);
    data=records2mat(data);
    seizmoverbose(verbose);
    
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

% fix method/center/npairs
if(ischar(pv.method))
    pv.method=lower(pv.method);
    [clat,clon]=arraycenter(st(:,1),st(:,2));
    switch pv.method
        case 'coarray'
            npairs=nrecs*(nrecs-1)/2;
        case {'full' 'capon'}
            npairs=nrecs*nrecs;
        case 'center'
            npairs=nrecs;
    end
else
    clat=pv.method(1);
    clon=pv.method(2);
    pv.method='user';
    npairs=nrecs;
end

% check weights again
if(~any(numel(pv.w)==[nrecs npairs]))
    error('seizmo:fss:badInput',...
        'Number of WEIGHTS must match the number of stations or pairs!');
end

% get time limits
[bmin,bmini]=min(timediff(butc(1,:),butc,'utc'));
[emax,emaxi]=max(timediff(eutc(1,:),eutc,'utc'));

% check nyquist
fnyq=1/(2*delta(1));
if(any(frng(:,1)>fnyq))
    error('seizmo:fss:badFRNG',...
        ['FRNG exceeds nyquist frequency (' num2str(fnyq) ')!']);
end

% longest record gives record length for fft
maxnpts=max(npts);

% points in tiles based on longest record
if(pv.ntiles>1)
    oldmax=maxnpts;
    maxnpts=ceil(maxnpts/pv.ntiles);
end

% get frequencies
nspts=2^nextpow2(maxnpts); % half xcorr
%nspts=2^nextpow2(2*maxnpts-1); % full xcorr for verification purposes
f=(0:nspts/2)/(delta(1)*nspts);  % only +freq
nf=numel(f);

% tiling records
if(pv.ntiles>1)
    data=[data; zeros(maxnpts*pv.ntiles-oldmax,nrecs)];
    data=permute(reshape(data.',[nrecs maxnpts pv.ntiles]),[2 1 3]);
end

% remove trend & taper
if(pv.prep)
    data(:,:)=detrend(data(:,:),'constant');
    data(:,:)=detrend(data(:,:));
    data(:,:)=repmat(hann(maxnpts),[1 size(data(:,:),2)]).*data(:,:);
end

% get fft
data=fft(data,nspts,1);

% trim -freq and permute
if(pv.ntiles>1)
    data=permute(data(1+(0:nspts/2),:,:),[2 1 3]);
else
    data=data(1+(0:nspts/2),:).';
end

% whiten data if desired
if(pv.whiten); data=data./abs(data); data(isnan(data))=0; end

% setup slowness grid
if(pv.polar)
    if(numel(smax)==1)
        if(numel(spts)==1); spts(2)=181; end % default # azimuthal points
        sx=(0:spts(2)-1)/(spts(2)-1)*360; % baz (wedge decided x/y)
        sy=(0:spts(1)-1).'/(spts(1)-1)*smax; % smag
    else
        if(numel(smax)==2); smax(3)=101; end % default # slow mag points
        if(numel(spts)==2); spts(3)=181; end % default # azimuthal points
        sx=linspace(spts(1),spts(2),spts(3)); % baz (wedge decided x/y)
        sy=linspace(smax(1),smax(2),smax(3)).'; % smag
        spts=[smax(3) spts(3)]; % save npts for later
    end
else
    if(numel(smax)==1)
        if(numel(spts)==1); spts(2)=spts(1); end % default # east points
        sx=-smax:2*smax/(spts(2)-1):smax; % east
        sy=fliplr(-smax:2*smax/(spts(1)-1):smax).'; % north
    else
        if(numel(smax)==2); smax(3)=101; end % default # north points
        if(numel(spts)==2); spts(3)=101; end % default # east points
        sx=linspace(spts(1),spts(2),spts(3)); % east
        sy=fliplr(linspace(smax(1),smax(2),smax(3))).'; % north
        spts=[smax(3) spts(3)]; % save npts for later
    end
end

% setup output
[s(1:nrng,1).nsta]=deal(nrecs);
[s(1:nrng,1).st]=deal(st);
[s(1:nrng,1).butc]=deal(butc(bmini,:));
[s(1:nrng,1).eutc]=deal(eutc(emaxi,:));
[s(1:nrng,1).delta]=deal(delta(1));
[s(1:nrng,1).npts]=deal(maxnpts);
[s(1:nrng,1).polar]=deal(pv.polar);
[s(1:nrng,1).x]=deal(sx);
[s(1:nrng,1).y]=deal(sy);
[s(1:nrng,1).freq]=deal([]);
[s(1:nrng,1).method]=deal(pv.method);
[s(1:nrng,1).npairs]=deal(npairs);
[s(1:nrng,1).center]=deal([clat clon]);
[s(1:nrng,1).whiten]=deal(pv.whiten);
[s(1:nrng,1).weights]=deal(pv.w);
[s(1:nrng,1).ntiles]=deal(pv.ntiles);
[s(1:nrng,1).fhwidth]=deal(pv.fhwidth);
[s(1:nrng,1).damping]=deal(pv.damping);
[s(1:nrng,1).spectra]=deal(zeros(0,'single'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of defaults, checks, & data setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% radius constants
d2r=pi/180;
d2km=6371*d2r;

% get number of slowness points
nslow=prod(spts);

% get utc static time shifts
% - case 'center' 'user' only does the line below
dt=timediff(butc(1,:),butc,'utc').'; % need row vector
switch pv.method
    case 'coarray'
        % row is master index, column is slave index
        [master,slave]=find(triu(true(nrecs),1));
        dt=dt(slave)-dt(master); % alter to pair shifts
    case {'full' 'capon'}
        [master,slave]=find(true(nrecs));
        dt=dt(slave)-dt(master); % alter to pair shifts
end

% normalize & expand weights
pv.w=pv.w(:);
switch pv.method
    case {'coarray' 'full' 'capon'}
        if(numel(pv.w)==nrecs); pv.w=pv.w(slave).*conj(pv.w(master)); end
end
pv.w=pv.w./sum(abs(pv.w));

% get relative positions for each pair (the coarray)
% r=(x  ,y  )
%     ij  ij
%
% x   is km east, y   is km north
%  ij              ij
%
% r is a 2xNPAIRS matrix
switch pv.method
    case {'coarray' 'full' 'capon'}
        % [ r   r   ... r
        %    11  12      1N
        %   r   r   ... r
        %    21  22      2N
        %    .   .  .    .
        %    .   .   .   .
        %    .   .    .  .
        %   r   r   ... r   ]
        %    N1  N2      NN
        [e,n]=geographic2enu(st(:,1),st(:,2),0,clat,clon,0);
        e=e(slave)-e(master);
        n=n(slave)-n(master);
    otherwise % {'center' 'user'}
        % each station relative to array center
        [e,n]=geographic2enu(st(:,1),st(:,2),0,clat,clon,0);
end
r=[e(:) n(:)]';
clear e n

% Get phasors for each wave speed+direction (slowness) and pair
% which allows us to delay and sum in the frequency domain.
% p=s*r
%
% where r was defined above
%       s is the slowness vector s=(s ,s ) and is NSLOWx2
%                                    x  y
%       and s  is in sec/km east, s  is in sec/km north
%            x                     y
%
% p is the projection of the slowness vectors s onto the
% spatial difference vectors r (called the coarray)
%
% p is NSLOWxNPAIRS and has units of seconds
if(pv.polar)
    sx=sx(ones(spts(1),1),:); % baz in degrees
    sy=sy/d2km; % smag in sec/km
    sy=sy(:,ones(spts(2),1));
    [sx,sy]=deal(sy.*sind(sx),sy.*cosd(sx)); % convert to cartesian
else % cartesian
    sx=sx/d2km;
    sx=sx(ones(spts(1),1),:);
    sy=sy/d2km;
    sy=sy(:,ones(spts(2),1));
end
p=[sx(:) sy(:)]*r+dt(ones(nslow,1),:);
clear r sx sy

% loop over frequency ranges
for a=1:nrng
    % get frequencies
    if(frng(a,1)==frng(a,2))
        % nearest single frequency
        [fidx,fidx]=min(abs(f-frng(a,1)));
    else
        fidx=find(f>=frng(a,1) & f<=frng(a,2));
    end
    s(a).freq=f(fidx);
    nfreq=numel(fidx);
    
    % preallocate spectra
    if(pv.avg); s(a).spectra=zeros(spts,'single');
    else s(a).spectra=zeros([spts nfreq],'single');
    end
    
    % error if no frequencies
    if(~nfreq)
        error('seizmo:fss:noFreqs',...
            'No frequencies within the range %g to %g Hz!',...
            frng(a,1),frng(a,2));
    end
    
    % detail message
    if(verbose)
        fprintf('Getting spectra %d: %g to %g Hz\n',...
            a,frng(a,1),frng(a,2));
        print_time_left(0,nfreq);
    end
    
    % proceed by type
    switch pv.method
        case {'coarray' 'full'}
            if(pv.avg)
                for b=1:nfreq
                    % get frequency range
                    newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                    newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                    
                    % compute average cross spectra for freq range & tiles
                    % then weighted beam and add
                    s(a).spectra=s(a).spectra+reshape(real(...
                        exp(-2*pi*1i*f(fidx(b))*p)...
                        *(mean(mean(data(slave,newfidx,:)...
                        .*conj(data(master,newfidx,:)),2),3).*pv.w)),spts);
                    if(verbose); print_time_left(b,nfreq); end
                end
            else
                for b=1:nfreq
                    % get frequency range
                    newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                    newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                    
                    % compute average cross spectra for freq range & tiles
                    % then weighted beam and insert
                    s(a).spectra(:,:,b)=reshape(real(...
                        exp(-2*pi*1i*f(fidx(b))*p)...
                        *(mean(mean(data(slave,newfidx,:)...
                        .*conj(data(master,newfidx,:)),2),3).*pv.w)),spts);
                    if(verbose); print_time_left(b,nfreq); end
                end
            end
        case 'capon'
            if(pv.avg)
                for b=1:nfreq
                    % get frequency range
                    newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                    newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                    
                    % get average cross spectral matrix
                    cs=reshape(mean(mean(data(slave,newfidx,:).*conj(...
                        data(master,newfidx,:)),2),3),[nrecs nrecs]);
                    
                    % damped inversion of cross spectral matrix
                    cs=pinv((1-pv.damping)*cs+pv.damping*eye(nrecs));
                    
                    % weighted capon beam and add
                    s(a).spectra=s(a).spectra+1./reshape(real(exp(...
                        -2*pi*1i*f(fidx(b))*p)*(cs(:).*pv.w)),spts);
                    if(verbose); print_time_left(b,nfreq); end
                end
            else
                for b=1:nfreq
                    % get frequency range
                    newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                    newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                    
                    % get average cross spectral matrix
                    cs=reshape(mean(mean(data(slave,newfidx,:).*conj(...
                        data(master,newfidx,:)),2),3),[nrecs nrecs]);
                    
                    % damped inversion of cross spectral matrix
                    cs=pinv((1-pv.damping)*cs+pv.damping*eye(nrecs));
                    
                    % weighted capon beam and insert
                    s(a).spectra(:,:,b)=1./reshape(real(exp(...
                        -2*pi*1i*f(fidx(b))*p)*(cs(:).*pv.w)),spts);
                    if(verbose); print_time_left(b,nfreq); end
                end
            end
        otherwise % {'center' 'user'}
            if(pv.avg)
                for b=1:nfreq
                    % get frequency range
                    newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                    newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                    newnf=numel(newfidx);
                    
                    for c=1:pv.ntiles
                        s(a).spectra=s(a).spectra+reshape(mean(abs(exp(...
                            -2*pi*1i*f(fidx(b))*p)*(data(:,newfidx,c)...
                            .*pv.w(:,ones(1,newnf)))).^2,2),spts);
                    end
                    if(verbose); print_time_left(b,nfreq); end
                end
            else
                for b=1:nfreq
                    % get frequency range
                    newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                    newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                    newnf=numel(newfidx);
                    
                    for c=1:pv.ntiles
                        s(a).spectra(:,:,b)=s(a).spectra(:,:,b)+reshape(...
                            mean(abs(exp(-2*pi*1i*f(fidx(b))*p)*(data(:,...
                            newfidx,c).*pv.w(:,ones(1,newnf)))).^2,2),...
                            spts);
                    end
                    if(verbose); print_time_left(b,nfreq); end
                end
            end
            s(a).spectra=s(a).spectra/pv.ntiles;
    end
    if(pv.avg); s(a).spectra=s(a).spectra/nfreq; end
end

end


function [pv]=parse_fss_pv_pairs(varargin)

% defaults
pv.avg=false;
pv.polar=false;
pv.method='center';
pv.whiten=true;
pv.w=[];
pv.damping=0.001; % only for capon
pv.ntiles=1;
pv.fhwidth=0;
pv.prep=false;

% require pv pairs
if(mod(nargin,2))
    error('seizmo:fss:badInput',...
        'Unpaired parameter/value!');
end

% parse pv pairs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:fss:badInput',...
        'Parameters must be specified using a string!');
end
for i=1:2:nargin
    switch lower(varargin{i})
        case {'polar' 'plr' 'pol' 'p'}
            pv.polar=varargin{i+1};
        case {'method' 'meth' 'm'}
            pv.method=varargin{i+1};
        case {'whiten' 'wh' 'white' 'normalize' 'n' 'coherency' 'c'}
            pv.whiten=varargin{i+1};
        case {'weights' 'w' 'wgt' 'wgts' 'weight'}
            pv.w=varargin{i+1};
        case {'damping' 'damp' 'd'}
            pv.damping=varargin{i+1};
        case {'ntiles' 'ntile' 'tiles' 'tile' 'nt' 't'}
            pv.ntiles=varargin{i+1};
        case {'fhwidth' 'fhwid' 'fh' 'fwid' 'f'}
            pv.fhwidth=varargin{i+1};
        case {'average' 'av' 'avg' 'a'}
            pv.avg=varargin{i+1};
        case {'prep'}
            pv.prep=varargin{i+1};
        otherwise
            error('seizmo:fss:badInput',...
                'Unknown parameter: %s !',varargin{i});
    end
end

% valid method strings
valid.METHOD={'center' 'coarray' 'full' 'capon'};

% check values
if(~isscalar(pv.polar) || (~islogical(pv.polar) && ~isnumeric(pv.polar)))
    error('seizmo:fss:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif((isnumeric(pv.method) ...
        && (~isreal(pv.method) || ~numel(pv.method)==2)) ...
        || (ischar(pv.method) && ~any(strcmpi(pv.method,valid.METHOD))))
    error('seizmo:fss:badInput',...
        ['METHOD must be one of the following:\n' ...
        '''CENTER'', ''COARRAY'', ''CAPON'', ''FULL'', or [LAT LON]!']);
elseif(~isscalar(pv.whiten) || ~islogical(pv.whiten))
    error('seizmo:fss:badInput',...
        'WHITEN must be TRUE or FALSE!');
elseif(~isnumeric(pv.w) || any(pv.w(:)<0))
    error('seizmo:fss:badInput',...
        'WEIGHTS must be positive!');
elseif(~isreal(pv.damping) || ~isscalar(pv.damping) || pv.damping<0)
    error('seizmo:fss:badInput',...
        'DAMPING must be a positive scalar!');
elseif(~isreal(pv.ntiles) || ~isscalar(pv.ntiles) ...
        || pv.ntiles<=0 || pv.ntiles~=fix(pv.ntiles))
    error('seizmo:fss:badInput',...
        'NTILES must be a positive integer!');
elseif(~isreal(pv.fhwidth) || ~isscalar(pv.fhwidth) ...
        || pv.fhwidth~=fix(pv.fhwidth))
    error('seizmo:fss:badInput',...
        'FHWIDTH must be an integer!');
elseif(~isscalar(pv.avg) || ~islogical(pv.avg))
    error('seizmo:fss:badInput',...
        'AVG must be TRUE or FALSE!');
elseif(~isscalar(pv.prep) || ~islogical(pv.prep))
    error('seizmo:fss:badInput',...
        'PREP must be TRUE or FALSE!');
end

end


function [lgc]=isxc(data)
[m,s]=getheader(data,'kuser0','kuser1');
if(all(strcmp(m,'MASTER') & strcmp(s,'SLAVE'))); lgc=true;
else lgc=false;
end
end
