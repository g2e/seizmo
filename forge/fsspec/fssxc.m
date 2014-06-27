function [s]=fssxc(xcdata,smax,spts,frng,varargin)
%FSSXC    Estimate frequency-slowness spectrum with cross-correlations
%
%    Usage:    s=fssxc(xcdata,smax,spts,frng)
%              s=fssxc(xcdata,smax,[nspts espts],frng)
%              s=fssxc(xcdata,[nsmin nsmax nspts],[esmin esmax espts],frng)
%              s=fssxc(...,'polar',true|false,...)
%              s=fssxc(...,'method',string,...)
%              s=fssxc(...,'whiten',true|false,...)
%              s=fssxc(...,'weights',w,...)
%              s=fssxc(...,'damping',d,...)
%              s=fssxc(...,'fhwidth',n,...)
%              s=fssxc(...,'avg',true|false,...)
%
%    Description:
%     S=FSSXC(XCDATA,SMAX,SPTS,FRNG) estimates the frequency-slowness
%     power spectra for an array by frequency-domain beamforming the
%     correlation dataset XCDATA in a cartesian grid of size SPTSxSPTS
%     extending from -SMAX to SMAX (sec/deg) in both North & East
%     horizontal slowness space.  The dataset XCDATA is a SEIZMO struct
%     containing array info and correlations.  FRNG gives the frequency
%     range as [FREQLOW FREQHIGH] in Hz (the individual frequencies are
%     determined by the fft of the data).  FRNG may also be given as a
%     single frequency in which case the nearest discrete frequency from
%     the fft operation is returned.  The output S is a struct containing
%     relevant info and the frequency-slowness spectra (with size
%     SPTSxSPTSxNFREQ).  The struct layout is:
%          .nsta     - number of stations (uses naming info to get this)
%          .st       - station positions [lat lon elev depth]
%          .butc     - UTC start time of data
%          .eutc     - UTC end time of data
%          .npts     - number of time points
%          .delta    - number of seconds between each time point
%          .polar    - true if slowness is sampled in polar coordinates
%          .x        - east slowness (sec/deg) or azimuth values (deg)
%          .y        - north or radial slowness values (sec/deg)
%          .freq     - frequency values (Hz)
%          .npairs   - number of pairs (aka correlations)
%          .method   - beamforming method ('xc' or 'caponxc')
%          .center   - array center as [LAT LON]
%          .whiten   - is spectra whitened? true/false
%          .weights  - weights used in beamforming
%          .ntiles   - number of timesection tiles
%          .fhwidth  - frequency smoothing halfwidth in samples
%          .damping  - damping for cross spectral matrix inversion
%          .spectra  - frequency-slowness spectra estimate
%
%     S=FSSXC(XCDATA,SMAX,[NSPTS ESPTS],FRNG) specifies the horizontal
%     slowness sample size separately.  For example [200 100] would give
%     twice the sampling in the North slowness space compared to that of
%     the East slowness space.
%
%     S=FSSXC(XCDATA,[NSMIN NSMAX NSPTS],[ESMIN ESMAX ESPTS],FRNG) gives
%     the range and sampling of the North & East slowness.  This allows
%     targeted sampling of slowness space to save computation and memory.
%     The 3rd element (the number of samples) may be omitted (will default
%     to 101).
%
%     S=FSSXC(...,'POLAR',TRUE|FALSE,...) specifies if the spectra is
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
%     S=FSSXC(...,'METHOD',STRING,...) specifies the beamforming method.
%     STRING may be either 'xc' or 'caponxc'.  The default is 'xc' which
%     uses the correlations as a substitute for the cross power spectral
%     density matrix in a conventional frequency-slowness estimation (eg.
%     the 'coarray' or 'full' methods of FSS).  The 'caponxc' method
%     requires that the correlation dataset can form the entire cross power
%     spectral density matrix.  This method can provide a higher resolution
%     estimate of the frequency-slowness spectra.
%
%     S=FSSXC(...,'WHITEN',TRUE|FALSE,...) whitens the cross power spectral
%     matrix elements before beamforming if WHITEN is TRUE.  The default is
%     TRUE.
%
%     S=FSSXC(...,'WEIGHTS',W,...) specifies the relative weights for each
%     correlogram in XCDATA (must match size of XCDATA).
%
%     S=FSSXC(...,'DAMPING',D,...) alters the dampening parameter used in
%     the inversion of the cross power spectral density matrix for the
%     'caponxc' method.  This is done by diagonal loading which means the
%     dampening value is added the elements along the diagonal.  The
%     default value of 0.001 likely needs to be adjusted for best results.
%
%     S=FSSXC(...,'FHWIDTH',N,...) sets the sliding window halfwidth for
%     averaging neighboring frequencies of the cross spectral matrix.  This
%     is mainly for stabilizing the inversion in the 'caponxc' method but
%     has utility in smoothing spectra.  The window is size 2N+1 so for
%     example N=2 averages each frequency with the closest 2 discrete
%     frequencies above and below (a 5 point sliding average).  The default
%     N=0 which does not average.
%
%     S=FSSXC(...,'AVG',TRUE|FALSE,...) indicates if the spectra is
%     averaged across frequency during computation.  This can save a
%     significant amount of memory.  The default is false.
%
%    Notes:
%     - Correlations in XCDATA must have equal and regular sample spacing.
%     - Correlations are expected to have station naming and location info
%       as stored in the header by the function CORRELATE.  This means the
%       headers of correlations done outside of SEIZMO will probably
%       require adjustment.
%
%    Examples:
%     % Show slowness spectra for an artificial dataset at 40s periods:
%     plotfss(fssxc(correlate(capon1970,'mcxc'),50,201,1/40));
%
%     % A full set of correlations is required for a Capon Estimator:
%     plotfss(fssxc(correlate(capon1970,capon1970,'mcxc'),50,201,1/40,...
%         'm','caponxc'));
%
%    See also: FSS, ARF, SNYQUIST, PLOTFSS, KXY2SLOWBAZ, SLOWBAZ2KXY,
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
%        Sep. 17, 2012 - altered from fss, fkxcvolume, & geofssxc
%        Sep. 20, 2012 - no need for splitpad or resampling beforehand
%        Sep. 21, 2012 - allow 0-Fnyq range (for full spectrum)
%        Sep. 22, 2012 - allow station or pair weights
%        Sep. 27, 2012 - pv pair inputs, capon method works with damping,
%                        doc update, error for no freq
%        Sep. 30, 2012 - avg option
%        Jan.  9, 2013 - allow options to be any case, notes on amp issue
%        Jan. 14, 2013 - update history
%        Mar. 31, 2014 - fhwidth option, single frequency selection,
%                        additional gridding flexibility
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 31, 2014 at 14:05 GMT

% todo:
% - discrepancy for full & partial capon input
%   - weights are part of the issue
%   - something else is also a problem
% - crazyiness for fhwidth>0
%
% - discrepancy between fss capon & fssxc caponxc
%   - the absolute dBs are off from each other and from other methods
%     - they vary as a function of frequency which means that addition is
%       broken between the two b/c the contribution varies per frequency
%   - relative dBs are very close but do vary sometimes off
%     - plotfss(fssxc(correlate(capon1970,capon1970,'mcxc'),50,201,1/40,'m','caponxc'),[-48 0]);
%     - plotfss(fss(capon1970,50,201,1/40,'m','capon'),[-48 0]);

% check nargin
error(nargchk(4,inf,nargin));

% check struct
error(seizmocheck(xcdata,'dep'));

% require xc dataset
if(~isxc(xcdata))
    error('seizmo:fssxc:badInput',...
        'XCDATA must be correlations with metadata as from CORRELATE!');
end

% number of pairings
npairs=numel(xcdata);

% check inputs
if(isscalar(smax))
    % (smax,spts) or (smax,[ypts xpts])
    %       ep=np or (smax,[npts epts])
    %      bp=181 or (smax,[spts bazpts])
    if(~isnumeric(smax) || ~isreal(smax) || smax<=0)
        error('seizmo:fssxc:badInput',...
            'SMAX must be a positive real scalar in sec/deg!');
    elseif(~isnumeric(spts) || ~isreal(spts) || ~isvector(spts) ...
            || isempty(spts) || any(spts<2) || any(fix(spts)~=spts))
        error('seizmo:fssxc:badInput',...
            'SPTS must be positive integer(s) >1!');
    elseif(numel(spts)>2)
        error('seizmo:fssxc:badInput',...
            'SPTS must contain 1 or 2 elements!');
    end
else
    % ([ymin ymax],[xmin xmax]) or ([ymin ymax ypts],[xmin xmax xpts])
    % ([nmin nmax],[emin emax]) with npts=101
    % ([smin smax],[bazmin bazmax]) with spts=101, bazpts=181
    if(~isnumeric(smax) || ~isreal(smax))
        error('seizmo:fssxc:badInput',...
            '[YMIN YMAX YPTS] input must be real-valued!');
    elseif(~isnumeric(spts) || ~isreal(spts))
        error('seizmo:fssxc:badInput',...
            '[XMIN XMAX XPTS] input must be real-valued!');
    end
    if(isempty(smax) || numel(smax)>3)
        error('seizmo:fssxc:badInput',...
            '2nd argument must be [YMIN YMAX YPTS] or [YMIN YMAX]!');
    elseif(isempty(smax) || numel(smax)>3)
        error('seizmo:fssxc:badInput',...
            '3rd argument must be [XMIN XMAX XPTS] or [XMIN XMAX]!');
    end
    if(numel(smax)==3 && (smax(3)<1 || smax(3)~=fix(smax(3))))
        error('seizmo:fssxc:badInput',...
            'YPTS must be a positive integer!');
    elseif(numel(spts)==3 && (spts(3)<1 || spts(3)~=fix(spts(3))))
        error('seizmo:fssxc:badInput',...
            'XPTS must be a positive integer!');
    end
end
sf=size(frng);
if(~isreal(frng) || numel(sf)~=2 || any(frng(:)<0))
    error('seizmo:fssxc:badInput',...
        'FRNG must be 1 or more positive real values in Hz!');
end
if(sf(2)==1); frng=frng(:,[1 1]); sf(2)=2; end
if(sf(2)~=2 || any(frng(:,1)>frng(:,2)))
    error('seizmo:fssxc:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
end
nrng=sf(1);

% parse options
pv=parse_fssxc_pv_pairs(varargin{:});

% defaults for optionals
if(isempty(pv.w)); pv.w=ones(npairs,1); end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    xcdata=checkheader(xcdata,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
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

% expand matrix if caponxc
if(strcmpi(pv.method,'caponxc'))
    [lgc,rev]=is_full_matrix_of_correlations(xcdata);
    if(~lgc)
        error('seizmo:fssxc:badInput',...
            'XCDATA does not have all pairs necessary for CAPONXC!');
    end
    
    % add reversed xc if necessary
    if(~isempty(rev))
        xcdata=[xcdata(:); reverse_correlations(xcdata(rev))];
        
        % expand weights if necessary
        if(numel(pv.w)==npairs)
            %pv.w(rev)=pv.w(rev)/2;
            pv.w=[pv.w(:); pv.w(rev)];
        end
        
        % update npairs
        npairs=npairs+numel(rev);
    end
end

% extract header info & xcdata
try
    % verbosity
    verbose=seizmoverbose;
    
    % grab necessary header info
    [b,npts,delta,autc,futc,st,ev,stnm,evnm,ntiles]=getheader(xcdata,...
        'b','npts','delta','a utc','f utc','st','ev','kname','kt','scale');
    autc=cell2mat(autc); futc=cell2mat(futc);
    
    % extract xcdata (silently)
    seizmoverbose(false);
    xcdata=records2mat(xcdata);
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

% get time limits from correlations
autc=autc(all(~isnan(autc),2),:);
futc=futc(all(~isnan(futc),2),:);
if(~isempty(autc))
    [ai,ai]=min(timediff(autc(1,:),autc,'utc'));
    autc=autc(ai,:);
else
    autc=[0 0 0 0 0];
end
if(~isempty(futc))
    [fi,fi]=max(timediff(futc(1,:),futc,'utc'));
    futc=futc(fi,:);
else
    futc=[0 0 0 0 0];
end

% check nyquist
fnyq=1/(2*delta(1));
if(any(frng(:,1)>fnyq))
    error('seizmo:fssxc:badFRNG',...
        ['FRNG exceeds nyquist frequency (' num2str(fnyq) ')!']);
end

% longest record
maxnpts=max(npts);

% get frequencies
nspts=2^nextpow2(maxnpts);
f=(0:nspts/2)/(delta(1)*nspts);  % only +freq
nf=numel(f);

% get fft
xcdata=fft(xcdata,nspts,1);

% trim -freq and permute
xcdata=xcdata(1+(0:nspts/2),:).';

% whiten data if desired
if(pv.whiten); xcdata=xcdata./abs(xcdata); xcdata(isnan(xcdata))=0; end

% use unique stations names to get number of stations & locations
stnm=strcat(stnm(:,1),'.',stnm(:,2),'.',stnm(:,3),'.',stnm(:,4));
evnm=strcat(evnm(:,1),'.',evnm(:,2),'.',evnm(:,3),'.',evnm(:,4));
[stnm,idx1,idx2]=unique([stnm;evnm]);
st=[st;ev];
st=st(idx1,:);
nrecs=numel(stnm);

% check weights again
if(~any(numel(pv.w)==[nrecs npairs]))
    error('seizmo:fssxc:badInput',...
        '# of WEIGHTS must match the # of stations or records in XCDATA!');
end

% slave/master indices
idx2=reshape(idx2,[],2);
slave=idx2(:,1);
master=idx2(:,2);
csidx=sub2ind([nrecs nrecs],master,slave);

% array center
[clat,clon]=arraycenter(st(:,1),st(:,2));

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
[s(1:nrng,1).butc]=deal(autc);
[s(1:nrng,1).eutc]=deal(futc);
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
[s(1:nrng,1).ntiles]=deal(ntiles); % NPAIRSx1
[s(1:nrng,1).fhwidth]=deal(pv.fhwidth);
[s(1:nrng,1).damping]=deal(pv.damping);
[s(1:nrng,1).spectra]=deal(zeros(0,'single'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of defaults, checks, & xcdata setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% radius constants
d2r=pi/180;
d2km=6371*d2r;

% get number of slowness points
nslow=prod(spts);

% get static time shifts
% - MUST shift to zero lag
dt=b.';

% normalize & expand weights
pv.w=pv.w(:);
if(numel(pv.w)==nrecs); pv.w=pv.w(slave).*conj(pv.w(master)); end
pv.w=pv.w./sum(abs(pv.w));

% get relative positions for each pair (the coarray)
% r=(x  ,y  )
%     ij  ij
%
% x   is km east, y   is km north
%  ij              ij
%
% r is a 2xNPAIRS matrix
[e,n]=geographic2enu(st(:,1),st(:,2),0,clat,clon,0);
e=e(slave)-e(master);
n=n(slave)-n(master);
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
% p is NSLOWxNPAIRS
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
    
    % warning if no frequencies
    if(~nfreq)
        error('seizmo:fssxc:noFreqs',...
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
        case 'xc'
            if(pv.avg)
                for b=1:nfreq
                    % get frequency range
                    newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                    newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                    
                    % compute average cross spectra for freq range
                    % then weighted beam and add
                    s(a).spectra=s(a).spectra+reshape(real(...
                        exp(-2*pi*1i*f(fidx(b))*p)...
                        *(mean(xcdata(:,newfidx),2).*pv.w)),spts);
                    if(verbose); print_time_left(b,nfreq); end
                end
            else
                for b=1:nfreq
                    % get frequency range
                    newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                    newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                    
                    % compute average cross spectra for freq range
                    % then weighted beam and insert
                    s(a).spectra(:,:,b)=reshape(real(...
                        exp(-2*pi*1i*f(fidx(b))*p)...
                        *(mean(xcdata(:,newfidx),2).*pv.w)),spts);
                    if(verbose); print_time_left(b,nfreq); end
                end
            end
        otherwise % 'caponxc'
            cs=zeros(nrecs,nrecs);
            if(pv.avg)
                for b=1:nfreq
                    % get frequency range
                    newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                    newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                    
                    % get average cross spectral matrix
                    cs(csidx)=mean(xcdata(:,newfidx),2);
                    
                    % damped inversion of cross spectral matrix
                    cs=pinv((1-pv.damping)*cs+pv.damping*eye(nrecs));
                    
                    % weighted capon beam and add
                    s(a).spectra=s(a).spectra+1./reshape(real(exp(...
                        -2*pi*1i*f(fidx(b))*p)*(cs(csidx).*pv.w)),spts);
                    if(verbose); print_time_left(b,nfreq); end
                end
            else
                for b=1:nfreq
                    % get frequency range
                    newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                    newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                    
                    % get average cross spectral matrix
                    cs(csidx)=mean(xcdata(:,newfidx),2);
                    
                    % damped inversion of cross spectral matrix
                    cs=pinv((1-pv.damping)*cs+pv.damping*eye(nrecs));
                    
                    % weighted capon beam and insert
                    s(a).spectra(:,:,b)=1./reshape(real(exp(...
                        -2*pi*1i*f(fidx(b))*p)*(cs(csidx).*pv.w)),spts);
                    if(verbose); print_time_left(b,nfreq); end
                end
            end
    end
    if(pv.avg); s(a).spectra=s(a).spectra/nfreq; end
end

end


function [pv]=parse_fssxc_pv_pairs(varargin)

% defaults
pv.avg=false;
pv.polar=false;
pv.method='xc';
pv.whiten=true;
pv.w=[];
pv.damping=0.001; % only for caponxc
pv.fhwidth=0;

% require pv pairs
if(mod(nargin,2))
    error('seizmo:fssxc:badInput',...
        'Unpaired parameter/value!');
end

% parse pv pairs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:fssxc:badInput',...
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
        case {'fhwidth' 'fhwid' 'fh' 'fwid' 'f'}
            pv.fhwidth=varargin{i+1};
        case {'average' 'av' 'avg' 'a'}
            pv.avg=varargin{i+1};
        otherwise
            error('seizmo:fssxc:badInput',...
                'Unknown parameter: %s !',varargin{i});
    end
end

% valid method strings
valid.METHOD={'xc' 'caponxc'};

% check values
if(~isscalar(pv.polar) || (~islogical(pv.polar) && ~isnumeric(pv.polar)))
    error('seizmo:fssxc:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif((isnumeric(pv.method) ...
        && (~isreal(pv.method) || ~numel(pv.method)==2)) ...
        || (ischar(pv.method) && ~any(strcmpi(pv.method,valid.METHOD))))
    error('seizmo:fssxc:badInput',...
        ['METHOD must be one of the following:\n' ...
        '''XC'' or ''CAPONXC''!']);
elseif(~isscalar(pv.whiten) || ~islogical(pv.whiten))
    error('seizmo:fssxc:badInput',...
        'WHITEN must be TRUE or FALSE!');
elseif(~isnumeric(pv.w) || any(pv.w(:)<0))
    error('seizmo:fssxc:badInput',...
        'WEIGHTS must be positive!');
elseif(~isreal(pv.damping) || ~isscalar(pv.damping) || pv.damping<0)
    error('seizmo:fssxc:badInput',...
        'DAMPING must be a positive scalar!');
elseif(~isreal(pv.fhwidth) || ~isscalar(pv.fhwidth) ...
        || pv.fhwidth~=fix(pv.fhwidth))
    error('seizmo:fssxc:badInput',...
        'FHWIDTH must be an integer!');
elseif(~isscalar(pv.avg) || ~islogical(pv.avg))
    error('seizmo:fssxc:badInput',...
        'AVG must be TRUE or FALSE!');
end

end


function [lgc]=isxc(data)
[m,s]=getheader(data,'kuser0','kuser1');
if(all(strcmp(m,'MASTER') & strcmp(s,'SLAVE'))); lgc=true;
else lgc=false;
end
end

