function [s]=geofss(data,ll,slow,frng,varargin)
%GEOFSS    Estimate frequency-slowness-position spectrum
%
%    Usage:    s=geofss(data,latlon,slow,frng)
%              s=geofss(...,'method',string,...)
%              s=geofss(...,'whiten',true|false,...)
%              s=geofss(...,'weights',w,...)
%              s=geofss(...,'damping',d,...)
%              s=geofss(...,'ntiles',nt,...)
%              s=geofss(...,'fhwidth',n,...)
%              s=geofss(...,'prep',true|false,...)
%              s=geofss(...,'avg',true|false,...)
%
%    Description:
%     S=GEOFSS(DATA,LATLON,SLOW,FRNG) computes an estimate of the
%     frequency-slowness-position power spectra for an array by frequency
%     domain beamforming the time series dataset DATA.  The dataset DATA
%     is a SEIZMO struct containing array info and time series recordings.
%     This function differs from FSS in that the waves are defined as
%     surface waves expanding and contracting on a sphere rather than plane
%     waves traveling across a planar surface.  This is essential for
%     characterizing surface wave sources that are within or near the array
%     (a rule of thumb: a source within one array aperture width has
%     significant near-field terms).  LATLON contains the latitude and
%     longitude wavefield beamforming positions and must be in units of
%     degrees and formatted as an Nx2 array of [LAT LON].  Latitudes are
%     assumed to be geographic (but distances are computed as on a sphere).
%     SLOW is the magnitude of the horizontal slowness of the waves in
%     sec/deg.  FRNG gives the frequency range as [FREQLOW FREQHIGH] in Hz
%     (the actual frequencies are predetermined by the fft of the data).
%     The output S is a struct containing relevant info and the frequency-
%     slowness-position spectra (with size NPOSxNSLOWxNFREQ).  The struct
%     layout is:
%          .nsta     - number of stations
%          .st       - station positions [lat lon elev depth]
%          .butc     - UTC start time of data
%          .eutc     - UTC end time of data
%          .npts     - number of time points
%          .delta    - number of seconds between each time point
%          .latlon   - latitude/longitude positions (deg)
%          .slow     - horizontal slowness (sec/deg)
%          .freq     - frequency values (Hz)
%          .npairs   - number of pairs (aka correlograms)
%          .method   - beamforming method ('center', 'coarray', etc)
%          .center   - array center as [LAT LON]
%          .whiten   - is spectra whitened? true/false
%          .weights  - weights used in beamforming
%          .spectra  - frequency-slowness-position spectra estimate
%
%     S=GEOFSS(...,'METHOD',STRING,...) defines the beaming method.  STRING
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
%     the best resolution out of those three methods but may be inferior to
%     the 'capon' method.  The 'capon' method is similar to the 'full'
%     method in computation time (the only additional operation is a
%     dampened inversion of the full cross power spectral density matrix)
%     but the resolution at every slowness is optimized in a least-squares
%     sense using the maximum likelyhood method.  Because of the reliance
%     of the 'capon' method on the actual data, the signal to incoherent
%     noise levels must be high otherwise the 'coarray' method is the best
%     choice.  Using [LAT LON] as a method is algorithmically the same as
%     the 'center' method but uses the defined coordinates as the center
%     for the array.
%
%     S=GEOFSS(...,'WHITEN',TRUE|FALSE,...) whitens the spectras before
%     beamforming if WHITEN is TRUE.  The default is TRUE.
%
%     S=GEOFSS(...,'WEIGHTS',W,...) specifies the relative weights for each
%     record in DATA (must match size of DATA) or pairing (if METHOD is
%     'coarray', 'capon', or 'full').  You may also specify weighting
%     based on pair geometry using 'ddiff' or 'azdiff' (but only if METHOD
%     is 'coarray' or 'full').  The default is even weighting.
%
%     S=GEOFSS(...,'DAMPING',D,...) alters the dampening parameter used in
%     the inversion of the cross power spectral density matrix for the
%     'capon' method.  This is done by "diagonal loading" which means the
%     dampening value is added the elements along the diagonal.  The
%     default value of 0.001 may need to be adjusted for better results.
%
%     S=GEOFSS(...,'NTILES',NT,...) sets the number of nonoverlapping
%     timesections the data should be split into.  This is mainly for the
%     'capon' method as it needs time/frequency averaging of the cross
%     power spectral density matrix for stability.  Setting NT>=N where N
%     is the number of records in DATA is recommended by Capon 1969.  The
%     default is 1 which provides the best frequency resolution.
%
%     S=GEOFSS(...,'FHWIDTH',N,...) sets the sliding window halfwidth used to
%     average neighboring frequencies of the cross spectral matrix.  This
%     is mainly for stabilizing the inversion in the 'capon' method but has
%     utility in smoothing spectra.  The window is size 2N+1 so for example
%     N=2 averages each frequency with the closest 2 discrete frequencies
%     above and below (a 5 point sliding average).  The default N=0 which
%     does not average.
%
%     S=GEOFSS(...,'PREP',TRUE|FALSE,...) prepares the tiled data by
%     detrending and tapering.  The default is FALSE.
%
%     S=GEOFSS(...,'AVG',TRUE|FALSE,...) indicates if the spectra is
%     averaged across frequency & slowness during computation.  This can
%     save a significant amount of memory.  The default is false.
%
%    Notes:
%     - Records in DATA must have equal and regular sample spacing.
%     - Attenuation correction is not included.
%     - Latitudes are assumed to be geographic!
%     - The slowness input may also be specified as a function that outputs
%       travel time given the following style of input:
%        tt=func([evla evlo ...],[stla stlo ...],freq)
%       where tt is the travel time between the locations.  The lat/lon
%       inputs are Nx2+ arrays and so tt is expected to be a Nx1 array.
%       This allows GEOFSS to beamform just about any simple wave.
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
%     % Compare resolution of a random source with a random array:
%     src=randlatlon;
%     st=randlatlon(30);
%     d=bseizmo([0 1 zeros(1,998)]);
%     d=changeheader(d(ones(size(st,1),1)),...
%       'stla',st(:,1),'stlo',st(:,2),...
%       'b',sphericalinv(src(:,1),src(:,2),st(:,1),st(:,2))*30);
%     snr=20;
%     d=solofun(d,@(x)x+(rand(1000,1)-.5)/snr);
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     ll=[lat(:) lon(:)]; clear lat lon;
%     frng=[.0001 .002];
%     plotgeofss(geofss(d,ll,30,frng,'a',true));
%     plotgeofss(geofss(d,ll,30,frng,'a',true,'m','coarray'));
%     plotgeofss(geofss(d,ll,30,frng,'a',true,'m','coarray','w','ddiff'));
%     plotgeofss(geofss(d,ll,30,frng,'a',true,'m','coarray','w','azdiff'));
%     plotgeofss(geofss(d,ll,30,frng,'a',true,'m','capon'));
%
%    See also: GEOFSSXC, FSSXC, GEOFSSAVG, GEOFSSSUB, PLOTGEOFSS, GEOARF
%              GEOFSSDBINFO, GEOFSSCORRCOEF, GEOFSSFRAMESLIDE,
%              GEOFSSFREQSLIDE, GEOFSSSLOWSLIDE, PLOTGEOARF, GEOFSSHORZ,
%              GEOFSSHORZXC, GEOTDSSXC, GEOTDSSHORZXC

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
%        June 13, 2012 - capon method added (needs work)
%        Aug. 30, 2012 - cleaned out some deprecated comments
%        Sep. 11, 2012 - handle latitudes properly
%        Oct.  1, 2012 - big rewrite: pv pair inputs, tt function input,
%                        special weighting schemes, avg & ntiles options,
%                        capon method with dampening, major code reordering
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  1, 2012 at 14:05 GMT

% todo:
% - examples
% - everything is antipodal!

% check nargin
error(nargchk(4,inf,nargin));

% check struct
error(seizmocheck(data,'dep'));

% require not a xc dataset
if(isxc(data))
    error('seizmo:geofss:badInput',...
        'Use GEOFSSXC for correlations!');
end

% number of records
nrecs=numel(data);

% need 2+ records
if(nrecs<2)
    error('seizmo:geofss:arrayTooSmall',...
        'DATA must have 2+ records!');
end

% check inputs
if(~isreal(ll) || ndims(ll)~=2 || size(ll,2)<2 || size(ll,1)<1)
    error('seizmo:geofss:badInput',...
        'LATLON must be a Nx2 real matrix of [LAT LON]!');
elseif(~isa(slow,'function_handle') ...
        && (~isreal(slow) || any(slow<=0) || numel(slow)<1))
    error('seizmo:geofss:badInput',...
        'SLOW must be positive real vector in s/deg!');
end
sf=size(frng);
if(~isreal(frng) || numel(sf)~=2 || any(frng(:)<0))
    error('seizmo:geofss:badInput',...
        'FRNG must be 1 or more positive real values in Hz!');
end
if(sf(2)==1); frng=frng(:,[1 1]); sf(2)=2; end
if(sf(2)~=2 || any(frng(:,1)>frng(:,2)))
    error('seizmo:geofss:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
end
nrng=sf(1);

% slowness function logical
sisafunc=isa(slow,'function_handle');

% parse options
pv=parse_geofss_pv_pairs(varargin{:});

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
    [clalo(1),clalo(2)]=arraycenter(st(:,1),st(:,2));
    clalo=[clalo mean(st(:,3:end))];
    switch pv.method
        case 'coarray'
            [master,slave]=find(triu(true(nrecs),1));
            npairs=nrecs*(nrecs-1)/2;
        case {'full' 'capon'}
            [master,slave]=find(true(nrecs));
            npairs=nrecs*nrecs;
        case 'center'
            npairs=nrecs;
    end
else
    clalo=pv.method;
    pv.method='user';
    npairs=nrecs;
end

% check weights again
if(isnumeric(pv.w) && ~any(numel(pv.w)==[nrecs npairs]))
    error('seizmo:geofss:badInput',...
        'Number of WEIGHTS must match the number of stations or pairs!');
end

% get time limits
[bmin,bmini]=min(timediff(butc(1,:),butc,'utc'));
[emax,emaxi]=max(timediff(eutc(1,:),eutc,'utc'));

% check nyquist
fnyq=1/(2*delta(1));
if(any(frng(:,1)>fnyq))
    error('seizmo:geofss:badFRNG',...
        ['FRNG exceeds nyquist frequency (' num2str(fnyq) ')!']);
end

% longest record
maxnpts=max(npts);

% tiling
if(pv.ntiles>1)
    oldmax=maxnpts;
    maxnpts=ceil(maxnpts/pv.ntiles);
end

% get frequencies
nspts=2^nextpow2(maxnpts); % half xcorr
%nspts=2^nextpow2(2*maxnpts-1); % full xcorr for verification
f=(0:nspts/2)/(delta(1)*nspts);  % only +freq
nf=numel(f);

% tiling
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

% column vector slownesses
if(sisafunc)
    nslow=1;
else % explicit
    slow=slow(:);
    nslow=numel(slow);
end

% fix lat/lon
[ll(:,1),ll(:,2)]=fixlatlon(ll(:,1),ll(:,2));
nll=size(ll,1);

% setup output
[s(1:nrng,1).nsta]=deal(nrecs);
[s(1:nrng,1).st]=deal(st);
[s(1:nrng,1).butc]=deal(butc(bmini,:));
[s(1:nrng,1).eutc]=deal(eutc(emaxi,:));
[s(1:nrng,1).delta]=deal(delta(1));
[s(1:nrng,1).npts]=deal(maxnpts);
[s(1:nrng,1).latlon]=deal(ll);
[s(1:nrng,1).slow]=deal(slow);
[s(1:nrng,1).freq]=deal([]);
[s(1:nrng,1).method]=deal(pv.method);
[s(1:nrng,1).npairs]=deal(npairs);
[s(1:nrng,1).center]=deal(clalo);
[s(1:nrng,1).whiten]=deal(pv.whiten);
[s(1:nrng,1).weights]=deal(pv.w);
[s(1:nrng,1).ntiles]=deal(pv.ntiles);
[s(1:nrng,1).fhwidth]=deal(pv.fhwidth);
[s(1:nrng,1).damping]=deal(pv.damping);
[s(1:nrng,1).spectra]=deal(zeros(0,'single'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of defaults, checks, & data setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert latitudes
ll(:,1)=geographic2geocentriclat(ll(:,1));
st(:,1)=geographic2geocentriclat(st(:,1));

% get utc static time shifts
dt=timediff(butc(1,:),butc,'utc').'; % need row vector
switch pv.method
    case {'coarray' 'full' 'capon'}
        dt=dt(slave)-dt(master); % alter to pair shifts
end

% expand & normalize weights
if(isnumeric(pv.w))
    pv.w=pv.w(:);
    switch pv.method
        case {'coarray' 'full' 'capon'}
            if(numel(pv.w)==nrecs)
                pv.w=pv.w(slave).*conj(pv.w(master));
            end
    end
    pv.w=pv.w./sum(abs(pv.w));
    z=1;
else
    switch pv.method
        case {'user' 'center'}
            error('seizmo:geofss:badInput',...
                'Cannot use a pair weighting scheme for solo methods!');
    end
    switch pv.w
        case 'azdiff'
            [az,az]=sphericalinv(...
                ll(:,ones(nrecs,1)),ll(:,2*ones(nrecs,1)),...
                st(:,ones(nll,1))',st(:,2*ones(nll,1))');
            pv.w=ones(npairs,1);
            z=abs(azdiff(az(:,master),az(:,slave)));
            z=z./repmat(sum(z,2),1,npairs);
        case 'ddiff'
            pv.w=sphericalinv(...
                st(master,1),st(master,2),st(slave,1),st(slave,2));
            pv.w=pv.w./sum(abs(pv.w));
            z=1;
    end
end

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
    if(pv.avg); s(a).spectra=zeros(nll,1,'single');
    else s(a).spectra=zeros(nll,nslow,nfreq,'single');
    end
    
    % error if no frequencies
    if(~nfreq)
        error('seizmo:geofss:noFreqs',...
            'No frequencies within the range %g to %g Hz!',...
            frng(a,1),frng(a,2));
    end
    
    % detail message
    if(verbose)
        fprintf('Getting spectra %d: %g to %g Hz\n',a,frng(a,1),frng(a,2));
        print_time_left(0,nfreq);
    end
    
    % loop over frequency
    for b=1:nfreq
        % Get distance or traveltime from each
        % position of interest to every station
        % tt is NLLxNPAIRS
        if((a==1 && b==1) || ttredo)
            [tt,ctt,shift,cshift,ttredo]=gettt(...
                ll,st,clalo,slow,f(fidx),b,nll,nrecs,sisafunc,pv);
            % get degree distance difference or traveltime difference
            switch pv.method
                case {'coarray' 'full' 'capon'}
                    tt=-tt(:,slave)+tt(:,master)+dt(ones(nll,1),:);
                    shift=shift(:,slave)-shift(:,master);
                otherwise % {'center' 'user'}
                    tt=-tt+ctt(:,ones(1,nrecs))+dt(ones(nll,1),:);
                    shift=shift-cshift(:,ones(1,nrecs));
            end
        end
        
        % loop over slowness
        for c=1:nslow
            % modify distance to traveltime?
            if(sisafunc); hs=1; else hs=slow(c); end
            
            % beamforming method
            switch pv.method
                case {'coarray' 'full'}
                    if(pv.avg)
                        % get frequency range
                        newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                        newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                        
                        % compute average cross spectra for freq range & tiles
                        % then weighted beam and add
                        s(a).spectra=s(a).spectra+real((z.*exp(-2*pi*1i*f(fidx(b))*(hs*tt)-1i*shift))...
                            *(mean(mean(data(slave,newfidx,:).*conj(data(master,newfidx,:)),2),3).*pv.w));
                    else
                        % get frequency range
                        newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                        newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                        
                        % compute average cross spectra for freq range & tiles
                        % then weighted beam and insert
                        s(a).spectra(:,c,b)=real((z.*exp(-2*pi*1i*f(fidx(b))*(hs*tt)-1i*shift))...
                            *(mean(mean(data(slave,newfidx,:).*conj(data(master,newfidx,:)),2),3).*pv.w));
                    end
                case 'capon'
                    if(pv.avg)
                        % get frequency range
                        newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                        newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                        
                        % get average cross spectral matrix
                        cs=reshape(mean(mean(data(slave,newfidx,:).*conj(...
                            data(master,newfidx,:)),2),3),[nrecs nrecs]);
                        
                        % damped inversion of cross spectral matrix
                        cs=pinv((1-pv.damping)*cs+pv.damping*eye(nrecs));
                        
                        % weighted capon beam and add
                        s(a).spectra=s(a).spectra+1./real((z.*exp(-2*pi*1i*f(fidx(b))*(hs*tt)-1i*shift))...
                            *(cs(:).*pv.w));
                    else
                        % get frequency range
                        newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                        newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                        
                        % get average cross spectral matrix
                        cs=reshape(mean(mean(data(slave,newfidx,:).*conj(...
                            data(master,newfidx,:)),2),3),[nrecs nrecs]);
                        
                        % damped inversion of cross spectral matrix
                        cs=pinv((1-pv.damping)*cs+pv.damping*eye(nrecs));
                        
                        % weighted capon beam and insert
                        s(a).spectra(:,c,b)=1./real((z.*exp(-2*pi*1i*f(fidx(b))*(hs*tt)-1i*shift))...
                            *(cs(:).*pv.w));
                    end
                otherwise  % {'center' 'user'}
                    if(pv.avg)
                        % get frequency range
                        newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                        newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                        newnf=numel(newfidx);
                        
                        for d=1:pv.ntiles
                            s(a).spectra=s(a).spectra+mean(abs((z.*exp(-2*pi*1i*f(fidx(b))*(hs*tt)-1i*shift))...
                                *(data(:,newfidx,d).*pv.w(:,ones(1,newnf)))).^2,2);
                        end
                    else
                        % get frequency range
                        newfidx=fidx(b)+(-2*pv.fhwidth:2*pv.fhwidth);
                        newfidx=newfidx(newfidx>=1 & newfidx<=nf);
                        newnf=numel(newfidx);
                        
                        for d=1:pv.ntiles
                            s(a).spectra(:,c,b)=s(a).spectra(:,c,b)+mean(abs((z.*exp(-2*pi*1i*f(fidx(b))*(hs*tt)-1i*shift))...
                                *(data(:,newfidx,d).*pv.w(:,ones(1,newnf)))).^2,2);
                        end
                    end
                    s(a).spectra=s(a).spectra/pv.ntiles;
            end
        end
        if(verbose); print_time_left(b,nfreq); end
    end
    if(pv.avg); s(a).spectra=s(a).spectra/(nfreq*nslow); end
end

end

function [tt,ctt,shift,cshift,ttredo]=gettt(...
    lalo,stlalo,clalo,slow,f0,a,nll,nrecs,sisafunc,pv)
% super overcomplex distance/traveltime retriever
ctt=[];
cshift=[];
ttredo=true;
if(a==1)
    % function gives traveltime but numeric input gives slowness
    if(sisafunc)
        % travel time from function
        if(isscalar(unique(f0)))
            % only one frequency so only need to get tt once
            tt=nan(nll,nrecs);
            shift=nan(nll,nrecs);
            for i=1:nrecs
                [tt(:,i),shift(:,i)]=slow(...
                    lalo,stlalo(i*ones(nll,1),:),f0(1));
            end
            switch pv.method
                case {'user' 'center'}
                    [ctt,cshift]=slow(lalo,clalo,f0(1));
            end
            %ttredo=false;  % MIGHT NEED TO REDO FOR ANOTHER FREQ RANGE
            ttredo=false;
        else % multiple frequencies
            try
                % see if function depends on f
                [tt,shift]=slow(lalo(1,:),stlalo(1,:));
                
                % still here?
                % guess that means frequency is ignored...
                tt=nan(nll,nrecs);
                shift=nan(nll,nrecs);
                for i=1:nrecs
                    [tt(:,i),shift(:,i)]=slow(...
                        lalo,stlalo(i*ones(nll,1),:));
                end
                switch pv.method
                    case {'user' 'center'}
                        [ctt,cshift]=slow(lalo,clalo);
                end
                ttredo=false;
            catch
                % slow is a function of frequency
                tt=nan(nll,nrecs);
                shift=nan(nll,nrecs);
                for i=1:nrecs
                    [tt(:,i),shift(:,i)]=slow(...
                        lalo,stlalo(i*ones(nll,1),:),f0(a));
                end
                switch pv.method
                    case {'user' 'center'}
                        [ctt,cshift]=slow(lalo,clalo,f0(a));
                end
            end
        end
    else
        % degree distance (scaled by the slowness gives the traveltime)
        tt=sphericalinv(...
            lalo(:,ones(nrecs,1)),lalo(:,2*ones(nrecs,1)),...
            stlalo(:,ones(nll,1))',stlalo(:,2*ones(nll,1))');
        shift=zeros(size(tt));
        ctt=sphericalinv(lalo(:,1),lalo(:,2),clalo(1),clalo(2));
        cshift=zeros(size(ctt));
        ttredo=false;
    end
else
    % slow is a function of frequency
    tt=nan(nll,nrecs);
    shift=nan(nll,nrecs);
    for i=1:nrecs
        [tt(:,i),shift(:,i)]=slow(lalo,stlalo(i*ones(nll,1),:),f0(a));
    end
    switch pv.method
        case {'user' 'center'}
            [ctt,cshift]=slow(lalo,clalo,f0(a));
    end
end

end


function [pv]=parse_geofss_pv_pairs(varargin)

% defaults
pv.avg=false;
pv.method='center';
pv.whiten=true;
pv.w=[];
pv.damping=0.001; % only for capon
pv.ntiles=1;
pv.fhwidth=0;
pv.prep=false;

% require pv pairs
if(mod(nargin,2))
    error('seizmo:geofss:badInput',...
        'Unpaired parameter/value!');
end

% parse pv pairs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:geofss:badInput',...
        'Parameters must be specified using a string!');
end
for i=1:2:nargin
    switch varargin{i}
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
            error('seizmo:geofss:badInput',...
                'Unknown parameter: %s !',varargin{i});
    end
end

% valid method strings
valid.METHOD={'center' 'coarray' 'full' 'capon'};
valid.W={'azdiff' 'ddiff'};

% check values
if((isnumeric(pv.method) ...
        && (~isreal(pv.method) || ~numel(pv.method)==2)) ...
        || (ischar(pv.method) && ~any(strcmpi(pv.method,valid.METHOD))))
    error('seizmo:geofss:badInput',...
        ['METHOD must be one of the following:\n' ...
        '''CENTER'', ''COARRAY'', ''CAPON'', ''FULL'', or [LAT LON]!']);
elseif(~isscalar(pv.whiten) || ~islogical(pv.whiten))
    error('seizmo:geofss:badInput',...
        'WHITEN must be TRUE or FALSE!');
elseif((~isnumeric(pv.w) && ~ischar(pv.w)) ...
        || (isnumeric(pv.w) && any(pv.w(:)<0)) ...
        || (ischar(pv.w) && ~any(strcmpi(pv.w,valid.W))))
    error('seizmo:geofss:badInput',...
        'WEIGHTS must be positive!');
elseif(~isreal(pv.damping) || ~isscalar(pv.damping) || pv.damping<0)
    error('seizmo:geofss:badInput',...
        'DAMPING must be a positive scalar!');
elseif(~isreal(pv.ntiles) || ~isscalar(pv.ntiles) ...
        || pv.ntiles<=0 || pv.ntiles~=fix(pv.ntiles))
    error('seizmo:geofss:badInput',...
        'NTILES must be a positive integer!');
elseif(~isreal(pv.fhwidth) || ~isscalar(pv.fhwidth) ...
        || pv.fhwidth~=fix(pv.fhwidth))
    error('seizmo:geofss:badInput',...
        'FHWIDTH must be an integer!');
elseif(~isscalar(pv.avg) || ~islogical(pv.avg))
    error('seizmo:geofss:badInput',...
        'AVG must be TRUE or FALSE!');
elseif(~isscalar(pv.prep) || ~islogical(pv.prep))
    error('seizmo:geofss:badInput',...
        'PREP must be TRUE or FALSE!');
end

end


function [lgc]=isxc(data)
[m,s]=getheader(data,'kuser0','kuser1');
if(all(strcmp(m,'MASTER') & strcmp(s,'SLAVE'))); lgc=true;
else lgc=false;
end
end
