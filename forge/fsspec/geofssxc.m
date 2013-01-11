function [s]=geofssxc(xcdata,ll,slow,frng,varargin)
%GEOFSSXC    Estimate frequency-slowness-position spectrum
%
%    Usage:    s=geofssxc(xcdata,latlon,slow,frng)
%              s=geofssxc(...,'method',string,...)
%              s=geofssxc(...,'whiten',true|false,...)
%              s=geofssxc(...,'weights',w,...)
%              s=geofssxc(...,'damping',d,...)
%              s=geofssxc(...,'avg',true|false,...)
%
%    Description:
%     S=GEOFSSXC(XCDATA,LATLON,SLOW,FRNG) computes an estimate of the
%     frequency-slowness-position power spectra for an array by frequency
%     domain beamforming the time series data XCDATA.  The dataset XCDATA
%     is a SEIZMO struct containing array info and precomputed cross
%     correlograms.  This function differs from FSSXC in that the waves are
%     defined as surface waves expanding and contracting on a sphere
%     rather than plane waves traveling across a planar surface.  This is
%     essential for characterizing surface wave sources that are within or
%     near the array (a rule of thumb: a source within one array aperture
%     width has significant near-field terms).  XCDATA is expected to be
%     correlograms with header formatting as output from CORRELATE.
%     LATLON contains the latitude and longitude wavefield beamforming
%     positions and must be in units of degrees and formatted as an Nx2
%     array of [LAT LON].  SLOW is the magnitude of the horizontal slowness
%     of the waves in sec/deg.  FRNG gives the frequency range as
%     [FREQLOW FREQHIGH] in Hz (the actual frequencies are predetermined by
%     the fft of the data).  The output S is a struct containing relevant
%     info and the frequency-slowness-position spectra (with size
%     NPOSxNSLOWxNFREQ).  The struct layout is:
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
%          .method   - beamforming method ('xc' or 'caponxc')
%          .center   - array center as [LAT LON]
%          .whiten   - is spectra whitened? true/false
%          .weights  - weights used in beamforming
%          .spectra  - frequency-slowness-position spectra estimate
%
%     S=GEOFSSXC(...,'METHOD',STRING,...) specifies the beamforming method.
%     STRING may be either 'xc' or 'caponxc'.  The default is 'xc' which
%     uses the correlations as a substitute for the cross power spectral
%     density matrix in a conventional frequency-slowness estimation (eg.
%     'coarray' or 'full').  The 'caponxc' method requires that the
%     correlation dataset can form the entire cross power spectral density
%     matrix (ie. this requires all N^2 pairs).  If possible then this
%     method can return a higher resolution estimate of the frequency-
%     slowness spectra.
%
%     S=GEOFSSXC(...,'WHITEN',TRUE|FALSE,...) whitens the cross power
%     spectral matrix elements before beamforming if WHITEN is TRUE.  The
%     default is TRUE.
%
%     S=GEOFSSXC(...,'WEIGHTS',W,...) sets the relative weights for each
%     correlogram in XCDATA (must match size of XCDATA).  You may also
%     specify weighting based on pair geometry using 'ddiff' or 'azdiff'
%     (but only if METHOD is 'xc').  The default is even weighting.
%
%     S=GEOFSSXC(...,'DAMPING',D,...) alters the dampening parameter used
%     in the inversion of the cross power spectral density matrix for the
%     'capon' method.  This is done by diagonal loading which means the
%     dampening value is added the elements along the diagonal.  The
%     default value of 0.001 may need to be adjusted for better results.
%
%     S=GEOFSSXC(...,'AVG',TRUE|FALSE,...) indicates if the spectra is
%     averaged across frequency during computation.  This can save a
%     significant amount of memory.  The default is false.
%
%    Notes:
%     - Correlations in XCDATA must have equal and regular sample spacing.
%     - Correlations are expected to have station naming and location info
%       as stored in the header by the function CORRELATE.  This means the
%       headers of correlations done outside of SEIZMO will probably
%       require adjustment.
%     - Best/quickest results for the 'xc' method are obtained when XCDATA
%       is only one "triangle" of the cross correlation matrix (no auto-
%       correlation terms).  This corresponds to the 'coarray' method of
%       GEOFSS.
%     - Latitudes are assumed to be geographic!
%     - The slowness input may also be specified as a function that outputs
%       travel time given the following style of input:
%        tt=func([evla evlo ...],[stla stlo ...],freq)
%       where tt is the travel time between the locations.  The lat/lon
%       inputs are Nx2+ arrays and so tt is expected to be a Nx1 array.
%       This allows GEOFSSXC to beamform just about any simple wave.
%     - References:
%        Capon 1969, High-Resolution Frequency-Wavenumber Spectrum
%         Analysis, Proc. IEEE, Vol. 57, No. 8, pp. 1408-1418
%
%    Examples:
%     % Artificial lake source detection by LASA:
%     src=[47.9 -106.5];
%     [st,nm]=lasa;
%     d=bseizmo([0 1 zeros(1,98)]);
%     d=changeheader(d(ones(size(st,1),1)),...
%       'stla',st(:,1),'stlo',st(:,2),'kstnm',nm,...
%       'b',sphericalinv(src(1),src(2),st(:,1),st(:,2))*30);
%     snr=20;
%     d=solofun(d,@(x)x+(rand(100,1)-.5)/snr);
%     [lat,lon]=meshgrid(45:.025:48.5,-110:.025:-104.5);
%     plotgeofss(geofssxc(correlate(d),[lat(:) lon(:)],30,[1/4 1/3],...
%         'a',true,'w','azdiff'),[-3 0],[],'g','i');
%
%    See also: GEOFSS, GEOFSSAVG, GEOFSSSUB, PLOTGEOFSS, GEOFSSHORZXC,
%              GEOARF, GEOFSSDBINFO, GEOFSSCORRCOEF, GEOFSSFRAMESLIDE,
%              GEOFSSFREQSLIDE, GEOFSSSLOWSLIDE, PLOTGEOARF, GEOFSSHORZ

%     Version History:
%        June 22, 2010 - initial version
%        July  6, 2010 - major update to struct, doc update
%        July  7, 2010 - removed deg to km conversions
%        Apr.  3, 2012 - use seizmocheck
%        May  30, 2012 - pow2pad=0 by default
%        June  4, 2012 - altered from geofkxcvolume
%        June 10, 2012 - no more conj (multiply shift by -1), verified test
%                        results
%        June 12, 2012 - add whiten option
%        Oct.  1, 2012 - big rewrite: pv pair inputs, tt function input,
%                        special weighting schemes, avg option, capon
%                        method with dampening, major code reordering
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  1, 2012 at 14:05 GMT

% todo:

% check nargin
error(nargchk(4,inf,nargin));

% check struct
error(seizmocheck(xcdata,'dep'));

% require xc dataset
if(~isxc(xcdata))
    error('seizmo:geofssxc:badInput',...
        'XCDATA must be correlations with metadata as from CORRELATE!');
end

% number of pairings
npairs=numel(xcdata);

% check inputs
sf=size(frng);
if(~isreal(ll) || ndims(ll)~=2 || size(ll,2)~=2)
    error('seizmo:geofssxc:badInput',...
        'LATLON must be a Nx2 real matrix of [LAT LON]!');
elseif(~isreal(slow) || ~isvector(slow) || any(slow<=0))
    error('seizmo:geofssxc:badInput',...
        'SLOW must be a positive real vector in sec/deg!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:geofssxc:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
end
nrng=sf(1);

% parse options
pv=parse_geofssxc_pv_pairs(varargin{:});

% slowness function logical
sisafunc=isa(slow,'function_handle');

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

% extract header info & xcdata
try
    % verbosity
    verbose=seizmoverbose;
    
    % grab necessary header info
    [b,npts,delta,autc,futc,st,ev,stnm,evnm]=getheader(xcdata,...
        'b','npts','delta','a utc','f utc','st','ev','kname','kt');
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

% get fft
xcdata=fft(xcdata,nspts,1);
xcdata=xcdata(1+(0:nspts/2),:).'; % trim -freq

% whiten data if desired
if(pv.whiten); xcdata=xcdata./abs(xcdata); end

% use unique stations names to get number of stations & locations
evnm=evnm(:,1:4);
stnm=strcat(stnm(:,1),'.',stnm(:,2),'.',stnm(:,3),'.',stnm(:,4));
evnm=strcat(evnm(:,1),'.',evnm(:,2),'.',evnm(:,3),'.',evnm(:,4));
[stnm,idx1,idx2]=unique([stnm;evnm]);
st=[st;ev];
st=st(idx1,:);
nrecs=numel(stnm);

% check weights again
if(isnumeric(pv.w) && ~any(numel(pv.w)==[nrecs npairs]))
    error('seizmo:geofssxc:badInput',...
        '# of WEIGHTS must match the # of stations or records in XCDATA!');
end

% slave/master indices
idx2=reshape(idx2,[],2);
slave=idx2(:,1);
master=idx2(:,2);

% check if full matrix for caponxc
if(strcmpi(pv.method,'caponxc'))
    tmp=false(nrecs);
    csidx=sub2ind([nrecs nrecs],master,slave);
    tmp(csidx)=true;
    if(npairs~=nrecs^2 || ~all(tmp(:)))
        error('seizmo:geofssxc:badInput',...
            'XCDATA does not have all pairs necessary for CAPONXC!');
    end
end

% array center
[clalo(1),clalo(2)]=arraycenter(st(:,1),st(:,2));
clalo=[clalo mean(st(:,3:end))];

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
[s(1:nrng,1).butc]=deal(autc);
[s(1:nrng,1).eutc]=deal(futc);
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
[s(1:nrng,1).spectra]=deal(zeros(0,'single'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of defaults, checks, & data setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert latitudes
ll(:,1)=geographic2geocentriclat(ll(:,1));
st(:,1)=geographic2geocentriclat(st(:,1));

% get static time shifts
% - MUST shift to zero lag
dt=b.';

% expand & normalize weights
if(isnumeric(pv.w))
    pv.w=pv.w(:);
    if(numel(pv.w)==nrecs); pv.w=pv.w(slave).*conj(pv.w(master)); end
    pv.w=pv.w./sum(abs(pv.w));
    z=1;
else
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
    fidx=find(f>=frng(a,1) & f<=frng(a,2));
    s(a).freq=f(fidx);
    nfreq=numel(fidx);
    
    % preallocate spectra
    if(pv.avg); s(a).spectra=zeros(nll,1,'single');
    else s(a).spectra=zeros(nll,nslow,nfreq,'single');
    end
    
    % error if no frequencies
    if(~nfreq)
        error('seizmo:geofssxc:noFreqs',...
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
            [tt,ttredo]=gettt(ll,st,slow,f(fidx),b,nll,nrecs,sisafunc);
            % get degree distance difference or traveltime difference
            tt=tt(:,slave)-tt(:,master);
        end
        
        % invert cross power spectral density matrix for capon
        switch pv.method
            case 'caponxc'
                cs=zeros(nrecs,nrecs);
                cs(csidx)=xcdata(:,fidx(b)); % populate matrix
                cs=pinv((1-pv.damping)*cs+pv.damping*eye(nrecs));
        end
        
        % loop over slowness
        for c=1:nslow
            % modify distance to traveltime?
            if(sisafunc); hs=1; else hs=slow(c); end
            
            % beamforming method
            switch pv.method
                case 'xc'
                    if(pv.avg)
                        s(a).spectra=s(a).spectra+real(z.*exp(...
                            2*pi*1i*f(fidx(b))*(hs*tt-dt(ones(nll,1),:)))...
                            *(xcdata(:,fidx(b)).*pv.w));
                    else
                        s(a).spectra(:,c,b)=real(z.*exp(...
                            2*pi*1i*f(fidx(b))*(hs*tt-dt(ones(nll,1),:)))...
                            *(xcdata(:,fidx(b)).*pv.w));
                    end
                case 'caponxc'
                    if(pv.avg)
                        s(a).spectra=s(a).spectra+1./real(z.*exp(...
                            2*pi*1i*f(fidx(b))*(hs*tt-dt(ones(nll,1),:)))...
                            *(cs(csidx).*pv.w));
                    else
                        s(a).spectra(:,c,b)=1./real(z.*exp(...
                            2*pi*1i*f(fidx(b))*(hs*tt-dt(ones(nll,1),:)))...
                            *(cs(csidx).*pv.w));
                    end
            end
        end
        if(verbose); print_time_left(b,nfreq); end
    end
    if(pv.avg); s(a).spectra=s(a).spectra/(nfreq*nslow); end
end

end

function [tt,ttredo]=gettt(lalo,stlalo,slow,f0,a,nll,nrecs,sisafunc)
% super overcomplex distance/traveltime retriever
ttredo=true;
if(a==1)
    % function gives traveltime but numeric input gives slowness
    if(sisafunc)
        % travel time from function
        if(isscalar(unique(f0)))
            % only one frequency so only need to get tt once
            tt=nan(nll,nrecs);
            for i=1:nrecs
                tt(:,i)=slow(lalo,stlalo(i*ones(nll,1),:),f0(1));
            end
            %ttredo=false;  % MIGHT NEED TO REDO FOR ANOTHER FREQ RANGE
        else % multiple frequencies
            try
                % see if function depends on f
                tt=slow(lalo(1,:),stlalo(1,:)); %#ok<NASGU>
                
                % still here?
                % guess that means frequency is ignored...
                tt=nan(nll,nrecs);
                for i=1:nrecs
                    tt(:,i)=slow(lalo,stlalo(i*ones(nll,1),:));
                end
                ttredo=false;
            catch
                % slow is a function of frequency
                tt=nan(nll,nrecs);
                for i=1:nrecs
                    tt(:,i)=slow(lalo,stlalo(i*ones(nll,1),:),f0(a));
                end
            end
        end
    else
        % degree distance (scaled by the slowness gives the traveltime)
        tt=sphericalinv(...
            lalo(:,ones(nrecs,1)),lalo(:,2*ones(nrecs,1)),...
            stlalo(:,ones(nll,1))',stlalo(:,2*ones(nll,1))');
        ttredo=false;
    end
else
    % slow is a function of frequency
    tt=nan(nll,nrecs);
    for i=1:nrecs; tt(:,i)=slow(lalo,stlalo(i*ones(nll,1),:),f0(a)); end
end

end


function [pv]=parse_geofssxc_pv_pairs(varargin)

% defaults
pv.avg=false;
pv.method='xc';
pv.whiten=true;
pv.w=[];
pv.damping=0.001; % only for capon

% require pv pairs
if(mod(nargin,2))
    error('seizmo:geofssxc:badInput',...
        'Unpaired parameter/value!');
end

% parse pv pairs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:geofssxc:badInput',...
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
        case {'average' 'av' 'avg' 'a'}
            pv.avg=varargin{i+1};
        otherwise
            error('seizmo:geofssxc:badInput',...
                'Unknown parameter: %s !',varargin{i});
    end
end

% valid method strings
valid.METHOD={'xc' 'caponxc'};
valid.W={'azdiff' 'ddiff'};

% check values
if((isnumeric(pv.method) ...
        && (~isreal(pv.method) || ~numel(pv.method)==2)) ...
        || (ischar(pv.method) && ~any(strcmpi(pv.method,valid.METHOD))))
    error('seizmo:geofssxc:badInput',...
        ['METHOD must be one of the following:\n' ...
        '''XC'' or ''CAPONXC''!']);
elseif(~isscalar(pv.whiten) || ~islogical(pv.whiten))
    error('seizmo:geofssxc:badInput',...
        'WHITEN must be TRUE or FALSE!');
elseif((~isnumeric(pv.w) && ~ischar(pv.w)) ...
        || (isnumeric(pv.w) && any(pv.w(:)<0)) ...
        || (ischar(pv.w) && ~any(strcmpi(pv.method,valid.METHOD))))
    error('seizmo:geofssxc:badInput',...
        'WEIGHTS must be positive!');
elseif(~isreal(pv.damping) || ~isscalar(pv.damping) || pv.damping<0)
    error('seizmo:geofssxc:badInput',...
        'DAMPING must be a positive scalar!');
elseif(~isscalar(pv.avg) || ~islogical(pv.avg))
    error('seizmo:geofssxc:badInput',...
        'AVG must be TRUE or FALSE!');
end

end


function [lgc]=isxc(data)
[m,s]=getheader(data,'kuser0','kuser1');
if(all(strcmp(m,'MASTER') & strcmp(s,'SLAVE'))); lgc=true;
else lgc=false;
end
end

