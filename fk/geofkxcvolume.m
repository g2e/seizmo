function [svol]=geofkxcvolume(data,ll,s,frng,w)
%GEOFKXCVOLUME    Geographic FK beamforming
%
%    Usage:    sgeo=geofkxcvolume(xcdata,latlon,horzslow,frng)
%              sgeo=geofkxcvolume(xcdata,latlon,horzslow,frng,weights)
%
%    Description: SGEO=GEOFKXCVOLUME(XCDATA,LATLON,HORZSLOW,FRNG) computes
%     the spherical wave coherency (beam strength) through an array as a
%     function of frequency, horizontal slowness & geographic location.
%     The array info and correlograms are contained in the SEIZMO struct
%     XCDATA.  This differs from FKXCVOLUME in that the waves are assumed
%     to be spherical rather than planar.  This is essential for handling
%     sources that are near the array (a good rule of thumb is that if it
%     is within 1 array aperture width of the array then a near-field
%     based fk analysis should be done).  XCDATA must contain correlograms
%     with header formatting as output from CORRELATE.  LATLON contains the
%     latitude and longitude wavefield source positions and must be in
%     units of degrees and formatted as an Nx2 array of [LAT LON].
%     HORZSLOW is the magnitude of the horizontal slowness of the waves in
%     sec/deg and should be a vector.  FRNG gives the frequency range as
%     [FREQLOW FREQHIGH] in Hz.  The output SGEO is a struct containing
%     relevant info and the frequency-slowness-position volume itself (with
%     size NPOSxNSLOWxNFREQ).  The struct layout is:
%          .response - frequency-slowness-position array response
%          .nsta     - number of stations utilized in making map
%          .stla     - station latitudes
%          .stlo     - station longitudes
%          .stel     - station elevations (surface)
%          .stdp     - station depths (from surface)
%          .butc     - UTC start time of data (set to all zeros)
%          .eutc     - UTC end time of data (set to all zeros)
%          .npts     - number of time points
%          .delta    - number of seconds between each time point
%          .latlon   - latitude/longitude source positions (deg)
%          .horzslow - horizontal slowness (sec/deg)
%          .freq     - frequency values (Hz)
%          .normdb   - what 0dB actually corresponds to
%          .volume   - true if frequency-slowness-position volume
%
%     SGEO=GEOFKXCVOLUME(DATA,LATLON,HORZSLOW,FRNG,WEIGHTS) specifies
%     weights for each correlogram in XCDATA (must match size of XCDATA).
%     The weights are normalized internally to sum to 1.
%
%    Notes:
%     - Records in XCDATA must be correlograms following the formatting
%       from CORRELATE.  The records must have the same lag range & sample
%       spacing.
%     - Best/quickest results are obtained when XCDATA is only one
%       "triangle" of the cross correlation matrix.
%
%    Examples:
%     Do you see the 26s microseism in your data?:
%      [lat,lon]=meshgrid();
%      hs=27:0.5:33;
%      frng=[1/27 1/26];
%      sgeo=geofkxcvolume(xcdata,[lat(:) lon(:)],hs,frng);
%
%    See also: FKXCVOLUME, GEOFKXCVOLUME, CHKGEOFKSTRUCT

%     Version History:
%        June 22, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 22, 2010 at 12:55 GMT

% todo:
% - geometrical spreading & Q would be nice

% check nargin
error(nargchk(4,5,nargin));
error(nargchk(1,1,nargout));

% define some constants
d2r=pi/180;
d2km=6371*d2r;

% check struct
versioninfo(data,'dep');

% number of correlograms
ncorr=numel(data);

% defaults for optionals
if(nargin<5 || isempty(w)); w=ones(numel(data),1); end

% check inputs
sf=size(frng);
if(~isreal(ll) || ndims(ll)~=2 || size(ll,2)~=2)
    error('seizmo:geofkxcvolume:badInput',...
        'LATLON must be a Nx2 positive real matrix of [LAT LON]!');
elseif(~isreal(s) || ~isvector(s) || any(s<=0))
    error('seizmo:geofkxcvolume:badInput',...
        'HORZSLOW must be a positive real vector in sec/deg!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:geofkxcvolume:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
elseif(numel(w)~=ncorr || any(w(:)<0) ||  ~isreal(w) || sum(w(:))==0)
    error('seizmo:geofkxcvolume:badInput',...
        'WEIGHTS must be equal sized with XCDATA & be positive numbers!');
end
nrng=sf(1);

% column vector slownesses (convert to sec/km and count)
s=s(:)/d2km;
nslow=numel(s);

% fix lat/lon
[ll(:,1),ll(:,2)]=fixlatlon(ll(:,1),ll(:,2));
nll=size(ll,1);

% convert weights to column vector
w=w(:);
sw=sum(w);

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
    error(lasterror)
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
    
    % check nyquist
    fnyq=1/(2*delta(1));
    if(any(frng>=fnyq))
        error('seizmo:geofkxcvolume:badFRNG',...
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
    [svol(1:nrng,1).volume]=deal(true(1,2));
    [svol(1:nrng,1).latlon]=deal(ll);
    [svol(1:nrng,1).horzslow]=deal(s*d2km);
    
    % get frequencies
    pow2pad=0; % 0 is the default
    nspts=2^(pow2pad+nextpow2(npts(1)));
    f=(0:nspts/2)/(delta(1)*nspts);  % only +freq
    
    % extract data (silently)
    seizmoverbose(false);
    data=splitpad(data,pow2pad);
    data=records2mat(data);
    seizmoverbose(verbose);
    
    % get fft (conjugate is b/c my xc is flipped?)
    % - this is the true cross spectra
    data=conj(fft(data,nspts,1));
    
    % distance difference for the phasors that steer the array
    % dd is NLLxNCORR
    ev=ev.'; st=st.';
    distm=sphericalinv(ll(:,ones(ncorr,1)),ll(:,2*ones(ncorr,1)),...
        ev(ones(nll,1),:),ev(2*ones(nll,1),:));
    dists=sphericalinv(ll(:,ones(ncorr,1)),ll(:,2*ones(ncorr,1)),...
        st(ones(nll,1),:),st(2*ones(nll,1),:));
    dd=2*pi*1i*(distm-dists)*d2km;
    
    % loop over frequency ranges
    for a=1:nrng
        % get frequencies
        fidx=find(f>=frng(a,1) & f<=frng(a,2));
        svol(a).freq=f(fidx);
        nfreq=numel(fidx);
        
        % preallocate fk space
        svol(a).response=zeros(nll,nslow,nfreq,'single');
        
        % warning if no frequencies
        if(~nfreq)
            warning('seizmo:geofkxcvolume:noFreqs',...
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
            % loop over slownesses
            for c=1:nslow
                % get response
                svol(a).response(:,c,b)=10*log10(abs(real(...
                    exp(f(fidx(b))*s(c)*dd)*cs(:,b)))/sw);
            end
            
            % detail message
            if(verbose); print_time_left(b,nfreq); end
        end
        
        % normalize so max peak is at 0dB
        svol(a).normdb=max(svol(a).response(:));
        svol(a).response=svol(a).response-svol(a).normdb;
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
