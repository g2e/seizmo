function [r,t]=fsshorz(n,e,smax,spts,frng,varargin)
%FSSHORZ    Estimate frequency-slowness spectrum of horizontals
%
%    Usage:   [r,t]=fsshorz(n,e,smax,spts,frng)
%             [r,t]=fsshorz(...,'polar',true|false,...)
%             [r,t]=fsshorz(...,'method',string,...)
%             [r,t]=fsshorz(...,'whiten',true|false,...)
%             [r,t]=fsshorz(...,'weights',w,...)
%             [r,t]=fsshorz(...,'ntiles',nt,...)
%             [r,t]=fsshorz(...,'avg',true|false,...)
%
%    Description:
%     [R,T]=FSSHORZ(N,E,SMAX,SPTS,FRNG) computes an estimate of the radial
%     and tranverse frequency-slowness power spectra in a cartesian grid
%     for a spatial array of stations recording horizontal motions by
%     frequency-domain beamforming the time series data in N & E.  The
%     datasets N & E are SEIZMO structs containing the array info and time
%     series recordings of motion in the North and East directions.  This
%     function differs from GEOFSSHORZ in that the waves are assumed to be
%     plane waves traveling on a planar surface rather than surface waves
%     expanding and contracting on a sphere.  The range of the horizontal
%     slowness grid is given by SMAX (sec/deg) and extends from -SMAX to
%     SMAX for both East/West and North/South directions.  SPTS controls
%     the number of slowness points for both directions (SPTSxSPTS grid).
%     FRNG gives the frequency range as [FREQLOW FREQHIGH] in Hz (the
%     individual frequencies are determined by the fft of the data).  The
%     outputs R & T are structs containing relevant info and the frequency-
%     slowness spectra for the radial and transverse wavefields (with size
%     SPTSxSPTSxNFREQ).  The struct layout is:
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
%          .spectra  - frequency-slowness spectra estimate
%          .orient   - 'radial' or 'transverse'
%
%     [R,T]=FSSHORZ(...,'POLAR',TRUE|FALSE,...) specifies if the slownesses
%     are sampled regularly in cartesian or polar coordinates.  Polar
%     coords are useful for slicing the spectra by azimuth (pie slice) or
%     slowness (rings).  Cartesian coords (the default) samples the
%     slowness space regularly in the East/West & North/South directions
%     and so exhibits less distortion in plots of the slowness space. If
%     POLAR=TRUE, SPTS may be given as [SPTS BAZPTS] to control the
%     azimuthal resolution (default is BAZPTS=181 points).
%
%     [R,T]=FSSHORZ(...,'METHOD',STRING,...) sets the beamforming method.
%     STRING may be 'center', 'coarray', 'full', or [LAT LON].  The default
%     is 'center' which is extremely fast for large arrays compared to the
%     other methods as it phase delays the auto-correlation of each record
%     using the distance to the array center (using ARRAYCENTER).  The 
%     'coarray' method utilizes all unique pairings of records to compute
%     the spectrum with a little less than half the cross power spectral
%     density matrix while the 'full' method uses the entire cross power
%     spectral density matrix.  The 'full' method is significantly slower
%     and gives degraded results compared to the 'coarray' method and so is
%     not recommended except in verification (as it gives the exact same
%     result of the 'center' method).  The 'coarray' method has the best
%     resolution out of the available methods.  Using [LAT LON] as a method
%     is algorithmically the same as the 'center' method but uses the
%     defined coordinates as the center for the array.
%
%     [R,T]=FSSHORZ(...,'WHITEN',TRUE|FALSE,...) whitens the spectras
%     before beamforming if WHITEN is TRUE.  The default is TRUE.
%
%     [R,T]=FSSHORZ(...,'WEIGHTS',W,...) specifies the relative weights for
%     each record (must match the size of N & E) or pairing (if METHOD is
%     'coarray' or 'full').
%
%     [R,T]=FSSHORZ(...,'NTILES',NT,...) sets how many nonoverlapping
%     timesections the data should be split into.  The default is 1 which
%     provides the best frequency resolution.
%
%     [R,T]=FSSHORZ(...,'AVG',TRUE|FALSE,...) indicates if the spectra is
%     averaged across frequency during computation.  This can save a
%     significant amount of memory.  The default is false.
%
%    Notes:
%     - Records in N & E must have equal and regular sample spacing.
%     - N & E should contain records in the same order so that N(3) is the
%       North component of the station which has an East component as E(3).
%     - Attenuation is ignored.
%     - References:
%        Rost & Thomas 2002, Array Seismology: Methods and Applications,
%         Rev of Geoph, Vol. 40, No. 3, doi:10.1029/2000RG000100
%        Haubrich & McCamy 1969, Microseisms: Coastal and Pelagic Sources,
%         Rev of Geoph, Vol.  7, No. 3, pp. 539-571
%
%    Examples:
%     % Horizontals need to be rotated to North/East prior to
%     % calling FSSHORZ.  Here is an example of how to do that:
%     data=rotate(data,'to',0,'kcmpnm1','N','kcmpnm2','E');
%     [r,t]=fsshorz(data(1:2:end),data(2:2:end),50,201,[.02 .03],'a',true);
%     plotfss(r);
%     plotfss(t);
%
%     % Compare some artificial horizontal data to ARFHORZ:
%     load fss_test_data;
%     dataW30n=cut(dataW30n,'x',1,'n',500,'fill',true);
%     dataW30e=cut(dataW30e,'x',1,'n',500,'fill',true);
%     [r,t]=fsshorz(dataW30n,dataW30e,50,201,[.02 .03],'a',true);
%     plotfss(r);
%     plotfss(t);
%     [r0,t0]=arfhorz(st,50,201,270,30,r.freq,false,'a',true);
%     plotfss(r0);
%     plotfss(t0);
%
%    See also: FSSHORZXC, ARFHORZ, SNYQUIST, PLOTFSS, KXY2SLOWBAZ, FSS,
%              SLOWBAZ2KXY, FSSAVG, FSSSUB, ARF, FSSXC, FSSDBINFO,
%              FSSFREQSLIDE, FSSFRAMESLIDE, PLOTARF, FSSCORRCOEF

%     Version History:
%        Sep. 22, 2012 - initial version
%        Sep. 27, 2012 - pv pair inputs, doc update, less memory usage,
%                        error for no freq, fixed frequency indexing bug
%        Sep. 30, 2012 - avg & ntiles options
%        Jan.  9, 2013 - allow options to be any case
%        Jan. 14, 2013 - update history
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 14, 2013 at 14:05 GMT

% todo:

% check nargin
error(nargchk(5,inf,nargin));

% check structs
error(seizmocheck(n,'dep'));
error(seizmocheck(e,'dep'));

% require not a xc dataset
if(isxc(e) || isxc(n))
    error('seizmo:fsshorz:badInput',...
        'Use FSSXC for correlations!');
end

% number of records
nrecs=numel(n);
if(nrecs~=numel(e))
    error('seizmo:fsshorz:badInput',...
        'N & E datasets are not the same size!');
end

% need 2+ records
if(nrecs<2)
    error('seizmo:fsshorz:arrayTooSmall',...
        'N & E datasets must have 2+ records!');
end

% check inputs
sf=size(frng);
if(~isreal(smax) || ~isscalar(smax) || smax<=0)
    error('seizmo:fsshorz:badInput',...
        'SMAX must be a positive real scalar in sec/deg!');
elseif(~any(numel(spts)==[1 2]) || any(fix(spts)~=spts) || any(spts<=2))
    error('seizmo:fsshorz:badInput',...
        'SPTS must be a positive scalar integer >2!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<0) ...
        || any(frng(1)>frng(2)))
    error('seizmo:fsshorz:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
end
nrng=sf(1);

% parse options
pv=parse_fsshorz_pv_pairs(varargin{:});

% defaults for optionals
if(isempty(pv.w)); pv.w=ones(nrecs,1); end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    n=checkheader(n,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'NONINTEGER_REFTIME','ERROR',...
        'UNSET_REFTIME','ERROR',...
        'OUTOFRANGE_REFTIME','ERROR',...
        'UNSET_ST_LATLON','ERROR');
    e=checkheader(e,...
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
    [butc,eutc,b,delta,npts,st,kname,cmp]=getheader(n,...
        'b utc','e utc','b','delta','npts','st','kname','cmp');
    [butc2,b2,delta2,npts2,st2,kname2,cmp2]=getheader(e,...
        'b utc','b','delta','npts','st','kname','cmp');
    butc=cell2mat(butc); eutc=cell2mat(eutc); butc2=cell2mat(butc2);
    
    % extract data (silently)
    seizmoverbose(false);
    n=records2mat(n);
    e=records2mat(e);
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

% need n & e to match on several fields:
% kname, st, b, b utc, delta, npts
if(~isequal(butc,butc2) || ~isequal(b,b2) || ~isequal(delta,delta2) ...
        || ~isequal(npts,npts2) || ~isequal(st(:,1:2),st2(:,1:2)) ...
        || ~isequal(kname(:,1:3),kname2(:,1:3)))
    error('seizmo:fsshorz:badInput',...
        'N & E are not synced!');
end
clear butc2 b b2 delta2 npts2 st2 kname kname2;

% require north & east orientation
if(any(cmp(:,1)~=90) || any(cmp2(:,1)~=90)) % cmpinc
    error('seizmo:fsshorz:badInput',...
        'N & E must have CMPINC==90!');
elseif(any(cmp(:,2)~=0))
    error('seizmo:fsshorz:badInput',...
        'N must have CMPAZ==0!');
elseif(any(cmp2(:,2)~=90))
    error('seizmo:fsshorz:badInput',...
        'E must have CMPAZ==90!');
end
clear cmp cmp2;

% fix method/center/npairs
if(ischar(pv.method))
    pv.method=lower(pv.method);
    [clat,clon]=arraycenter(st(:,1),st(:,2));
    switch pv.method
        case 'coarray'
            npairs=nrecs*(nrecs-1)/2;
        case 'full'
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
    error('seizmo:fsshorz:badInput',...
        'Number of WEIGHTS must match the number of stations or pairs!');
end

% get time limits
[bmin,bmini]=min(timediff(butc(1,:),butc,'utc'));
[emax,emaxi]=max(timediff(eutc(1,:),eutc,'utc'));

% check nyquist
fnyq=1/(2*delta(1));
if(any(frng(:,1)>fnyq))
    error('seizmo:fsshorz:badFRNG',...
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

% tiling
if(pv.ntiles>1)
    n=[n; zeros(maxnpts*pv.ntiles-oldmax,nrecs)];
    e=[e; zeros(maxnpts*pv.ntiles-oldmax,nrecs)];
    n=permute(reshape(n.',[nrecs maxnpts pv.ntiles]),[2 1 3]);
    e=permute(reshape(e.',[nrecs maxnpts pv.ntiles]),[2 1 3]);
end

% get fft
n=fft(n,nspts,1);
e=fft(e,nspts,1);

% trim -freq and permute
if(pv.ntiles>1)
    n=permute(n(1+(0:nspts/2),:,:),[2 1 3]);
    e=permute(e(1+(0:nspts/2),:,:),[2 1 3]);
else
    n=n(1+(0:nspts/2),:).';
    e=e(1+(0:nspts/2),:).';
end

% whiten data if desired
if(pv.whiten)
    n=n./sqrt(n.*conj(n)+e.*conj(e));
    e=e./sqrt(n.*conj(n)+e.*conj(e));
end

% setup slowness grid
if(pv.polar)
    if(numel(spts)==1); spts(2)=181; end % default # azimuthal points
    sx=(0:spts(2)-1)/(spts(2)-1)*360; % baz (wedge decided x/y)
    sy=(0:spts(1)-1).'/(spts(1)-1)*smax; % smag
else
    spts(2)=spts(1);
    sx=-smax:2*smax/(spts(1)-1):smax; % east
    sy=fliplr(sx).'; % north
end

% setup output
[r(1:nrng,1).nsta]=deal(nrecs);
[r(1:nrng,1).st]=deal(st);
[r(1:nrng,1).butc]=deal(butc(bmini,:));
[r(1:nrng,1).eutc]=deal(eutc(emaxi,:));
[r(1:nrng,1).delta]=deal(delta(1));
[r(1:nrng,1).npts]=deal(maxnpts);
[r(1:nrng,1).polar]=deal(pv.polar);
[r(1:nrng,1).x]=deal(sx);
[r(1:nrng,1).y]=deal(sy);
[r(1:nrng,1).freq]=deal([]);
[r(1:nrng,1).method]=deal(pv.method);
[r(1:nrng,1).npairs]=deal(npairs);
[r(1:nrng,1).center]=deal([clat clon]);
[r(1:nrng,1).whiten]=deal(pv.whiten);
[r(1:nrng,1).weights]=deal(pv.w);
[r(1:nrng,1).spectra]=deal(zeros(0,'single'));
[r(1:nrng,1).orient]=deal('radial');
t=r;
[t(1:nrng,1).orient]=deal('transverse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of defaults, checks, & data setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% radius constants
d2r=pi/180;
d2km=6371*d2r;

% get number of slowness points
nslow=prod(spts);

% get utc static time shifts
dt=timediff(butc(1,:),butc,'utc').'; % need row vector
switch pv.method
    case 'coarray'
        % row is master index, column is slave index
        [master,slave]=find(triu(true(nrecs),1));
        dt=dt(slave)-dt(master); % alter to pair shifts
    case 'full'
        [master,slave]=find(true(nrecs));
        dt=dt(slave)-dt(master); % alter to pair shifts
end
dt=dt(ones(nslow,1),:);

% normalize & expand weights
pv.w=pv.w(:);
switch pv.method
    case {'coarray' 'full'}
        if(numel(pv.w)==nrecs); pv.w=pv.w(slave).*conj(pv.w(master)); end
end
pv.w=pv.w./sum(abs(pv.w));

% project locations onto tangent plane (tangent to earth at array center)
[ekm,nkm]=geographic2enu(st(:,1),st(:,2),0,clat,clon,0);

% get the angular difference between true and
% grid (aka tangent plane) north for each station
% - THIS IS SIGNIFICANT AT HIGH LATITUDES
dkm=max(abs([ekm;nkm])); % maximum north/east deviation from center
dn=dkm/d2km/1000;        % how far north from station to get difference
[ekm1,nkm1]=geographic2enu(st(:,1)+dn,st(:,2),0,clat,clon,0);
truenorth=atan2(ekm1-ekm,nkm1-nkm).'/d2r; % 1xNRECS

% get relative positions for each pair (the coarray)
% r=(x  ,y  )
%     ij  ij
%
% x   is km east, y   is km north
%  ij              ij
%
% r is a 2xNPAIRS matrix
switch pv.method
    case {'coarray' 'full'}
        % [ r   r   ... r
        %    11  12      1N
        %   r   r   ... r
        %    21  22      2N
        %    .   .  .    .
        %    .   .   .   .
        %    .   .    .  .
        %   r   r   ... r   ]
        %    N1  N2      NN
        ekm=ekm(slave)-ekm(master);
        nkm=nkm(slave)-nkm(master);
end
rkm=[ekm(:) nkm(:)]';
clear ekm nkm

% Get rotators for each planewave and station.
% These can be used to rotate from true North & East (which is not North &
% East on the tangent plane) to each planewave's radial & transverse
% directions.
%
% u=cos(baz+180-truenorth)=-cos(baz-truenorth)
% v=sin(baz+180-truenorth)=-sin(baz-truenorth)
%
% where
%  baz is the backazimuth of the planewave
%  truenorth is the azimuth of true north from the tangent plane north
%
% u & v are NSLOWxNRECS
if(pv.polar); baz=sx(ones(spts(1),1),:);
else baz=atan2(sx(ones(spts(1),1),:),sy(:,ones(spts(1),1)))/d2r;
end
baz=baz(:);
u=cosd(baz(:,ones(nrecs,1))-truenorth(ones(nslow,1),:));
v=sind(baz(:,ones(nrecs,1))-truenorth(ones(nslow,1),:));

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
    sy=sy(:,ones(spts(2),1))/d2km; % smag in sec/km
    [sx,sy]=deal(sy.*sind(sx),sy.*cosd(sx));
    p=[sx(:) sy(:)]*rkm;
else % cartesian
    sx=sx/d2km;
    sx=sx(ones(spts(1),1),:);
    sy=fliplr(sx)';
    p=[sx(:) sy(:)]*rkm;
end
clear rkm sx sy

% loop over frequency ranges
for a=1:nrng
    % get frequencies
    fidx=find(f>=frng(a,1) & f<=frng(a,2));
    r(a).freq=f(fidx);
    t(a).freq=f(fidx);
    nfreq=numel(fidx);
    
    % preallocate spectra
    if(pv.avg)
        r(a).spectra=zeros(spts,'single');
        t(a).spectra=zeros(spts,'single');
    else
        r(a).spectra=zeros([spts nfreq],'single');
        t(a).spectra=zeros([spts nfreq],'single');
    end
    
    % error if no frequencies
    if(~nfreq)
        error('seizmo:fsshorz:noFreqs',...
            'No data within the range %g to %g Hz!',...
            frng(a,1),frng(a,2));
    end
    
    % detail message
    if(verbose)
        fprintf('Getting horizontal spectra %d: %g to %g Hz\n',...
            a,frng(a,1),frng(a,2));
        print_time_left(0,nfreq);
    end
    
    % get spectra based upon method
    switch pv.method
        case {'coarray' 'full'}
            if(pv.avg)
                for b=1:nfreq
                    r(a).spectra=r(a).spectra+reshape(real(...
                         (u(:,master).*u(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(n(master,fidx(b),:)).*n(slave,fidx(b),:),3).*pv.w) ...
                        +(u(:,master).*v(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(n(master,fidx(b),:)).*e(slave,fidx(b),:),3).*pv.w) ...
                        +(v(:,master).*u(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(e(master,fidx(b),:)).*n(slave,fidx(b),:),3).*pv.w) ...
                        +(v(:,master).*v(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(e(master,fidx(b),:)).*e(slave,fidx(b),:),3).*pv.w)),spts);
                    t(a).spectra=t(a).spectra+reshape(real(...
                         (v(:,master).*v(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(n(master,fidx(b),:)).*n(slave,fidx(b),:),3).*pv.w) ...
                        -(v(:,master).*u(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(n(master,fidx(b),:)).*e(slave,fidx(b),:),3).*pv.w) ...
                        -(u(:,master).*v(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(e(master,fidx(b),:)).*n(slave,fidx(b),:),3).*pv.w) ...
                        +(u(:,master).*u(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(e(master,fidx(b),:)).*e(slave,fidx(b),:),3).*pv.w)),spts);
                    if(verbose); print_time_left(b,nfreq); end
                end
            else
                for b=1:nfreq
                    r(a).spectra(:,:,b)=reshape(real(...
                         (u(:,master).*u(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(n(master,fidx(b),:)).*n(slave,fidx(b),:),3).*pv.w) ...
                        +(u(:,master).*v(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(n(master,fidx(b),:)).*e(slave,fidx(b),:),3).*pv.w) ...
                        +(v(:,master).*u(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(e(master,fidx(b),:)).*n(slave,fidx(b),:),3).*pv.w) ...
                        +(v(:,master).*v(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(e(master,fidx(b),:)).*e(slave,fidx(b),:),3).*pv.w)),spts);
                    t(a).spectra(:,:,b)=reshape(real(...
                         (v(:,master).*v(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(n(master,fidx(b),:)).*n(slave,fidx(b),:),3).*pv.w) ...
                        -(v(:,master).*u(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(n(master,fidx(b),:)).*e(slave,fidx(b),:),3).*pv.w) ...
                        -(u(:,master).*v(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(e(master,fidx(b),:)).*n(slave,fidx(b),:),3).*pv.w) ...
                        +(u(:,master).*u(:,slave).*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(mean(conj(e(master,fidx(b),:)).*e(slave,fidx(b),:),3).*pv.w)),spts);
                    if(verbose); print_time_left(b,nfreq); end
                end
            end
        otherwise % {'center' 'user'}
            if(pv.avg)
                for b=1:nfreq
                    for c=1:pv.ntiles
                        r(a).spectra=r(a).spectra+reshape(abs(...
                             (u.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(n(:,fidx(b),c).*pv.w) ...
                            +(v.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(e(:,fidx(b),c).*pv.w)).^2,spts);
                        t(a).spectra=t(a).spectra+reshape(abs(...
                            -(v.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(n(:,fidx(b),c).*pv.w) ...
                            +(u.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(e(:,fidx(b),c).*pv.w)).^2,spts);
                    end
                    if(verbose); print_time_left(b,nfreq); end
                end
            else
                for b=1:nfreq
                    for c=1:pv.ntiles
                        r(a).spectra(:,:,b)=r(a).spectra(:,:,b)+reshape(abs(...
                             (u.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(n(:,fidx(b),c).*pv.w) ...
                            +(v.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(e(:,fidx(b),c).*pv.w)).^2,spts);
                        t(a).spectra(:,:,b)=t(a).spectra(:,:,b)+reshape(abs(...
                            -(v.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(n(:,fidx(b),c).*pv.w) ...
                            +(u.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(e(:,fidx(b),c).*pv.w)).^2,spts);
                    end
                    if(verbose); print_time_left(b,nfreq); end
                end
            end
            r(a).spectra=r(a).spectra/pv.ntiles;
            t(a).spectra=t(a).spectra/pv.ntiles;
    end
    if(pv.avg)
        r(a).spectra=r(a).spectra/nfreq;
        t(a).spectra=t(a).spectra/nfreq;
    end
end

end


function [pv]=parse_fsshorz_pv_pairs(varargin)

% defaults
pv.avg=false;
pv.polar=false;
pv.method='center';
pv.whiten=true;
pv.w=[];
pv.ntiles=1;

% require pv pairs
if(mod(nargin,2))
    error('seizmo:fsshorz:badInput',...
        'Unpaired parameter/value!');
end

% parse pv pairs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:fsshorz:badInput',...
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
        case {'ntiles' 'ntile' 'tiles' 'tile' 'nt' 't'}
            pv.ntiles=varargin{i+1};
        case {'average' 'av' 'avg' 'a'}
            pv.avg=varargin{i+1};
        otherwise
            error('seizmo:fsshorz:badInput',...
                'Unknown parameter: %s !',varargin{i});
    end
end

% valid method strings
valid.METHOD={'center' 'coarray' 'full'};

% check values
if(~isscalar(pv.polar) || (~islogical(pv.polar) && ~isnumeric(pv.polar)))
    error('seizmo:fsshorz:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif((isnumeric(pv.method) ...
        && (~isreal(pv.method) || ~numel(pv.method)==2)) ...
        || (ischar(pv.method) && ~any(strcmpi(pv.method,valid.METHOD))))
    error('seizmo:fsshorz:badInput',...
        ['METHOD must be one of the following:\n' ...
        '''CENTER'', ''COARRAY'', ''FULL'', or [LAT LON]!']);
elseif(~isscalar(pv.whiten) || ~islogical(pv.whiten))
    error('seizmo:fsshorz:badInput',...
        'WHITEN must be TRUE or FALSE!');
elseif(~isnumeric(pv.w) || any(pv.w(:)<0))
    error('seizmo:fsshorz:badInput',...
        'WEIGHTS must be positive!');
elseif(~isreal(pv.ntiles) || ~isscalar(pv.ntiles) ...
        || pv.ntiles<=0 || pv.ntiles~=fix(pv.ntiles))
    error('seizmo:fsshorz:badInput',...
        'NTILES must be a positive integer!');
elseif(~isscalar(pv.avg) || ~islogical(pv.avg))
    error('seizmo:fsshorz:badInput',...
        'AVG must be TRUE or FALSE!');
end

end


function [lgc]=isxc(data)
[m,s]=getheader(data,'kuser0','kuser1');
if(all(strcmp(m,'MASTER') & strcmp(s,'SLAVE'))); lgc=true;
else lgc=false;
end
end

