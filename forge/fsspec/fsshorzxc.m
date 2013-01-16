function [r,t]=fsshorzxc(nn,ne,en,ee,smax,spts,frng,varargin)
%FSSHORZXC    Estimate frequency-slowness spectrum of horizontals from xc
%
%    Usage:   [r,t]=fsshorzxc(nn,ne,en,ee,smax,spts,frng)
%             [r,t]=fsshorzxc(rr,rt,tr,tt,smax,spts,frng)
%             [r,t]=fsshorzxc(...,'polar',true|false,...)
%             [r,t]=fsshorzxc(...,'whiten',true|false,...)
%             [r,t]=fsshorzxc(...,'weights',w,...)
%             [r,t]=fsshorzxc(...,'avg',true|false,...)
%
%    Description:
%     [R,T]=FSSHORZXC(NN,NE,EN,EE,SMAX,SPTS,FRNG) computes an estimate of
%     the radial and tranverse frequency-slowness power spectra in a
%     cartesian grid for a spatial array of stations recording horizontal
%     motions by frequency-domain beamforming the cross correlations in NN,
%     NE, EN & EE.  The datasets NN, NE, EN & EE are SEIZMO structs
%     containing the array info and correlations of motion in the North and
%     East directions.  This function differs from GEOFSSHORZXC in that the
%     waves are defined as plane waves traveling on a planar surface rather
%     than surface waves expanding and contracting on a sphere.  The range
%     of the horizontal slowness grid is given by SMAX (sec/deg) and
%     extends from -SMAX to SMAX for both East/West and North/South
%     directions.  SPTS controls the number of slowness points for both
%     directions (SPTSxSPTS grid).  FRNG gives the frequency range as
%     [FREQLOW FREQHIGH] in Hz (the individual frequencies are determined
%     by the fft of the data).  The outputs R & T are structs containing
%     relevant info and the frequency-slowness spectra for the radial and
%     transverse wavefields (with size SPTSxSPTSxNFREQ).  The struct layout
%     is:
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
%          .method   - beamforming method ('center', 'coarray', or 'user')
%          .center   - array center as [LAT LON]
%          .whiten   - is spectra whitened? true/false
%          .weights  - weights used in beamforming
%          .spectra  - frequency-slowness spectra estimate
%          .orient   - 'radial' or 'transverse'
%
%     [R,T]=FSSHORZXC(RR,RT,TR,TT,SMAX,SPTS,FRNG) provides a convenient
%     alternative for inputing horizontal correlograms.  The correlograms
%     are expected to be in pairwise radial and transverse directions.
%     ROTATE_CORRELATIONS gives correlograms in this orientation (and takes
%     them in the E/N system).
%
%     [R,T]=FSSHORZXC(...,'POLAR',TRUE|FALSE,...) toggles if the slownesses
%     are sampled regularly in cartesian or polar coordinates.  Polar
%     coords are useful for slicing the spectra by azimuth (pie slice) or
%     slowness (rings).  Cartesian coords (the default) samples the
%     slowness space regularly in the East/West & North/South directions
%     and so exhibits less distortion in plots of the slowness space. If
%     POLAR=TRUE, SPTS may be given as [SPTS BAZPTS] to control the
%     azimuthal resolution (default is BAZPTS=181 points).
%
%     [R,T]=FSSHORZXC(...,'WHITEN',TRUE|FALSE,...) whitens the spectras
%     before beamforming if WHITEN is TRUE.  The default is TRUE.
%
%     [R,T]=FSSHORZXC(...,'WEIGHTS',W,...) specifies the relative weights
%     for each record (must match dataset size).  That is WEIGHTS, NN, NE,
%     EN, & EE must be the same size.
%
%     [R,T]=FSSHORZXC(...,'AVG',TRUE|FALSE,...) indicates if the spectra is
%     averaged across frequency during computation.  This can save a
%     significant amount of memory.  The default is false.
%
%    Notes:
%     - Correlations in NN, NE, EN & EE must have equal and regular sample
%       spacing.
%     - Correlations are expected to have station naming and location info
%       as stored in the header by the function CORRELATE.  This means the
%       headers of correlations done outside of SEIZMO will probably
%       require adjustment.
%     - A specific record of NN, NE, EN, & EE should correspond to the same
%       pair of stations.  So NN(3), NE(3), EN(3) & EE(3) should correspond
%       to the correlation of horizontals from the same two stations.
%     - Best/quickest results are obtained when using only correlations
%       of unique station pairs (order independent).
%     - Attenuation is ignored.
%     - References:
%        Rost & Thomas 2002, Array Seismology: Methods and Applications,
%         Rev of Geoph, Vol. 40, No. 3, doi:10.1029/2000RG000100
%        Haubrich & McCamy 1969, Microseisms: Coastal and Pelagic Sources,
%         Rev of Geoph, Vol.  7, No. 3, pp. 539-571
%
%    Examples:
%     % Horizontals need to be rotated to North/East and then correlated
%     % prior to calling FSSHORZXC.  Here is an example of how to do that:
%     data=rotate(data,'to',0,'kcmpnm1','N','kcmpnm2','E');
%     xc=correlate(data,'normxc',false);
%     [m,s]=getheader(xc,'user0','user1'); % retreive indexing
%     xc(mod(m,2) & m+1==s)=[]; % remove correlations for same station
%     [m,s]=getheader(xc,'user0','user1'); % retreive new indexing
%     nn=xc( mod(m,2) &  mod(s,2));
%     ne=xc( mod(m,2) & ~mod(s,2));
%     en=xc(~mod(m,2) &  mod(s,2));
%     ee=xc(~mod(m,2) & ~mod(s,2));
%     [r,t]=fsshorzxc(nn,ne,en,ee,50,201,[1/50 1/40]);
%     plotfss(fssavg(r));
%     plotfss(fssavg(t));
%
%    See also: FSSHORZ, ARFHORZ, SNYQUIST, PLOTFSS, KXY2SLOWBAZ, FSS,
%              SLOWBAZ2KXY, FSSAVG, FSSSUB, ARF, FSSXC, FSSDBINFO,
%              FSSFREQSLIDE, FSSFRAMESLIDE, PLOTARF, FSSCORRCOEF

%     Version History:
%        Sep. 22, 2012 - initial version
%        Sep. 27, 2012 - pv pair inputs, doc update, less memory usage,
%                        error for no freq, fixed frequency indexing bug
%        Sep. 30, 2012 - avg option
%        Jan.  9, 2013 - allow options to be any case
%        Jan. 14, 2013 - update history
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 14, 2013 at 14:05 GMT

% todo:

% check nargin
error(nargchk(7,inf,nargin));

% check structs
error(seizmocheck(nn,'dep'));
error(seizmocheck(ne,'dep'));
error(seizmocheck(en,'dep'));
error(seizmocheck(ee,'dep'));

% require xc dataset
if(~isxc(nn) || ~isxc(ne) || ~isxc(en) || ~isxc(ee))
    error('seizmo:fsshorzxc:badInput',...
        'NN, NE, EN & EE must be correlations!');
end

% number of pairings
npairs=numel(nn);
if(~isequal(npairs,numel(ne),numel(en),numel(ee)))
    error('seizmo:fsshorzxc:badInput',...
        'NN, NE, EN & EE datasets are not the same size!');
end

% check inputs
sf=size(frng);
if(~isreal(smax) || ~isscalar(smax) || smax<=0)
    error('seizmo:fsshorzxc:badInput',...
        'SMAX must be a positive real scalar in sec/deg!');
elseif(~any(numel(spts)==[1 2]) || any(fix(spts)~=spts) || any(spts<=2))
    error('seizmo:fsshorzxc:badInput',...
        'SPTS must be a positive scalar integer >2!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<0) ...
        || any(frng(1)>frng(2)))
    error('seizmo:fsshorzxc:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
end
nrng=sf(1);

% parse options
pv=parse_fsshorzxc_pv_pairs(varargin{:});

% defaults for optionals
if(isempty(pv.w)); pv.w=ones(npairs,1); end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    nn=checkheader(nn,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    ne=checkheader(ne,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    en=checkheader(en,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    ee=checkheader(ee,...
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

% extract header info & data
try
    % verbosity
    verbose=seizmoverbose;
    
    % grab necessary header info
    [b,npts,delta,autc,futc,st,ev,stnm,evnm,scmp,ecmp]=getheader(nn,'b',...
       'npts','delta','a utc','f utc','st','ev','kname','kt','cmp','user');
    [b2,npts2,delta2,st2,ev2,stnm2,evnm2,scmp2,ecmp2]=getheader(ne,'b',...
       'npts','delta','st','ev','kname','kt','cmp','user');
    [b3,npts3,delta3,st3,ev3,stnm3,evnm3,scmp3,ecmp3]=getheader(en,'b',...
       'npts','delta','st','ev','kname','kt','cmp','user');
    [b4,npts4,delta4,st4,ev4,stnm4,evnm4,scmp4,ecmp4]=getheader(ee,'b',...
       'npts','delta','st','ev','kname','kt','cmp','user');
    autc=cell2mat(autc); futc=cell2mat(futc);
    ecmp=ecmp(:,3:4); ecmp2=ecmp2(:,3:4);
    ecmp3=ecmp3(:,3:4); ecmp4=ecmp4(:,3:4);
    
    % for rt case
    [az,gcp]=getheader(nn,'az','gcp');
    
    % extract data (silently)
    seizmoverbose(false);
    nn=records2mat(nn);
    ne=records2mat(ne);
    en=records2mat(en);
    ee=records2mat(ee);
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

% need nn, ne, en & ee to match on several fields
if(~isequal(b,b2,b3,b4) ...
        || ~isequal(delta,delta2,delta3,delta4) ...
        || ~isequal(npts,npts2,npts3,npts4) ...
        || ~isequal(st(:,1:2),st2(:,1:2),st3(:,1:2),st4(:,1:2)) ...
        || ~isequal(ev(:,1:2),ev2(:,1:2),ev3(:,1:2),ev4(:,1:2)) ...
        || ~isequal(stnm(:,1:3),stnm2(:,1:3),stnm3(:,1:3),stnm4(:,1:3)) ...
        || ~isequal(evnm(:,1:3),evnm2(:,1:3),evnm3(:,1:3),evnm4(:,1:3)))
    error('seizmo:fsshorzxc:badInput',...
        'NN, NE, EN & EE are not synced!');
end
clear b2 b3 b4 delta2 delta3 delta4 npts2 npts3 npts4;
clear st2 st3 st4 ev2 ev3 ev4 stnm2 stnm3 stnm4 evnm2 evnm3 evnm4;

% require north & east orientation or pairwise radial & transverse
if(any([scmp(:,1);scmp2(:,1);scmp3(:,1);scmp4(:,1)]~=90) ...
        || any([ecmp(:,1);ecmp2(:,1);ecmp3(:,1);ecmp4(:,1)]~=90))
    error('seizmo:fsshorzxc:badInput',...
        'All datasets must have CMPINC==90!');
end
try
    isne=true;
    if(any([scmp(:,2);ecmp(:,2)]~=0))
        error('seizmo:fsshorzxc:badInput',...
            'NN must have CMPAZ==0 & USER3==0!');
    elseif(any(scmp2(:,2)~=90) || any(ecmp2(:,2)~=0))
        error('seizmo:fsshorzxc:badInput',...
            'NE must have CMPAZ==90 & USER3==0!');
    elseif(any(scmp3(:,2)~=0) || any(ecmp3(:,2)~=90))
        error('seizmo:fsshorzxc:badInput',...
            'EN must have CMPAZ==0 & USER3==90!');
    elseif(any([scmp4(:,2);ecmp4(:,2)]~=90))
        error('seizmo:fsshorzxc:badInput',...
            'EE must have CMPAZ==90 & USER3==90!');
    end
catch
    isne=false; % RT?
    if(any(abs(azdiff(scmp(:,2),gcp))>1) ...
            || any(abs(azdiff(ecmp(:,2),az))>1))
        error('seizmo:fsshorzxc:badInput',...
            'RR must have CMPAZ==GCP & USER3==AZ!');
    elseif(any(abs(azdiff(scmp2(:,2),gcp+90))>1) ...
            || any(abs(azdiff(ecmp2(:,2),az))>1))
        error('seizmo:fsshorzxc:badInput',...
            'RT must have CMPAZ==GCP+90 & USER3==AZ!');
    elseif(any(abs(azdiff(scmp3(:,2),gcp))>1) ...
            || any(abs(azdiff(ecmp3(:,2),az+90))>1))
        error('seizmo:fsshorzxc:badInput',...
            'TR must have CMPAZ==GCP & USER3==AZ+90!');
    elseif(any(abs(azdiff(scmp4(:,2),gcp+90))>1) ...
            || any(abs(azdiff(ecmp4(:,2),az+90))>1))
        error('seizmo:fsshorzxc:badInput',...
            'TT must have CMPAZ==GCP+90 & USER3==AZ+90!');
    end
end
clear scmp2 scmp3 scmp4 ecmp2 ecmp3 ecmp4;

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
    error('seizmo:fsshorzxc:badFRNG',...
        ['FRNG exceeds nyquist frequency (' num2str(fnyq) ')!']);
end

% longest record
maxnpts=max(npts);

% get frequencies
nspts=2^nextpow2(maxnpts);
f=(0:nspts/2)/(delta(1)*nspts);  % only +freq

% get fft
nn=fft(nn,nspts,1);
ne=fft(ne,nspts,1);
en=fft(en,nspts,1);
ee=fft(ee,nspts,1);
nn=nn(1+(0:nspts/2),:).'; % trim -freq
ne=ne(1+(0:nspts/2),:).'; % trim -freq
en=en(1+(0:nspts/2),:).'; % trim -freq
ee=ee(1+(0:nspts/2),:).'; % trim -freq

% whiten data if desired
if(pv.whiten)
    tmp=sqrt(nn.*conj(nn)+ne.*conj(ne)+en.*conj(en)+ee.*conj(ee));
    nn=nn./tmp; ne=ne./tmp; en=en./tmp; ee=ee./tmp;
    clear tmp;
end

% use unique stations names to get number of stations & locations
evnm=evnm(:,1:4);
stnm=strcat(stnm(:,1),'.',stnm(:,2),'.',stnm(:,3),'.',stnm(:,4));
evnm=strcat(evnm(:,1),'.',evnm(:,2),'.',evnm(:,3),'.',evnm(:,4));
[stnm,idx1,idx2]=unique([stnm;evnm]);
st=[st;ev];
st=st(idx1,:);
nrecs=numel(stnm);

% check weights again
if(~any(numel(pv.w)==[nrecs npairs]))
    error('seizmo:fsshorzxc:badInput',...
        '# of WEIGHTS must match the # of records in N & E or pairs!');
end

% slave/master indices
idx2=reshape(idx2,[],2);
slave=idx2(:,1);
master=idx2(:,2);

% array center
[clat,clon]=arraycenter(st(:,1),st(:,2));

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
[r(1:nrng,1).butc]=deal(autc);
[r(1:nrng,1).eutc]=deal(futc);
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

% get static time shifts
% - MUST shift to zero lag
dt=b(:,ones(nslow,1)).';

% normalize & expand weights
pv.w=pv.w(:);
if(numel(pv.w)==nrecs); pv.w=pv.w(slave).*conj(pv.w(master)); end
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
ekm=ekm(slave)-ekm(master);
nkm=nkm(slave)-nkm(master);
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
% Rotation from pairwise radial & transverse to each planewave's radial &
% transverse directions is complicated by orientations dependent on the
% pairing.  So we have to expand to pairs.
%
% us=cos(baz+180-truenorth-gcp)=-cos(baz-truenorth-gcp)
% um=cos(baz+180-truenorth-az)=-cos(baz-truenorth-az)
% vs=sin(baz+180-truenorth-gcp)=-sin(baz-truenorth-gcp)
% vm=sin(baz+180-truenorth-az)=-sin(baz-truenorth-az)
%
% where
%  gcp is the direction to the master station + 180 deg
%  az is the direction to the slave station
%
% NOTE WE DROP THE SIGNS B/C THEY CANCEL OUT
%
% u & v are NSLOWxNPAIRS
if(pv.polar); baz=sx(ones(spts(1),1),:);
else baz=atan2(sx(ones(spts(1),1),:),sy(:,ones(spts(1),1)))/d2r;
end
baz=baz(:);
if(isne)
    us=cosd(baz(:,ones(nrecs,1))-truenorth(ones(nslow,1),:));
    vs=sind(baz(:,ones(nrecs,1))-truenorth(ones(nslow,1),:));
    um=us(:,master); vm=vs(:,master);
    us=us(:,slave); vs=vs(:,slave);
else % rt
    us=cosd(baz(:,ones(npairs,1))-truenorth(ones(nslow,1),slave)...
        -gcp(:,ones(nslow,1)).');
    vs=sind(baz(:,ones(npairs,1))-truenorth(ones(nslow,1),slave)...
        -gcp(:,ones(nslow,1)).');
    um=cosd(baz(:,ones(npairs,1))-truenorth(ones(nslow,1),master)...
        -az(:,ones(nslow,1)).');
    vm=sind(baz(:,ones(npairs,1))-truenorth(ones(nslow,1),master)...
        -az(:,ones(nslow,1)).');
end

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
        error('seizmo:fsshorzxc:noFreqs',...
            'No data within the range %g to %g Hz!',...
            frng(a,1),frng(a,2));
    end
    
    % detail message
    if(verbose)
        fprintf('Getting horizontal spectra %d: %g to %g Hz\n',...
            a,frng(a,1),frng(a,2));
        print_time_left(0,nfreq);
    end
    
    % get spectra
    if(pv.avg)
        for b=1:nfreq
            r(a).spectra=r(a).spectra+reshape(real(...
                 (um.*us.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(nn(:,fidx(b)).*pv.w) ...
                +(um.*vs.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(ne(:,fidx(b)).*pv.w) ...
                +(vm.*us.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(en(:,fidx(b)).*pv.w) ...
                +(vm.*vs.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(ee(:,fidx(b)).*pv.w)),spts);
            t(a).spectra=t(a).spectra+reshape(real(...
                 (vm.*vs.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(nn(:,fidx(b)).*pv.w) ...
                -(vm.*us.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(ne(:,fidx(b)).*pv.w) ...
                -(um.*vs.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(en(:,fidx(b)).*pv.w) ...
                +(um.*us.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(ee(:,fidx(b)).*pv.w)),spts);
            if(verbose); print_time_left(b,nfreq); end
        end
    else
        for b=1:nfreq
            r(a).spectra(:,:,b)=reshape(real(...
                 (um.*us.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(nn(:,fidx(b)).*pv.w) ...
                +(um.*vs.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(ne(:,fidx(b)).*pv.w) ...
                +(vm.*us.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(en(:,fidx(b)).*pv.w) ...
                +(vm.*vs.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(ee(:,fidx(b)).*pv.w)),spts);
            t(a).spectra(:,:,b)=reshape(real(...
                 (vm.*vs.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(nn(:,fidx(b)).*pv.w) ...
                -(vm.*us.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(ne(:,fidx(b)).*pv.w) ...
                -(um.*vs.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(en(:,fidx(b)).*pv.w) ...
                +(um.*us.*exp(-2*pi*1i*f(fidx(b))*(p+dt)))*(ee(:,fidx(b)).*pv.w)),spts);
            if(verbose); print_time_left(b,nfreq); end
        end
    end
    if(pv.avg)
        r(a).spectra=r(a).spectra/nfreq;
        t(a).spectra=t(a).spectra/nfreq;
    end
end

end


function [pv]=parse_fsshorzxc_pv_pairs(varargin)

% defaults
pv.avg=false;
pv.polar=false;
pv.method='xc';
pv.whiten=true;
pv.w=[];

% require pv pairs
if(mod(nargin,2))
    error('seizmo:fsshorzxc:badInput',...
        'Unpaired parameter/value!');
end

% parse pv pairs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:fsshorzxc:badInput',...
        'Parameters must be specified using a string!');
end
for i=1:2:nargin
    switch lower(varargin{i})
        case {'polar' 'plr' 'pol' 'p'}
            pv.polar=varargin{i+1};
        case {'whiten' 'wh' 'white' 'normalize' 'n' 'coherency' 'c'}
            pv.whiten=varargin{i+1};
        case {'weights' 'w' 'wgt' 'wgts' 'weight'}
            pv.w=varargin{i+1};
        case {'average' 'av' 'avg' 'a'}
            pv.avg=varargin{i+1};
        otherwise
            error('seizmo:fsshorzxc:badInput',...
                'Unknown parameter: %s !',varargin{i});
    end
end

% check values
if(~isscalar(pv.polar) || (~islogical(pv.polar) && ~isnumeric(pv.polar)))
    error('seizmo:fsshorzxc:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif(~isscalar(pv.whiten) || ~islogical(pv.whiten))
    error('seizmo:fsshorzxc:badInput',...
        'WHITEN must be TRUE or FALSE!');
elseif(~isnumeric(pv.w) || any(pv.w(:)<0))
    error('seizmo:fsshorzxc:badInput',...
        'WEIGHTS must be positive!');
elseif(~isscalar(pv.avg) || ~islogical(pv.avg))
    error('seizmo:fsshorzxc:badInput',...
        'AVG must be TRUE or FALSE!');
end

end


function [lgc]=isxc(data)
[m,s]=getheader(data,'kuser0','kuser1');
if(all(strcmp(m,'MASTER') & strcmp(s,'SLAVE'))); lgc=true;
else lgc=false;
end
end

