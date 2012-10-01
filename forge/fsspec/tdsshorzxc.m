function [r,t]=tdsshorzxc(nn,ne,en,ee,smax,spts,varargin)
%TDSSHORZXC    Estimate slowness spectrum of horizontals from xc
%
%    Usage:   [r,t]=tdsshorzxc(nn,ne,en,ee,smax,spts)
%             [r,t]=tdsshorzxc(rr,rt,tr,tt,smax,spts)
%             [r,t]=tdsshorzxc(...,'polar',true|false,...)
%             [r,t]=tdsshorzxc(...,'weights',w,...)
%
%    Description:
%     [R,T]=TDSSHORZXC(NN,NE,EN,EE,SMAX,SPTS) computes an estimate of
%     the radial and tranverse slowness power spectra in a cartesian grid
%     for a spatial array of stations recording horizontal motions by time-
%     domain beamforming the cross correlations in NN, NE, EN & EE.  The
%     datasets NN, NE, EN & EE are SEIZMO structs containing the array info
%     and correlations of motion in the North and East directions.  This
%     function differs from GEOTDSSHORZXC in that the waves are defined as
%     plane waves traveling on a planar surface rather than surface waves
%     expanding and contracting on a sphere.  The range of the horizontal
%     slowness grid is given by SMAX (sec/deg) and extends from -SMAX to
%     SMAX for both East/West and North/South directions.  SPTS controls
%     the number of slowness points for both directions (SPTSxSPTS grid).
%     The outputs R & T are structs containing relevant info and the
%     slowness spectra for the radial and transverse wavefields (with size
%     SPTSxSPTS).  The struct layout is:
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
%          .method   - beamforming method ('tdxc')
%          .center   - array center as [LAT LON]
%          .whiten   - is spectra whitened? (false)
%          .weights  - weights used in beamforming
%          .spectra  - slowness spectra estimate
%          .orient   - 'radial' or 'transverse'
%
%     [R,T]=TDSSHORZXC(RR,RT,TR,TT,SMAX,SPTS) provides a convenient
%     alternative for inputing horizontal correlograms.  The correlograms
%     are expected to be in pairwise radial and transverse directions.
%     ROTATE_CORRELATIONS gives correlograms in this orientation (and takes
%     them in the E/N system).
%
%     [R,T]=TDSSHORZXC(...,'POLAR',TRUE|FALSE,...) sets if the slownesses
%     are sampled regularly in cartesian or polar coordinates.  Polar
%     coords are useful for slicing the spectra by azimuth (pie slice) or
%     slowness (rings).  Cartesian coords (the default) samples the
%     slowness space regularly in the East/West & North/South directions
%     and so exhibits less distortion in plots of the slowness space. If
%     POLAR=TRUE, SPTS may be given as [SPTS BAZPTS] to control the
%     azimuthal resolution (default is BAZPTS=181 points).
%
%     [R,T]=TDSSHORZXC(...,'WEIGHTS',W,...) specifies the relative weights
%     for each record (must match dataset size).  That is WEIGHTS, NN, NE,
%     EN, & EE must be the same size.
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
%     - USE 2-PASS FILTERS ON CORRELATIONS SO THE PHASE INFORMATION IS
%       PRESERVED!  SINGLE PASS FILTERS WILL GIVE INCORRECT RESULTS!
%     - References:
%        Brzak et al 2009, Migration imaging and forward modeling of
%         microseismic noise sources near southern Italy, G3, Vol. 10,
%         No. 1, Q01012, doi:10.1029/2008GC002234
%        Burtin et al 2010, Location of river-induced seismic signal from
%         noise correlation functions, GJI, Vol. 182, pp. 1161-1173,
%         doi:10.1111/j.1365-246X.2010.04701.x
%        Gu & Shen 2012, Microseismic Noise from Large Ice-Covered Lakes?,
%         BSSA, Vol. 102, No. 3, pp. 1155-1166, doi:101785/0120100010
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
%     [r,t]=tdsshorzxc(nn,ne,en,ee,50,201);
%     plotfss(fssavg(r));
%     plotfss(fssavg(t));
%
%     % Filter to 20-25s and compare to FSSHORZXC output:
%     tic;
%     [rf,tf]=fsshorzxc(nn,ne,en,ee,50,201,[1/25 1/20],'avg',true);
%     toc; tic;
%     nn=iirfilter(nn,'bp','b','o',4,'c',[1/25 1/20],'p',2);
%     ne=iirfilter(ne,'bp','b','o',4,'c',[1/25 1/20],'p',2);
%     en=iirfilter(en,'bp','b','o',4,'c',[1/25 1/20],'p',2);
%     ee=iirfilter(ee,'bp','b','o',4,'c',[1/25 1/20],'p',2);
%     [rt,tt]=tdsshorzxc(nn,ne,en,ee,50,201);
%     toc;
%     plotfss(rf); plotfss(tf);
%     plotfss(rt); plotfss(tt);
%
%    See also: TDSSXC, GEOTDSSXC, GEOTDSSHORZXC, FSSXC, GEOFSSXC,
%              ARF, ARFHORZ, GEOARF, GEOARFHORZ

%     Version History:
%        Sep. 30, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 30, 2012 at 14:05 GMT

% todo:

% check nargin
error(nargchk(6,inf,nargin));

% check structs
error(seizmocheck(nn,'dep'));
error(seizmocheck(ne,'dep'));
error(seizmocheck(en,'dep'));
error(seizmocheck(ee,'dep'));

% require xc dataset
if(~isxc(nn) || ~isxc(ne) || ~isxc(en) || ~isxc(ee))
    error('seizmo:tdsshorzxc:badInput',...
        'NN, NE, EN & EE must be correlations!');
end

% number of pairings
npairs=numel(nn);
if(~isequal(npairs,numel(ne),numel(en),numel(ee)))
    error('seizmo:tdsshorzxc:badInput',...
        'NN, NE, EN & EE datasets are not the same size!');
end

% check inputs
if(~isreal(smax) || ~isscalar(smax) || smax<=0)
    error('seizmo:tdsshorzxc:badInput',...
        'SMAX must be a positive real scalar in sec/deg!');
elseif(~any(numel(spts)==[1 2]) || any(fix(spts)~=spts) || any(spts<=2))
    error('seizmo:tdsshorzxc:badInput',...
        'SPTS must be a positive scalar integer >2!');
end

% parse options
pv=parse_tdsshorzxc_pv_pairs(varargin{:});

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
    error('seizmo:tdsshorzxc:badInput',...
        'NN, NE, EN & EE are not synced!');
end
clear b2 b3 b4 delta2 delta3 delta4 npts2 npts3 npts4;
clear st2 st3 st4 ev2 ev3 ev4 stnm2 stnm3 stnm4 evnm2 evnm3 evnm4;

% require north & east orientation or pairwise radial & transverse
if(any([scmp(:,1);scmp2(:,1);scmp3(:,1);scmp4(:,1)]~=90) ...
        || any([ecmp(:,1);ecmp2(:,1);ecmp3(:,1);ecmp4(:,1)]~=90))
    error('seizmo:tdsshorzxc:badInput',...
        'All datasets must have CMPINC==90!');
end
try
    isne=true;
    if(any([scmp(:,2);ecmp(:,2)]~=0))
        error('seizmo:tdsshorzxc:badInput',...
            'NN must have CMPAZ==0 & USER3==0!');
    elseif(any(scmp2(:,2)~=90) || any(ecmp2(:,2)~=0))
        error('seizmo:tdsshorzxc:badInput',...
            'NE must have CMPAZ==90 & USER3==0!');
    elseif(any(scmp3(:,2)~=0) || any(ecmp3(:,2)~=90))
        error('seizmo:tdsshorzxc:badInput',...
            'EN must have CMPAZ==0 & USER3==90!');
    elseif(any([scmp4(:,2);ecmp4(:,2)]~=90))
        error('seizmo:tdsshorzxc:badInput',...
            'EE must have CMPAZ==90 & USER3==90!');
    end
catch
    isne=false; % RT?
    if(any(abs(azdiff(scmp(:,2),gcp))>1) ...
            || any(abs(azdiff(ecmp(:,2),az))>1))
        error('seizmo:tdsshorzxc:badInput',...
            'RR must have CMPAZ==GCP & USER3==AZ!');
    elseif(any(abs(azdiff(scmp2(:,2),gcp+90))>1) ...
            || any(abs(azdiff(ecmp2(:,2),az))>1))
        error('seizmo:tdsshorzxc:badInput',...
            'RT must have CMPAZ==GCP+90 & USER3==AZ!');
    elseif(any(abs(azdiff(scmp3(:,2),gcp))>1) ...
            || any(abs(azdiff(ecmp3(:,2),az+90))>1))
        error('seizmo:tdsshorzxc:badInput',...
            'TR must have CMPAZ==GCP & USER3==AZ+90!');
    elseif(any(abs(azdiff(scmp4(:,2),gcp+90))>1) ...
            || any(abs(azdiff(ecmp4(:,2),az+90))>1))
        error('seizmo:tdsshorzxc:badInput',...
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

% longest record
maxnpts=max(npts);

% get frequencies
nspts=2^nextpow2(maxnpts);
f=(0:nspts/2)/(delta(1)*nspts);  % only +freq
nfreq=numel(f);

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
    error('seizmo:tdsshorzxc:badInput',...
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
r.nsta=nrecs;
r.st=st;
r.butc=autc;
r.eutc=futc;
r.delta=delta(1);
r.npts=maxnpts;
r.polar=pv.polar;
r.x=sx;
r.y=sy;
r.freq=f;
r.method=pv.method;
r.npairs=npairs;
r.center=[clat clon];
r.whiten=false;
r.weights=pv.w;
r.spectra=zeros(spts,'single');
r.orient='radial';
t=r;
t.orient='transverse';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of defaults, checks, & data setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% radius constants
d2r=pi/180;
d2km=6371*d2r;

% get number of slowness points
nslow=prod(spts);

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

% convert phase delay in seconds to sample number
% - note: -p b/c r==s means plane wave is coming from that direction
x=round((-p-b(:,ones(1,nslow))')/delta(1));

% detail message
if(verbose)
    disp('Getting time-domain slowness spectra of horizontals');
    print_time_left(0,npairs);
end

% loop over each correlation
for a=1:npairs
    ok=x(:,a)>0 & x(:,a)<=npts(a);
    r.spectra(ok)=r.spectra(ok)+...
        (um(ok,a).*us(ok,a).*nn(a).dep(x(ok,a)) ...
        +um(ok,a).*vs(ok,a).*ne(a).dep(x(ok,a)) ...
        +vm(ok,a).*us(ok,a).*en(a).dep(x(ok,a)) ...
        +vm(ok,a).*vs(ok,a).*ee(a).dep(x(ok,a)))*pv.w(a);
    t.spectra(ok)=t.spectra(ok)+...
        (vm(ok,a).*vs(ok,a).*nn(a).dep(x(ok,a)) ...
        -vm(ok,a).*us(ok,a).*ne(a).dep(x(ok,a)) ...
        -um(ok,a).*vs(ok,a).*en(a).dep(x(ok,a)) ...
        +um(ok,a).*us(ok,a).*ee(a).dep(x(ok,a)))*pv.w(a);
    if(verbose); print_time_left(a,npairs); end
end

% scale by nfreq
r.spectra=r.spectra/nfreq;
t.spectra=t.spectra/nfreq;

end


function [pv]=parse_tdsshorzxc_pv_pairs(varargin)

% defaults
pv.polar=false;
pv.method='tdxc';
pv.w=[];

% require pv pairs
if(mod(nargin,2))
    error('seizmo:tdsshorzxc:badInput',...
        'Unpaired parameter/value!');
end

% parse pv pairs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:tdsshorzxc:badInput',...
        'Parameters must be specified using a string!');
end
for i=1:2:nargin
    switch varargin{i}
        case {'polar' 'plr' 'pol' 'p'}
            pv.polar=varargin{i+1};
        case {'weights' 'w' 'wgt' 'wgts' 'weight'}
            pv.w=varargin{i+1};
        otherwise
            error('seizmo:tdsshorzxc:badInput',...
                'Unknown parameter: %s !',varargin{i});
    end
end

% check values
if(~isscalar(pv.polar) || (~islogical(pv.polar) && ~isnumeric(pv.polar)))
    error('seizmo:tdsshorzxc:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif(~isnumeric(pv.w) || any(pv.w(:)<0))
    error('seizmo:tdsshorzxc:badInput',...
        'WEIGHTS must be positive!');
end

end


function [lgc]=isxc(data)
[m,s]=getheader(data,'kuser0','kuser1');
if(all(strcmp(m,'MASTER') & strcmp(s,'SLAVE'))); lgc=true;
else lgc=false;
end
end

