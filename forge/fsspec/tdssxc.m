function [s]=tdssxc(xcdata,smax,spts,varargin)
%TDSSXC    Estimate slowness spectrum with cross-correlations
%
%    Usage:    s=tdssxc(xcdata,smax,spts)
%              s=tdssxc(...,'polar',true|false,...)
%              s=tdssxc(...,'weights',w,...)
%
%    Description:
%     S=TDSSXC(XCDATA,SMAX,SPTS) computes an estimate of the slowness
%     spectra for an array by time-domain beamforming the correlation
%     dataset XCDATA in a cartesian grid.  The dataset XCDATA is a SEIZMO
%     struct containing array info and correlations.  This function differs
%     from GEOTDSSXC in that the delays are based on plane waves traveling
%     on a planar surface rather than surface waves expanding and
%     contracting on a sphere.  The range of the horizontal slowness grid
%     is given by SMAX (sec/deg) and extends from -SMAX to SMAX for both
%     East/West and North/South directions.  SPTS controls the number of
%     slowness points for both directions (SPTSxSPTS grid).  The output S
%     is a struct containing relevant info and the slowness spectra (with
%     size SPTSxSPTSx1).  The struct layout is:
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
%
%     S=TDSSXC(...,'POLAR',TRUE|FALSE,...) specifies if the spectra is
%     sampled regularly in cartesian or polar coordinates.  Polar coords
%     are useful for slicing the spectra by azimuth (pie slice) or slowness
%     (rings).  Cartesian coords (the default) samples the slowness space
%     regularly in the East/West & North/South directions and so exhibits
%     less distortion in plots of the slowness space. If POLAR=TRUE, SPTS
%     may be given as [SPTS BAZPTS] to control the azimuthal resolution
%     (default is BAZPTS=181 points).
%
%     S=TDSSXC(...,'WEIGHTS',W,...) specifies the relative weights for each
%     correlogram in XCDATA (must match size of XCDATA).
%
%    Notes:
%     - Correlations in XCDATA must have equal and regular sample spacing.
%     - Correlations are expected to have station naming and location info
%       as stored in the header by the function CORRELATE.  This means the
%       headers of correlations done outside of SEIZMO will probably
%       require adjustment.
%     - Attenuation is ignored.
%     - Brzak et al 2009 method uses a Hilbert transform (a 90deg phase
%       shift) that we do not include as it gives incorrect results.  Using
%       a Hilbert transform to get the envelope of the correlations (eg.
%       Burtin et al 2010, Gu & Shen 2012) gives an estimate of the group
%       velocity but at a far lower resolution.
%     - References:
%        Brzak et al 2009, Migration imaging and forward modeling of
%         microseismic noise sources near southern Italy, G3, Vol. 10,
%         No. 1, Q01012, doi:10.1029/2008GC002234
%        Burtin et al 2010, Location of river-induced seismic signal from
%         noise correlation functions, GJI, Vol. 182, pp. 1161-1173,
%         doi:10.1111/j.1365-246X.2010.04701.x
%        Gu & Shen 2012, Microseismic Noise from Large Ice-Covered Lakes?,
%         BSSA, Vol. 102, No. 3, pp. 1155-1166, doi:101785/0120100010
%
%    Examples:
%     % Show slowness spectra for a dataset:
%     plotfss(tdssxc(xcdata,50,201))
%
%     % Filter to 20-25s and compare to FSSXC output:
%     xcdata2=iirfilter(xcdata,'bp','b','o',4,'c',[1/25 1/20],'p',2);
%     plotfss(tdssxc(xcdata2,50,201));
%     plotfss(fssavg(fssxc(xcdata,50,201,[1/25 1/20])))
%
%    See also: TDSSHORZXC, GEOTDSSXC, GEOTDSSHORZXC, FSSXC, GEOFSSXC,
%              ARF, ARFHORZ, GEOARF, GEOARFHORZ

%     Version History:
%        Sep. 28, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 28, 2012 at 14:05 GMT

% todo:

% check nargin
error(nargchk(3,inf,nargin));

% check struct
error(seizmocheck(xcdata,'dep'));

% require xc dataset
if(~isxc(xcdata))
    error('seizmo:tdssxc:badInput',...
        'XCDATA must be correlations with metadata as from CORRELATE!');
end

% number of pairings
npairs=numel(xcdata);

% check inputs
if(~isreal(smax) || ~isscalar(smax) || smax<=0)
    error('seizmo:tdssxc:badInput',...
        'SMAX must be a positive real scalar in sec/deg!');
elseif(~any(numel(spts)==[1 2]) || any(fix(spts)~=spts) || any(spts<=2))
    error('seizmo:tdssxc:badInput',...
        'SPTS must be a positive scalar integer >2!');
end

% parse options
pv=parse_tdssxc_pv_pairs(varargin{:});

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
    error('seizmo:fssxc:badInput',...
        '# of WEIGHTS must match the # of stations or records in XCDATA!');
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
s.nsta=nrecs;
s.st=st;
s.butc=autc;
s.eutc=futc;
s.delta=delta(1);
s.npts=maxnpts;
s.polar=pv.polar;
s.x=sx;
s.y=sy;
s.freq=f;
s.method=pv.method;
s.npairs=npairs;
s.center=[clat clon];
s.whiten=false;
s.weights=pv.w;
s.spectra=zeros(spts,'single');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of defaults, checks, & xcdata setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% radius constants
d2r=pi/180;
d2km=6371*d2r;

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
% spatial difference vectors r (called the coarray) and
% is the phase delay in seconds
%
% p is NSLOWxNPAIRS
if(pv.polar)
    sx=sx(ones(spts(1),1),:); % baz in degrees
    sy=sy(:,ones(spts(2),1))/d2km; % smag in sec/km
    [sx,sy]=deal(sy.*sind(sx),sy.*cosd(sx));
    p=[sx(:) sy(:)]*r;
else % cartesian
    sx=sx/d2km;
    sx=sx(ones(spts(1),1),:);
    sy=fliplr(sx)';
    p=[sx(:) sy(:)]*r;
end
clear r sx sy

% convert phase delay in seconds to sample number
% - note: -p b/c r==s means plane wave is coming from that direction
x=round((-p-b(:,ones(1,prod(spts)))')/delta(1));

% hilbert testing (Brzak et al 2009)
%xcdata=hilbrt(xcdata);

% envelope testing (Burtin et al 2010, Gu & Shen 2012)
%xcdata=envelope(xcdata);

% detail message
if(verbose)
    disp('Getting time-domain slowness spectra');
    print_time_left(0,npairs);
end

% loop over each correlation
for a=1:npairs
    ok=x(:,a)>0 & x(:,a)<=npts(a);
    s.spectra(ok)=s.spectra(ok)+xcdata(a).dep(x(ok,a))*pv.w(a);
    if(verbose); print_time_left(a,npairs); end
end

% scale by nfreq
s.spectra=s.spectra/nfreq;

end


function [pv]=parse_tdssxc_pv_pairs(varargin)

% defaults
pv.polar=false;
pv.method='tdxc';
pv.w=[];

% require pv pairs
if(mod(nargin,2))
    error('seizmo:tdssxc:badInput',...
        'Unpaired parameter/value!');
end

% parse pv pairs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:tdssxc:badInput',...
        'Parameters must be specified using a string!');
end
for i=1:2:nargin
    switch varargin{i}
        case {'polar' 'plr' 'pol' 'p'}
            pv.polar=varargin{i+1};
        case {'method' 'meth' 'm'}
            pv.method=varargin{i+1};
        case {'weights' 'w' 'wgt' 'wgts' 'weight'}
            pv.w=varargin{i+1};
        otherwise
            error('seizmo:tdssxc:badInput',...
                'Unknown parameter: %s !',varargin{i});
    end
end

% valid method strings
valid.METHOD={'tdxc'};

% check values
if(~isscalar(pv.polar) || (~islogical(pv.polar) && ~isnumeric(pv.polar)))
    error('seizmo:fssxc:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif((isnumeric(pv.method) ...
        && (~isreal(pv.method) || ~numel(pv.method)==2)) ...
        || (ischar(pv.method) && ~any(strcmpi(pv.method,valid.METHOD))))
    error('seizmo:tdssxc:badInput',...
        ['METHOD must be one of the following:\n' ...
        '''TDXC'' !']);
elseif(~isnumeric(pv.w) || any(pv.w(:)<0))
    error('seizmo:tdssxc:badInput',...
        'WEIGHTS must be positive!');
end

end


function [lgc]=isxc(data)
[m,s]=getheader(data,'kuser0','kuser1');
if(all(strcmp(m,'MASTER') & strcmp(s,'SLAVE'))); lgc=true;
else lgc=false;
end
end

