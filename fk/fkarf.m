function [varargout]=fkarf(stla,stlo,smax,spts,s0,baz0,f0,polar,method,w)
%FKARF    Returns the fk array response function for a seismic array
%
%    Usage:    fkarf(stla,stlo,smax)
%              fkarf(stla,stlo,smax,spts)
%              fkarf(stla,stlo,smax,spts,s,baz)
%              fkarf(stla,stlo,smax,spts,s,baz,f)
%              fkarf(stla,stlo,smax,spts,s,baz,f,polar)
%              fkarf(stla,stlo,smax,spts,s,baz,f,polar,method)
%              fkarf(stla,stlo,smax,spts,s,baz,f,polar,method,weights)
%              arf=fkarf(...)
%
%    Description:
%     FKARF(STLA,STLO,SMAX) shows the array response function (ARF) for a
%     1Hz plane wave incoming from directly below an array of sensors at
%     geographic positions given by latitudes STLA and longitudes STLO.
%     The ARF is shown in East/West and North/South slowness space from
%     -SMAX to SMAX where SMAX is expected to be in seconds per degree.  It
%     is sampled 101 times in both the East/West and North/South directions
%     giving a 101x101 image.  ARF is plotted in units of decibels.
%
%     FKARF(STLA,STLO,SMAX,SPTS) changes the number of samples in slowness
%     space to SPTSxSPTS.  Higher values will provide better resolution at
%     the cost of more computation (ie time).  The default value is 101.
%
%     FKARF(STLA,STLO,SMAX,SPTS,S,BAZ) plots the ARF for a plane wave given
%     by slowness S and backazimuth BAZ.  S is expected in seconds per
%     degree.  BAZ is in degrees clockwise from North and indicates the
%     direction from the array to the wave (180 degree from the direction
%     the plane wave is heading through the array).  The default values are
%     S=0 & BAZ=0.  You may specify multiple plane waves by giving vectors
%     of values for S & BAZ.
%
%     FKARF(STLA,STLO,SMAX,SPTS,S,BAZ,F) alters the frequency of the plane
%     wave to F.  F is expected to be in Hz.  The default value is 1.
%
%     FKARF(STLA,STLO,SMAX,SPTS,S,BAZ,F,POLAR) decides if the slowness map
%     is sampled regularly in cartesian or polar coordinates.  Polar coords
%     are useful for slicing the volume by azimuth (pie slice) or slowness
%     (rings).  Cartesian coords (the default) samples the slowness space
%     regularly in the East/West & North/South directions and so exhibits
%     less distortion of the slowness space.
%
%     FKARF(STLA,STLO,SMAX,SPTS,S,BAZ,F,POLAR,METHOD) defines the
%     beamforming method.  METHOD may be 'center', 'coarray', 'full', or
%     [LAT LON].  The default is 'coarray' which utilizes information from
%     all unique record pairings in the beamforming and is the default.
%     The 'full' method will utilize all possible pairings including
%     pairing records with themselves and pairing records as (1st, 2nd) &
%     (2nd, 1st) making this method quite redundant and slow.  The 'center'
%     option only pairs each record against the array center (found using
%     ARRAYCENTER) and is extremely fast for large arrays compared to the
%     'coarray' & 'full' methods.  Both 'center' and 'full' methods give
%     slightly degraded results compared to 'coarray'.  Using [LAT LON] for
%     method is essentially the same as the 'center' method but uses the
%     defined coordinates as the center for the array.
%
%     FKARF(STLA,STLO,SMAX,SPTS,S,BAZ,F,POLAR,METHOD,WEIGHTS) specifies
%     weights for each station pair in the array (this depends on the
%     method) for use in beamforming.  For example, the 'center' method 
%     requires N weights whereas 'coarray' requires (N*N-N)/2 weights.  The
%     weights are normalized internally to sum to 1.  See the Examples
%     section for coarray weight indexing.
%
%     ARF=FKARF(...) returns the array response function in struct ARF
%     without plotting it.  ARF has the following fields:
%      ARF.beam      --  the array beam response function
%      ARF.nsta      --  number of stations
%      ARF.stla      --  station latitudes
%      ARF.stlo      --  station longitudes
%      ARF.x         --  east/west slowness or azimuth values
%      ARF.y         --  north/south or radial slowness values
%      ARF.npw       --  number of plane waves
%      ARF.s         --  plane wave slownesses
%      ARF.baz       --  plane wave backazimuths
%      ARF.f         --  plane wave frequencies
%      ARF.polar     --  true if slowness is sampled in polar coordinates
%      ARF.npairs    --  number of pairs
%      ARF.method    --  beamforming method (center, coarray, full, user)
%      ARF.center    --  array center as [LAT LON]
%      ARF.normdb    --  what 0dB actually corresponds to
%      ARF.weights   --  weights used in beam response function
%
%    Notes:
%
%    Examples:
%     % Show the array response function for +/-5s/deg for a dataset:
%     [stla,stlo]=getheader(data,'stla','stlo');
%     fkarf(stla,stlo,5);
%
%     % Get a multi-plane wave response:
%     arf=fkarf(stla,stlo,50,201,[20 10 20],[0 0 45],0.03);
%
%     % Indices (into station weights) for coarray weights:
%     [i,j]=find(triu(true(nsta),1));
%     w=ones(nsta,1); % start with uniform weighting
%     w(2:2:end)=2;   % even index stations are weighted up
%     w=w(i).*w(j);   % get weights for coarray
%
%    See also: PLOTFKARF, FKMAP, KXY2SLOWBAZ, SLOWBAZ2KXY, SNYQUIST

%     Version History:
%        May   1, 2010 - initial version
%        May   3, 2010 - show slowness nyquist ring, doc update
%        May   4, 2010 - circle is its own function now
%        May   7, 2010 - only doing one triangle gives better beam and
%                        takes less than half the time
%        May  10, 2010 - added in options available to FKMAP
%        May  18, 2010 - minor doc touch
%        June 16, 2010 - fixed nargchk, improved see also section and notes
%        July  6, 2010 - major update to struct, doc update, high latitude
%                        fix, arf data is now single precision
%        July  7, 2010 - center/user method for multi-plane wave ARF now is
%                        fixed
%        Nov. 16, 2010 - added weighting (franklin witnessed)
%        Nov. 17, 2010 - forgot .weights field
%        Nov. 18, 2010 - coarray weight indexing example, fixed scaling bug
%        Apr. 27, 2012 - spts now allows input as [spts bazpts]
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 27, 2012 at 14:40 GMT

% todo:

% check nargin
error(nargchk(3,10,nargin));

% constants
d2r=pi/180;
d2km=6371*d2r;

% defaults
if(nargin<4 || isempty(spts)); spts=101; end
if(nargin<5 || isempty(s0)); s0=0; end
if(nargin<6 || isempty(baz0)); baz0=0; end
if(nargin<7 || isempty(f0)); f0=1; end
if(nargin<8 || isempty(polar)); polar=false; end
if(nargin<9 || isempty(method)); method='coarray'; end
if(nargin<10); w=[]; end

% valid method strings
valid.METHOD={'center' 'coarray' 'full'};

% check inputs
if(~isreal(stla) || ~isreal(stlo) || ~isequalsizeorscalar(stla,stlo))
    error('seizmo:fkarf:badInput',...
        'STLA & STLO must be equal sized real arrays in degrees!');
elseif(~isreal(smax) || ~isscalar(smax) || smax<=0)
    error('seizmo:fkarf:badInput',...
        'SMAX must be a positive real scalar in s/deg!');
elseif(~any(numel(spts)==[1 2]) || any(fix(spts)~=spts) || any(spts<=2))
    error('seizmo:fkarf:badInput',...
        'SPTS must be a positive scalar integer >2!');
elseif(~isreal(s0) || any(s0<0))
    error('seizmo:fkarf:badInput',...
        'S must be a positive slowness in s/deg!');
elseif(~isreal(baz0))
    error('seizmo:fkarf:badInput',...
        'BAZ must be in degrees!');
elseif(~isreal(f0) || any(f0<=0))
    error('seizmo:fkarf:badInput',...
        'F must be positive frequency in Hz!');
elseif(~isequalsizeorscalar(s0,baz0,f0))
    error('seizmo:fkarf:badInput',...
        'S, BAZ, F must be equal sized or scalar!');
elseif(~isscalar(polar) || (~islogical(polar) && ~isnumeric(polar)))
    error('seizmo:fkarf:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif((isnumeric(method) && (~isreal(method) || ~numel(method)==2)) ...
        || (ischar(method) && ~any(strcmpi(method,valid.METHOD))))
    error('seizmo:fkarf:badInput',...
        'METHOD must be [LAT LON], ''CENTER'', ''COARRAY'' or ''FULL''');
elseif(~isempty(w) && (any(w(:)<0) || ~isreal(w) || sum(w(:))==0))
    error('seizmo:fkarf:badInput',...
        'WEIGHTS must be positive real values!');
end

% convert weights to column vector and normalize
w=w(:);
w=w./sum(w);

% number of stations
nrecs=max(numel(stla),numel(stlo));
if(isscalar(stla)); stla(1:nrecs,1)=stla; end
if(isscalar(stlo)); stlo(1:nrecs,1)=stlo; end

% require 2+ stations
if(nrecs<2)
    error('seizmo:fkarf:arrayTooSmall',...
            'FKARF requires 2+ station locations!');
end

% expand plane wave details
npw=max([numel(s0) numel(baz0) numel(f0)]);
if(isscalar(s0)); s0(1:npw,1)=s0; end
if(isscalar(baz0)); baz0(1:npw,1)=baz0; end
if(isscalar(f0)); f0(1:npw,1)=f0; end
s0=s0(:);
baz0=baz0(:);
f0=f0(:);

% create arf
arf.nsta=nrecs;
arf.stla=stla;
arf.stlo=stlo;
arf.npw=npw;
arf.s=s0;
arf.baz=baz0;
arf.f=f0;
arf.polar=polar;

% fix method/center
if(ischar(method))
    method=lower(method);
    [clat,clon]=arraycenter(stla,stlo);
    switch method
        case 'coarray'
            npairs=nrecs*(nrecs-1)/2;
        case 'full'
            npairs=nrecs*nrecs;
        case 'center'
            npairs=nrecs;
    end
else
    clat=method(1);
    clon=method(2);
    method='user';
    npairs=nrecs;
end
arf.npairs=npairs;
arf.method=method;
arf.center=[clat clon];

% check that number of weights is equal to number of pairs
if(isempty(w)); w=ones(npairs,1)/npairs; end
if(numel(w)~=npairs)
    error('seizmo:fkarf:badInput',...
        ['WEIGHTS must have ' num2str(npairs) ...
        ' elements for method: ' method '!']);
end
arf.weights=w;

% get relative positions for each pair
% r=(x  ,y  )
%     ij  ij
%
% x is km east
% y is km north
%
% r is a 2xNPAIRS matrix
switch method
    case 'coarray'
        % [ r   r   ... r
        %    11  12      1N
        %   r   r   ... r
        %    21  22      2N
        %    .   .  .    .
        %    .   .   .   .
        %    .   .    .  .
        %   r   r   ... r   ]
        %    N1  N2      NN
        %
        % then we just use the
        % upper triangle of that
        [e,n]=geographic2enu(stla,stlo,0,clat,clon,0);
        e=e(:,ones(nrecs,1))'-e(:,ones(nrecs,1));
        n=n(:,ones(nrecs,1))'-n(:,ones(nrecs,1));
        idx=triu(true(nrecs),1);
        e=e(idx);
        n=n(idx);
    case 'full'
        % retain full coarray
        [e,n]=geographic2enu(stla,stlo,0,clat,clon,0);
        e=e(:,ones(nrecs,1))'-e(:,ones(nrecs,1));
        n=n(:,ones(nrecs,1))'-n(:,ones(nrecs,1));
    case 'center'
        % each record relative to array center
        [e,n]=geographic2enu(stla,stlo,0,clat,clon,0);
    otherwise % user
        % each record relative to defined array center
        [e,n]=geographic2enu(stla,stlo,0,clat,clon,0);
end
r=[e(:) n(:)]';
clear e n

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
    arf.y=smag'*d2km;
    smag=smag(ones(bazpts,1),:)';
    baz=(0:bazpts-1)/(bazpts-1)*360*d2r;
    arf.x=baz/d2r;
    baz=baz(ones(spts,1),:);
    p=2*pi*1i*[smag(:).*sin(baz(:)) smag(:).*cos(baz(:))]*r;
    clear smag baz
    arf.beam=zeros(spts,bazpts,'single');
else % cartesian
    spts=spts(1);
    sx=-smax:2*smax/(spts-1):smax;
    arf.x=sx*d2km;
    arf.y=fliplr(sx*d2km)';
    sx=sx(ones(spts,1),:);
    sy=fliplr(sx)';
    p=2*pi*1i*[sx(:) sy(:)]*r;
    clear sx sy
    arf.beam=zeros(spts,'single');
end
ns=size(p,1);

% slowness offsets based on plane wave locations
s=2*pi*1i*[s0(:).*sin(d2r*baz0(:)) s0(:).*cos(d2r*baz0(:))]/d2km;

% detail message
verbose=seizmoverbose;
if(verbose)
    fprintf('Getting fk ARF for plane wave source:\n');
end

% loop over plane waves
for a=1:npw
    % detail message
    if(verbose)
        fprintf('SLOWNESS: %gs/deg, BAZ: %gdeg, FREQ: %gHz\n',...
            s0(a),baz0(a),f0(a));
    end
    
    % projection of plane wave
    p0=s(a,:)*r;
    p0=p0(ones(ns,1),:);

    % get beam
    switch method
        case {'full' 'coarray'}
            arf.beam(:)=arf.beam(:)+exp(f0(a)*(p-p0))*w;
        otherwise
            arf.beam(:)=arf.beam(:)+abs(exp(f0(a)*(p-p0))*w).^2;
    end
end

% convert to dB
switch method
    case {'full' 'coarray'}
        % using full and real here gives the exact plots of
        % Koper, Seats, and Benz 2010 in BSSA
        arf.beam=10*log10(abs(real(arf.beam))/npw);
    otherwise
        arf.beam=10*log10(arf.beam/npw);
end

% normalize so max peak is at 0dB
arf.normdb=max(arf.beam(:));
arf.beam=arf.beam-arf.normdb;

% return if output
if(nargout)
    varargout{1}=arf;
    return;
else
    plotfkarf(arf);
end

end
