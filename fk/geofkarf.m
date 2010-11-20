function [varargout]=geofkarf(stlalo,lalo,s,lalo0,s0,f0,method,w)
%GEOFKARF    Returns the geofk array response function for a seismic array
%
%    Usage:    geofkarf(stlalo,latlon,s,latlon0,s0)
%              geofkarf(stlalo,latlon,s,latlon0,s0,f0)
%              geofkarf(stlalo,latlon,s,latlon0,s0,f0,method)
%              geofkarf(stlalo,latlon,s,latlon0,s0,f0,method,weights)
%              arf=geofkarf(...)
%
%    Description:
%     GEOFKARF(STLALO,LATLON,S,lalo0,s0) computes the array
%     response function (ARF) for an array at the locations in STLALO with
%     spherical waves passing through the array defined by latlon0 & s0.
%     The ARF is specifically computed at the positions and slownesses
%     given by LATLON & S.  Latitudes & longitudes are expected in degrees
%     with the position arrays organized as [LAT LON].  Slownesses are
%     horizontal slowness magnitudes given in sec/deg.  All spherical waves
%     are assumed to have a frequency of 1Hz.  The ARF is plotted in units
%     of decibels.
%
%     GEOFKARF(STLALO,LATLON,S,latlon0,s0,f0) alters the frequency of the
%     spherical wave(s) to f0.  f0 is expected to be in Hz.  The default
%     value is 1 Hz.
%
%     GEOFKARF(STLALO,LATLON,S,latlon0,s0,f0,METHOD) defines the
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
%     GEOFKARF(STLALO,LATLON,S,latlon0,s0,f0,METHOD,WEIGHTS) specifies
%     weights for each station pair in the array (this depends on the
%     method) for use in beamforming.  For example, the 'center' method 
%     requires N weights whereas 'coarray' requires (N*N-N)/2 weights.  The
%     weights are normalized internally to sum to 1.  See the Examples
%     section of FKARF for coarray weight indexing.
%
%     ARF=GEOFKARF(...)  returns the array response function in struct ARF
%     without plotting it.  ARF has the following fields:
%      ARF.beam      --  the array beam response function
%      ARF.nsta      --  number of stations
%      ARF.stla      --  station latitudes
%      ARF.stlo      --  station longitudes
%      ARF.latlon    --  latitude/longitude positions (deg)
%      ARF.horzslow  --  horizontal slownesses (sec/deg)
%      ARF.nsw       --  number of spherical waves
%      ARF.latlon0   --  spherical wave latitude/longitude positions (deg)
%      ARF.horzslow0 --  spherical wave slownesses (sec/deg)
%      ARF.freq0     --  spherical wave frequencies (Hz)
%      ARF.npairs    --  number of station pairs
%      ARF.method    --  beamforming method (center, coarray, full, user)
%      ARF.center    --  array center as [LAT LON]
%      ARF.normdb    --  what 0dB actually corresponds to
%      ARF.volume    --  [true false] (slowness volume, not freq volume)
%      ARF.weights   --  weights used in beam response function
%
%    Notes:
%
%    Examples:
%
%    See also: PLOTGEOFKARF, GEOFKXCVOLUME, SNYQUIST, CHKGEOFKARFSTRUCT,
%              GEOFKARF2MAP, GEOFKSUBARF, UPDATEGEOFKARF, GEOFKARFSLOWSLIDE

%     Version History:
%        July  7, 2010 - initial version
%        Nov. 18, 2010 - add weighting
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 18, 2010 at 17:00 GMT

% todo:

% check nargin
error(nargchk(5,8,nargin));

% defaults
if(nargin<6 || isempty(f0)); f0=1; end
if(nargin<7 || isempty(method)); method='coarray'; end
if(nargin<8); w=[]; end

% valid method strings
valid.METHOD={'center' 'coarray' 'full'};

% check inputs
if(~isreal(stlalo) || ndims(stlalo)~=2 || size(stlalo,2)~=2)
    error('seizmo:geofkarf:badInput',...
        'STLALO must be a Nx2 real array of [STLA STLO] in deg!');
elseif(~isreal(lalo) || ndims(lalo)~=2 || size(lalo,2)~=2)
    error('seizmo:geofkarf:badInput',...
        'LATLON must be a Nx2 real array of [LAT LON] in deg!');
elseif(~isreal(s) || any(s<=0))
    error('seizmo:geofkarf:badInput',...
        'HORZSLOW must be positive real vector in s/deg!');
elseif(~isreal(s0) || any(s0<=0))
    error('seizmo:geofkarf:badInput',...
        's0 must be positive slowness in s/deg!');
elseif(~isreal(lalo0) || ndims(lalo0)~=2 || size(lalo0,2)~=2)
    error('seizmo:geofkarf:badInput',...
        'latlon0 must be a Nx2 real array of [LAT LON] in deg!');
elseif(~isreal(f0) || any(f0<=0))
    error('seizmo:geofkarf:badInput',...
        'f0 must be positive frequency in Hz!');
elseif(~isequalsizeorscalar(s0(:),lalo0(:,1),f0(:)))
    error('seizmo:geofkarf:badInput',...
        's0(:), latlon0(:,1), f0(:) must be equal sized or scalar!');
elseif((isnumeric(method) && (~isreal(method) || ~numel(method)==2)) ...
        || (ischar(method) && ~any(strcmpi(method,valid.METHOD))))
    error('seizmo:geofkarf:badInput',...
        'METHOD must be [LAT LON], ''CENTER'', ''COARRAY'' or ''FULL''');
elseif(~isempty(w) && (any(w(:)<0) || ~isreal(w) || sum(w(:))==0))
    error('seizmo:geofkarf:badInput',...
        'WEIGHTS must be positive real values!');
end

% convert weights to column vector and normalize
w=w(:);
w=w./sum(w);

% number of stations
nrecs=size(stlalo,1);

% require 2+ stations
if(nrecs<2)
    error('seizmo:geofkarf:arrayTooSmall',...
            'GEOFKARF requires 2+ station locations!');
end

% expand spherical wave details
nsw=max([numel(s0) size(lalo0,1) numel(f0)]);
if(isscalar(s0)); s0(1:nsw,1)=s0; end
if(isscalar(f0)); f0(1:nsw,1)=f0; end
if(size(lalo0,1)==1); lalo0=lalo0(ones(nsw,1),:); end
s0=s0(:);
f0=f0(:);

% fix method/center
if(ischar(method))
    method=lower(method);
    [clat,clon]=arraycenter(stlalo(:,1),stlalo(:,2));
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

% check that number of weights is equal to number of pairs
if(isempty(w)); w=ones(npairs,1)/npairs; end
if(numel(w)~=npairs)
    error('seizmo:geofkarf:badInput',...
        ['WEIGHTS must have ' num2str(npairs) ...
        ' elements for method: ' method '!']);
end

% fix lat/lon
[stlalo(:,1),stlalo(:,2)]=fixlatlon(stlalo(:,1),stlalo(:,2));
[lalo0(:,1),lalo0(:,2)]=fixlatlon(lalo0(:,1),lalo0(:,2));
[lalo(:,1),lalo(:,2)]=fixlatlon(lalo(:,1),lalo(:,2));
nll=size(lalo,1);

% column vector slownesses
s=s(:);
nslow=numel(s);

% create geo arf struct
arf.nsta=nrecs;
arf.stla=stlalo(:,1);
arf.stlo=stlalo(:,2);
arf.latlon=lalo;
arf.horzslow=s;
arf.nsw=nsw;
arf.latlon0=lalo0;
arf.horzslow0=s0;
arf.freq0=f0;
arf.npairs=npairs;
arf.method=method;
arf.center=[clat clon];
arf.weights=w;
arf.volume=[true false]; % summed across freq, but not slowness
arf.beam=zeros(nll,nslow,'single');

% - get distance between each station and the emitters
% - then get distance difference for each pair for each emitter to get the
%   phasors that steer the array (beamforming)
% dd is NLLxNPAIRS
stlalo=stlalo.';
dist=sphericalinv(lalo(:,ones(nrecs,1)),lalo(:,2*ones(nrecs,1)),...
    stlalo(ones(nll,1),:),stlalo(2*ones(nll,1),:));
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
        [row,col]=find(triu(true(nrecs),1));
        dd=dist(:,row)-dist(:,col);
    case 'full'
        % retain full coarray
        [row,col]=find(true(nrecs));
        dd=dist(:,row)-dist(:,col);
    case {'center' 'user'}
        % relative to array center
        cdist=sphericalinv(lalo(:,1),lalo(:,2),clat,clon);
        dd=cdist(:,ones(nrecs,1))-dist;
end
dd=2*pi*1i*dd;
clear dist cdist

% detail message
verbose=seizmoverbose;
if(verbose)
    fprintf('Getting geofk ARF for spherical wave source:\n');
end

% loop over spherical wave emitters
for a=1:nsw
    % detail message
    if(verbose)
        fprintf(['SLOWNESS: %gs/deg, LAT: %gdeg, ' ...
            'LON: %gdeg, FREQ: %gHz\n'],...
            s0(a),lalo0(a,1),lalo0(a,2),f0(a));
    end
    
    % spherical wave phasor offsets
    % dd0 is a 1xNPAIRS vector
    d0=sphericalinv(lalo0(a,ones(nrecs,1)),lalo0(a,2*ones(nrecs,1)),...
        stlalo(1,:),stlalo(2,:));
    switch method
        case 'coarray'
            [row,col]=find(triu(true(nrecs),1));
            dd0=d0(row)-d0(col);
        case 'full'
            [row,col]=find(true(nrecs));
            dd0=d0(row)-d0(col);
        case {'center' 'user'}
            cd0=sphericalinv(lalo0(a,1),lalo0(a,2),clat,clon);
            dd0=cd0(ones(1,nrecs))-d0;
    end
    p0=2*pi*1i*dd0*s0(a);
    
    % loop over slownesses
    for b=1:nslow
        % get beam
        switch method
            case {'full' 'coarray'}
                arf.beam(:,b)=arf.beam(:,b)...
                    +exp(f0(a)*(s(b)*dd-p0(ones(nll,1),:)))*w;
            otherwise
                arf.beam(:,b)=arf.beam(:,b)...
                    +abs(exp(f0(a)*(s(b)*dd-p0(ones(nll,1),:)))*w).^2;
        end
    end
end

% convert to dB
switch method
    case {'full' 'coarray'}
        % using full and real here gives the exact plots of
        % Koper, Seats, and Benz 2010 in BSSA
        arf.beam=10*log10(abs(real(arf.beam))/nsw);
    otherwise
        arf.beam=10*log10(arf.beam/nsw);
end

% normalize so max peak is at 0dB
arf.normdb=max(arf.beam(:));
arf.beam=arf.beam-arf.normdb;

% return if output
if(nargout)
    varargout{1}=arf;
    return;
else
    plotgeofkarf(geofkarf2map(arf));
end

end
