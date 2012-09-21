function [s]=geoarf(stlalo,lalo,slow,lalo0,slow0,f0,method,w)
%GEOARF    Calculate the geofss array response function for an array
%
%    Usage:    s=geoarf(stlalo,latlon,slow,latlon0,slow0,f0)
%              s=geoarf(stlalo,latlon,slow,latlon0,slow0,f0,method)
%              s=geoarf(stlalo,latlon,slow,latlon0,slow0,f0,method,weights)
%
%    Description:
%     S=GEOARF(STLALO,LATLON,SLOW,lalo0,slow0,f0) computes the array
%     response function (ARF) for an array at the locations in STLALO with
%     spherical waves passing through the array defined by positions
%     latlon0, horizontal slownesses slow0, and frequencies f0.  The ARF is
%     specifically computed at the positions LATLON (this is the grid of
%     latitude and longitudes that shows the aliasing as a function of
%     position) & horizontal slownesses SLOW at the given f0.  Latitudes
%     & longitudes are expected in degrees with the position arrays
%     organized as [LAT LON].  Slownesses are the magnitude of the
%     horizontal slowness in sec/deg.  The ARF is created using the
%     'center' method as this is faster (see next Usage form to use
%     another method).  The output struct S has the same fields as GEOFSS
%     output but also includes a S.source field which has the following
%     layout:
%      .source.nsrc   - number of sources
%      .source.latlon - latitude/longitude positions as [LAT LON] (deg)
%      .source.slow   - horizontal slowness (sec/deg)
%      .source.freq   - frequency (Hz)
%
%     S=GEOARF(STLALO,LATLON,SLOW,latlon0,slow0,f0,METHOD) defines the
%     beamforming method.  METHOD may be 'center', 'coarray', 'full', or
%     [LAT LON].  The default is 'center' which is extremely fast for large
%     arrays compared to the 'coarray' method as it only pairs each record
%     against the array center (found using ARRAYCENTER) to compute the
%     real-valued spectrum.  The 'coarray' method utilizes information from
%     all unique pairings of records to compute the complex slowness
%     spectrum while the 'full' method uses every possible pairing to do
%     the same.  The 'full' method is significantly slower and gives
%     degraded results compared to the 'coarray' method and so is not
%     recommended except in verification.  The 'center' method gives
%     results that are the same as the 'full' method but does it far
%     faster.  Using [LAT LON] for method is algorithmically the same as
%     the 'center' method but uses the defined coordinates as the center
%     for the array.
%
%     S=GEOARF(STLALO,LATLON,SLOW,latlon0,slow0,f0,METHOD,WEIGHTS)
%     specifies weights for each station in the array for use in the
%     beamforming.  The weights are normalized internally to sum to 1.
%
%    Notes:
%
%    Examples:
%     % Global ARF for a random array:
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     st=randlatlon(100);
%     s=geoarf(st,[lat(:) lon(:)],25,[0 0],25,1/500);
%     plotgeoarf(s);
%
%     % Multi-source ARF for a global array:
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     st=randlatlon(100);
%     s=geoarf(st,[lat(:) lon(:)],25,randlatlon(5),25,1/500);
%     plotgeoarf(geofssavg(s));
%
%    See also: PLOTGEOARF, GEOFSS, GEOFSSXC, GEOFSSHORZ, GEOFSSHORZXC,
%              GEOFSSAVG, GEOFSSSUB, PLOTGEOFSS, GEOFSSFRAMESLIDE,
%              GEOFSSFREQSLIDE, GEOFSSSLOWSLIDE, GEOFSSCORR, GEOFSSDBINFO

%     Version History:
%        July  7, 2010 - initial version
%        Nov. 18, 2010 - add weighting
%        June 12, 2012 - adapt from geofkarf, alter output format, default
%                        to center method, always output struct, don't
%                        average the sources
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2012 at 17:00 GMT

% todo:

% check nargin
error(nargchk(6,8,nargin));

% defaults
if(nargin<7 || isempty(method)); method='center'; end
if(nargin<8); w=[]; end

% valid method strings
valid.METHOD={'center' 'coarray' 'full'};

% check inputs
if(~isreal(stlalo) || ndims(stlalo)~=2 || size(stlalo,2)~=2)
    error('seizmo:geoarf:badInput',...
        'STLALO must be a Nx2 real array of [STLA STLO] in deg!');
elseif(~isreal(lalo) || ndims(lalo)~=2 ...
        || size(lalo,2)~=2 || size(lalo,1)<1)
    error('seizmo:geoarf:badInput',...
        'LATLON must be a Nx2 real array of [LAT LON] in deg!');
elseif(~isreal(slow) || any(slow<=0) || numel(slow)<1)
    error('seizmo:geoarf:badInput',...
        'SLOW must be positive real vector in s/deg!');
elseif(~isreal(slow0) || any(slow0<=0) || numel(slow0)<1)
    error('seizmo:geoarf:badInput',...
        'slow0 must be positive slowness in s/deg!');
elseif(~isreal(lalo0) || ndims(lalo0)~=2 ...
        || size(lalo0,2)~=2 || size(lalo0,1)<1)
    error('seizmo:geoarf:badInput',...
        'latlon0 must be a Nx2 real array of [LAT LON] in deg!');
elseif(~isreal(f0) || any(f0<=0) || numel(f0)<1)
    error('seizmo:geoarf:badInput',...
        'f0 must be positive frequency in Hz!');
elseif(~isequalsizeorscalar(slow0(:),lalo0(:,1),f0(:)))
    error('seizmo:geoarf:badInput',...
        'slow0(:), latlon0(:,1), f0(:) must be equal sized or scalar!');
elseif((isnumeric(method) && (~isreal(method) || ~numel(method)==2)) ...
        || (ischar(method) && ~any(strcmpi(method,valid.METHOD))))
    error('seizmo:geoarf:badInput',...
        'METHOD must be [LAT LON], ''CENTER'', ''COARRAY'' or ''FULL''');
elseif(~isempty(w) && (any(w(:)<0) || ~isreal(w) || sum(w(:))==0))
    error('seizmo:geoarf:badInput',...
        'WEIGHTS must be positive real values!');
end

% number of stations
nsta=size(stlalo,1);

% require 2+ stations
if(nsta<2)
    error('seizmo:geoarf:arrayTooSmall',...
            'GEOARF requires 2+ station locations!');
end

% default weighting
if(isempty(w)); w=ones(nsta,1); end

% convert weights to column vector and normalize
w=w(:);
w=w./sum(w);

% expand spherical wave details
nsrc=max([numel(slow0) size(lalo0,1) numel(f0)]);
if(isscalar(slow0)); slow0(1:nsrc,1)=slow0; end
if(isscalar(f0)); f0(1:nsrc,1)=f0; end
if(size(lalo0,1)==1); lalo0=lalo0(ones(nsrc,1),:); end
slow0=slow0(:); % force column vector
f0=f0(:); % force column vector

% fix method/center
if(ischar(method))
    method=lower(method);
    [clat,clon]=arraycenter(stlalo(:,1),stlalo(:,2));
    switch method
        case 'coarray'
            npairs=nsta*(nsta-1)/2;
        case 'full'
            npairs=nsta*nsta;
        case 'center'
            npairs=nsta;
    end
else
    clat=method(1);
    clon=method(2);
    method='user';
    npairs=nsta;
end

% fix lat/lon
[stlalo(:,1),stlalo(:,2)]=fixlatlon(stlalo(:,1),stlalo(:,2));
[lalo0(:,1),lalo0(:,2)]=fixlatlon(lalo0(:,1),lalo0(:,2));
[lalo(:,1),lalo(:,2)]=fixlatlon(lalo(:,1),lalo(:,2));
nll=size(lalo,1);

% column vector slownesses
slow=slow(:);
nslow=numel(slow);

% create geoARF struct
s.nsta=nsta;
s.stla=stlalo(:,1);
s.stlo=stlalo(:,2);
s.stel=zeros(nsta,1);
s.stdp=zeros(nsta,1);
s.butc=[0 0 0 0 0];
s.eutc=[0 0 0 0 0];
s.delta=nan;
s.npts=0;
s.vector=[nsrc~=1 nslow~=1];
s.latlon=lalo;
s.slow=slow;
s.freq=f0;
s.method=method;
s.npairs=npairs;
s.center=[clat clon];
s.weights=w;
s.source.nsrc=nsrc;
s.source.latlon=lalo0;
s.source.slow=slow0;
s.source.freq=f0;
s.spectra=zeros(nll,nslow,nsrc,'single');

% - get distance between each station and the emitters
% - then get distance difference for each pair for each emitter to get the
%   phasors that steer the array (frequency domain beamforming)
% - also expand weights if necessary
% dd is NLLxNPAIRS
stlalo=stlalo.';
dist=sphericalinv(lalo(:,ones(nsta,1)),lalo(:,2*ones(nsta,1)),...
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
        [row,col]=find(triu(true(nsta),1));
        w=w(col).*w(row);
        w=w./sum(w);
        dd=dist(:,row)-dist(:,col);
    case 'full'
        % retain full coarray
        [row,col]=find(true(nsta));
        w=w(col).*w(row);
        w=w./sum(w);
        dd=dist(:,row)-dist(:,col);
    case {'center' 'user'}
        % relative to array center
        cdist=sphericalinv(lalo(:,1),lalo(:,2),clat,clon);
        dd=cdist(:,ones(nsta,1))-dist;
end
dd=2*pi*1i*dd;
clear dist cdist

% detail message
verbose=seizmoverbose;
if(verbose)
    fprintf('Getting geoARF for %d spherical wave sources\n',nsrc);
    print_time_left(0,nsrc);
end

% loop over spherical wave emitters
for a=1:nsrc
    % spherical wave phasor offsets
    % dd0 is a 1xNPAIRS vector
    d0=sphericalinv(lalo0(a,ones(nsta,1)),lalo0(a,2*ones(nsta,1)),...
        stlalo(1,:),stlalo(2,:));
    switch method
        case 'coarray'
            [row,col]=find(triu(true(nsta),1));
            dd0=d0(row)-d0(col);
        case 'full'
            [row,col]=find(true(nsta));
            dd0=d0(row)-d0(col);
        case {'center' 'user'}
            cd0=sphericalinv(lalo0(a,1),lalo0(a,2),clat,clon);
            dd0=cd0(ones(1,nsta))-d0;
    end
    p0=2*pi*1i*dd0*slow0(a);
    
    % loop over slownesses
    for b=1:nslow
        % beamforming method
        switch method
            case {'full' 'coarray'}
                s.spectra(:,b,a)=...
                    +exp(f0(a)*(slow(b)*dd-p0(ones(nll,1),:)))*w;
            otherwise
                s.spectra(:,b,a)=...
                    +abs(exp(f0(a)*(slow(b)*dd-p0(ones(nll,1),:)))*w).^2;
        end
    end
    if(verbose); print_time_left(a,nsrc); end
end

end
