function [varargout]=fkarf(stla,stlo,smax,spts,s0,baz0,f0,polar,center)
%FKARF    Returns the fk array response function for a seismic array
%
%    Usage:    fkarf(stla,stlo,smax)
%              fkarf(stla,stlo,smax,spts)
%              fkarf(stla,stlo,smax,spts,s,baz)
%              fkarf(stla,stlo,smax,spts,s,baz,f)
%              fkarf(stla,stlo,smax,spts,s,baz,f,polar)
%              fkarf(stla,stlo,smax,spts,s,baz,f,polar,center)
%              arf=fkarf(...)
%
%    Description: FKARF(STLA,STLO,SMAX) shows the array response function
%     (ARF) for a 1Hz plane wave incoming from directly below an array of
%     sensors at geographic positions given by latitudes STLA and
%     longitudes STLO.  The ARF is shown in East/West and North/South
%     slowness space from -SMAX to SMAX where SMAX is expected to be in
%     seconds per degree.  It is sampled 101 times in both the East/West
%     and North/South directions giving a 101x101 image.  ARF is plotted
%     in units of decibels.
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
%     FKARF(STLA,STLO,SMAX,SPTS,S,BAZ,F,POLAR,CENTER) defines the array
%     center.  CENTER may be [LAT LON], 'center', 'coarray', or 'full'.
%     The default is 'coarray'.  The 'center' option finds the center
%     position of the array by averaging the station positions (using
%     ARRAYCENTER).  Both 'coarray' and 'full' are essentially centerless
%     methods using the relative positioning between every possible pairing
%     of stations in the array.  The 'full' method includes redundant and
%     same station pairings (and will always give poorer results compared
%     to 'coarray').
%
%     ARF=FKARF(...) returns the array response function in struct ARF
%     without plotting it.  This is useful for stacking ARFs for a multiple
%     plane wave response.  ARF has the following fields:
%      ARF.response  --  the array response
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
%      ARF.center    --  array center or method
%      ARF.normdb    --  what 0dB actually corresponds to
%
%    Notes:
%     - The circles of the bull's eye in the plot correspond to several
%       surface/diffracted/head seismic phases (which appear depends on the
%       plot's maximum slowness):
%        Phase:       Rg    Lg    Sn    Pn    Sdiff  Pdiff  PKPcdiff
%        Vel (km/s):  3.0   4.0   4.5   8.15  13.3   25.1   54.0
%        S (s/deg):   37.0  27.8  24.7  13.6  8.36   4.43   2.06
%       The radial lines correspond to 30deg steps in backazimuth.
%     - The red circle around the expected location of the plane wave on
%       the slowness diagram indicates the range of slownesses within the
%       nyquist slowness from that location.
%
%    Examples:
%     Show the array response function for +/-5s/deg for a dataset:
%      [stla,stlo]=getheader(data,'stla','stlo');
%      fkarf(stla,stlo,5);
%
%     Get a multi-plane wave response:
%      arf=fkarf(stla,stlo,50,201,[20 10 20],[0 0 45],0.03);
%
%    See also: FKMAP, PLOTFKMAP, FKVOLUME, PLAYFKVOLUME, FK4D
%              KXY2SLOWBAZ, SLOWBAZ2KXY, SNYQUIST

%     Version History:
%        May   1, 2010 - initial version
%        May   3, 2010 - show slowness nyquist ring, doc update
%        May   4, 2010 - circle is its own function now
%        May   7, 2010 - only doing one triangle gives better response and
%                        takes less than half the time
%        May  10, 2010 - added in options available to FKMAP
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  10, 2010 at 17:30 GMT

% todo:

% check nargin
msg=nargchk(3,10,nargin);
if(~isempty(msg)); error(msg); end

% constants
d2r=pi/180;
d2km=6371*d2r;

% defaults
if(nargin<4 || isempty(spts)); spts=101; end
if(nargin<5 || isempty(s0)); s0=0; end
if(nargin<6 || isempty(baz0)); baz0=0; end
if(nargin<7 || isempty(f0)); f0=1; end
if(nargin<8 || isempty(polar)); polar=false; end
if(nargin<9 || isempty(center)); center='coarray'; end

% valid center strings
valid.CENTER={'center' 'coarray' 'full'};

% check inputs
if(~isreal(stla) || ~isreal(stlo) || ~isequalsizeorscalar(stla,stlo))
    error('seizmo:fkarf:badInput',...
        'STLA & STLO must be equal sized real arrays in degrees!');
elseif(~isreal(smax) || ~isscalar(smax) || smax<=0)
    error('seizmo:fkarf:badInput',...
        'SMAX must be a positive real scalar in s/deg!');
elseif(~isscalar(spts) || fix(spts)~=spts || spts<=0)
    error('seizmo:fkarf:badInput',...
        'SPTS must be a positive scalar integer!');
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
elseif((isnumeric(center) && (~isreal(center) || ~numel(center)==2)) ...
        || (ischar(center) && ~any(strcmpi(center,valid.CENTER))))
    error('seizmo:fkarf:badInput',...
        'CENTER must be [LAT LON], ''CENTER'', ''COARRAY'' or ''FULL''');
end

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
arf.center=center;

% fix center
if(ischar(center))
    center=lower(center);
else
    clat=center(1);
    clon=center(2);
    center='user';
end

% get relative positions from center
% r=(x  ,y  )
%     ij  ij
%
% x is km east
% y is km north
%
% r is 2xNR
switch center
    case 'coarray'
        % centerless (make coarray)
        % [ r   r   ... r
        %    11  12      1N
        %   r   r   ... r
        %    21  22      2N
        %    .   .  .    .
        %    .   .   .   .
        %    .   .    .  .
        %   r   r   ... r   ]
        %    N1  N2      NN
        [dist,az]=vincentyinv(...
            stla(:,ones(nrecs,1))',stlo(:,ones(nrecs,1))',...
            stla(:,ones(nrecs,1)),stlo(:,ones(nrecs,1)));
        idx=triu(true(nrecs),1);
        dist=dist(idx);
        az=az(idx);
    case 'full'
        % centerless too
        [dist,az]=vincentyinv(...
            stla(:,ones(nrecs,1))',stlo(:,ones(nrecs,1))',...
            stla(:,ones(nrecs,1)),stlo(:,ones(nrecs,1)));
    case 'center'
        % get array center
        [clat,clon]=arraycenter(stla,stlo);
        [dist,az]=vincentyinv(clat,clon,stla,stlo);
    otherwise % user
        % array center was specified
        [dist,az]=vincentyinv(clat,clon,stla,stlo);
end
az=az*d2r;
r=[dist(:).*sin(az(:)) dist(:).*cos(az(:))]';
nidx=size(r,2);
clear dist az

% make projection arrays
% p=2*pi*i*s*r
%
% where s is the slowness vector s=(s ,s ) and is NSx2
%                                    x  y
% p is NSxNR
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
    arf.response=zeros(spts,bazpts);
else % cartesian
    spts=spts(1);
    sx=-smax:2*smax/(spts-1):smax;
    arf.x=sx*d2km;
    arf.y=fliplr(sx*d2km)';
    sx=sx(ones(spts,1),:);
    sy=fliplr(sx)';
    p=2*pi*1i*[sx(:) sy(:)]*r;
    clear sx sy
    arf.response=zeros(spts);
end
ns=size(p,1);

% offset projection by plane wave locations
s=2*pi*1i*[s0(:).*sin(d2r*baz0(:)) s0(:).*cos(d2r*baz0(:))]/d2km;

% detail message
verbose=seizmoverbose;
if(verbose)
    fprintf('Getting fk Map for plane wave at:\n');
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

    % get response
    switch center
        case {'full' 'coarray'}
            arf.response(:)=arf.response(:)...
                +sum(exp(f0(a)*(p-p0)),2);
        otherwise
            arf.response(:)=arf.response(:)...
                +10*log10(abs(sum(exp(f0(a)*(p-p0)),2)).^2/nidx);
    end
end

% convert to dB
switch center
    case {'full' 'coarray'}
        % using full and real here gives the exact plots of
        % Koper, Seats, and Benz 2010 in BSSA
        arf.response=...
            10*log10(abs(real(arf.response))/(npw*nidx));
    otherwise
        arf.response=arf.response/npw;
end

% normalize so max peak is at 0dB
arf.normdb=max(arf.response(:));
arf.response=arf.response-arf.normdb;

% return if output
if(nargout)
    varargout{1}=arf;
    return;
else
    figure;
    plotfkarf(arf);
end

end
