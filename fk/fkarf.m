function [varargout]=fkarf(stla,stlo,smax,spts,s,baz,f,arf)
%FKARF    Returns the fk array response function for a seismic array
%
%    Usage:    fkarf(stla,stlo,smax)
%              fkarf(stla,stlo,smax,spts)
%              fkarf(stla,stlo,smax,spts,s,baz)
%              fkarf(stla,stlo,smax,spts,s,baz,f)
%              arf=fkarf(...)
%              fkarf(stla,stlo,smax,spts,s,baz,f,arf)
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
%     S=0 & BAZ=0.
%
%     FKARF(STLA,STLO,SMAX,SPTS,S,BAZ,F) alters the frequency of the plane
%     wave to F.  F is expected to be in Hz.  The default value is 1.
%
%     ARF=FKARF(...) returns the array response function in struct ARF
%     without plotting it.  This is useful for stacking ARFs for a multiple
%     plane wave response.  ARF has the following fields:
%      ARF.response  --  the array response
%      ARF.nsta      --  number of stations
%      ARF.stla      --  station latitudes
%      ARF.stlo      --  station longitudes
%      ARF.smax      --  max slowness in the array response
%      ARF.npw       --  number of plane waves
%      ARF.s         --  plane wave slownesses
%      ARF.baz       --  plane wave backazimuths
%      ARF.f         --  plane wave frequencies
%
%     FKARF(STLA,STLO,SMAX,SPTS,S,BAZ,F,ARF) combines a previous array
%     response function ARF with the current run's.  See the previous usage
%     form for the structure of ARF.  Note that STLO, STLA, SMAX & SPTS
%     should be equal to that in ARF.
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
%      arf=fkarf(stla,stlo,50,201,20,0,0.03);
%      arf=fkarf(stla,stlo,50,201,10,0,0.03);
%      fkarf(stla,stlo,50,201,20,45,0.03);
%
%    See also: FKMAP, SNYQUIST, PLOTFKMAP, KXY2SLOWBAZ, SLOWBAZ2KXY

%     Version History:
%        May   1, 2010 - initial version
%        May   3, 2010 - show slowness nyquist ring, doc update
%        May   4, 2010 - circle is its own function now
%        May   7, 2010 - only doing one triangle gives better response and
%                        takes less than half the time
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   7, 2010 at 02:30 GMT

% todo:
% - ability to define an array center vs coarray
% - full coarray vs 1 triangle

% check nargin
msg=nargchk(3,8,nargin);
if(~isempty(msg)); error(msg); end

% constants
d2r=pi/180;
d2km=6371*d2r;

% defaults
if(nargin<4 || isempty(spts)); spts=101; end
if(nargin<5 || isempty(s)); s=0; end
if(nargin<6 || isempty(baz)); baz=0; end
if(nargin<7 || isempty(f)); f=1; end

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
elseif(~isreal(s) || ~isscalar(s) || s<0)
    error('seizmo:fkarf:badInput',...
        'S must be a positive real scalar in s/deg!');
elseif(~isreal(baz) || ~isscalar(baz))
    error('seizmo:fkarf:badInput',...
        'BAZ must be a real scalar in degrees!');
elseif(~isreal(f) || ~isscalar(f) || f<=0)
    error('seizmo:fkarf:badInput',...
        'F must be a positive real scalar in Hz!');
end

% number of stations
nsta=max(numel(stla),numel(stlo));
if(isscalar(stla)); stla(1:nsta,1)=stla; end
if(isscalar(stlo)); stlo(1:nsta,1)=stlo; end

% handle arf
arffields={'response' 'nsta' 'stla' 'stlo' 'smax' 'npw' 's' 'baz' 'f'};
if(nargin<8 || isempty(arf))
    % create arf
    arf.response=zeros(spts);
    arf.nsta=nsta;
    arf.stla=stla;
    arf.stlo=stlo;
    arf.smax=smax;
    arf.npw=1;
    arf.s=s;
    arf.baz=baz;
    arf.f=f;
else
    % check old arf
    if(~isstruct(arf) || any(~ismember(arffields,fieldnames(arf))))
        error('seizmo:fkarf:badInput',...
            'ARF is not a proper struct!');
    elseif(arf.nsta~=nsta || ~isequal(arf.stla,stla) ...
            || ~isequal(arf.stlo,stlo))
        warning('seizmo:fkarf:badInput',...
            'Station list has changed!');
    elseif(~isequal(size(arf.response),[spts spts]) ...
            || arf.smax~=smax || arf.npw~=fix(arf.npw) || arf.npw<0 ...
            || isinf(arf.npw) || isnan(arf.npw) ...
            || numel(arf.s)~=arf.npw || numel(arf.baz)~=arf.npw ...
            || numel(arf.f)~=arf.npw)
        error('seizmo:fkarf:badInput',...
            'ARF does not match current inputs or is corrupt!');
    end
    arf.response=10.^(arf.response/10)*arf.npw*nsta^2;
    arf.npw=arf.npw+1;
    arf.s=[arf.s(:); s];
    arf.baz=[arf.baz(:); baz];
    arf.f=[arf.f(:); f];
end

% calculate coarray (relative positions)
% r=(x  ,y  )
%     ij  ij
%
% x is km east
% y is km north
%
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
[dist,az]=vincentyinv(stla(:,ones(nsta,1)),stlo(:,ones(nsta,1)),...
    stla(:,ones(nsta,1))',stlo(:,ones(nsta,1))');
az=az*d2r;
x=dist.*sin(az);
y=dist.*cos(az);

% break down plane wave slowness/baz to sx/sy
s=s/d2km; % s/km
baz=baz*d2r;
sx0=s*sin(baz);
sy0=s*cos(baz);

% make wavenumber arrays
smax=smax/d2km;
sx=-smax:2*smax/(spts-1):smax;
sy=fliplr(sx)';
kk0x=2*pi*f*1i*(sx(ones(spts,1),:)-sx0);
kk0y=2*pi*f*1i*(sy(:,ones(spts,1))-sy0);

% get indices to go through
% (upper/lower/both triangles give the same result)
idx=find(triu(true(nsta),1))';   % upper triangle
%idx=find(tril(true(nsta),-1))'; % lower triangle
%idx=find(~eye(nsta))';          % both triangles (no diagonal)
%idx=find(true(nsta))';          % everything
nidx=numel(idx);

% detail message
verbose=seizmoverbose;
if(verbose)
    disp('Generating Array Response Function');
    print_time_left(0,nidx);
end

% get array response
%  2*pi*i*(k-k0)r
% e
for i=1:nidx
    arf.response=arf.response+exp(kk0x*x(idx(i))+kk0y*y(idx(i)));
    if(verbose); print_time_left(i,nidx); end
end

% convert to dB
% - note the use of abs to handle slightly negative terms
arf.response=10*log10(abs(abs(arf.response))/(arf.npw*nidx));

% return if output
if(nargout); varargout{1}=arf; return; end

% get nyquist slowness (plot in a bit)
[clat,clon]=arraycenter(stla,stlo);
[e,n]=geographic2enu(stla,stlo,0,clat,clon,0);
tri=delaunay(e,n);
friends=[tri(:,1:2); tri(:,2:3); tri(:,[3 1])];
friends=unique([min(friends,[],2) max(friends,[],2)],'rows');
dist=vincentyinv(stla(friends(:,1)),stlo(friends(:,1)),...
                 stla(friends(:,2)),stlo(friends(:,2)));
snyq=snyquist(min(dist),f); % closest 2 stations

% plotting slowness space

% first plot the slowness map
figure;
imagesc(sx*d2km,fliplr(sx*d2km),arf.response);
set(gca,'ydir','normal');
set(gca,'clim',[-12 0]);
hold on

% next plot the bull's eye
% Phase:       Rg    Lg    Sn    Pn    Sdiff  Pdiff  PKPcdiff
% Vel (km/s):  3.0   4.0   4.5   8.15  13.3   25.1   54.0
% S (s/deg):   37.0  27.8  24.7  13.6  8.36   4.43   2.06
smax=smax*d2km;
if(smax>=37)
    ph=[37 27.8 24.7 13.6 8.36 4.43];
elseif(smax>=28)
    ph=[27.8 24.7 13.6 8.36 4.43];
elseif(smax>=25)
    ph=[24.7 13.6 8.36 4.43];
elseif(smax>=14)
    ph=[13.6 8.36 4.43 2.06];
elseif(smax>=8.5)
    ph=[8.36 4.43 2.06];
elseif(smax>=4.5)
    ph=[4.43 2.06];
else
    ph=2.06;
end
[x,y]=circle(ph(1),12);
[x2,y2]=circle(ph(end),12);
plot([x; x2],[y; y2],'w','linewidth',1,'tag','bullseye');
for i=ph
    [x,y]=circle(i);
    plot(x,y,'w','linewidth',1,'tag','bullseye');
end

% last plot the nyquist ring about the plane wave location
[x,y]=circle(snyq);
x=x+sx0*d2km; y=y+sy0*d2km;
plot(x,y,'r','linewidth',2,'tag','nyquist ring');
hold off

% finally take care of coloring/labels/etc
title(['Array Response Function @ ' ...
    num2str(f) 'Hz (' num2str(1/f) 's)'],'fontweight','bold');
xlabel('East/West Slowness (s/deg)','fontweight','bold');
ylabel('North/South Slowness (s/deg)','fontweight','bold');
colormap(ritz);
c=colorbar('eastoutside','fontweight','bold');
set(c,'xaxislocation','top');
set(gca,'fontweight','bold');
xlabel(c,'dB')
axis equal
xlim([-smax smax])
ylim([-smax smax])

end
